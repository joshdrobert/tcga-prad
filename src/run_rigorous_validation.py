from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import proportional_hazard_test
from lifelines.utils import concordance_index
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score, roc_curve
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESEARCH_V2_DIR = PROJECT_ROOT / "research" / "v2"
RANDOM_STATE = 20260626


def build_pipeline(numeric_cols: list[str], categorical_cols: list[str]) -> Pipeline:
    transformers = []
    if numeric_cols:
        transformers.append(
            ("num", Pipeline([("imputer", SimpleImputer(strategy="median")), ("scaler", StandardScaler())]), numeric_cols)
        )
    if categorical_cols:
        transformers.append(
            ("cat", Pipeline([("imputer", SimpleImputer(strategy="most_frequent")), ("onehot", OneHotEncoder(handle_unknown="ignore"))]), categorical_cols)
        )
    preprocessor = ColumnTransformer(transformers=transformers, remainder="drop")
    return Pipeline(
        [
            ("preprocess", preprocessor),
            ("classifier", LogisticRegression(max_iter=2000, random_state=RANDOM_STATE)),
        ]
    )


def compute_bootstrap_intervals(y: np.ndarray, p: np.ndarray, repetitions: int = 1500) -> dict[str, float]:
    rng = np.random.default_rng(RANDOM_STATE)
    base_auroc = float(roc_auc_score(y, p))
    base_auprc = float(average_precision_score(y, p))
    base_brier = float(brier_score_loss(y, p))

    aurocs, auprcs, briers = [], [], []
    n = len(y)
    for _ in range(repetitions):
        idx = rng.integers(0, n, n)
        if len(np.unique(y[idx])) < 2:
            continue
        aurocs.append(roc_auc_score(y[idx], p[idx]))
        auprcs.append(average_precision_score(y[idx], p[idx]))
        briers.append(brier_score_loss(y[idx], p[idx]))

    return {
        "auroc": base_auroc,
        "auroc_ci_low": float(np.quantile(aurocs, 0.025)),
        "auroc_ci_high": float(np.quantile(aurocs, 0.975)),
        "auprc": base_auprc,
        "auprc_ci_low": float(np.quantile(auprcs, 0.025)),
        "auprc_ci_high": float(np.quantile(auprcs, 0.975)),
        "brier": base_brier,
        "brier_ci_low": float(np.quantile(briers, 0.025)),
        "brier_ci_high": float(np.quantile(briers, 0.975)),
    }


def compute_paired_difference(
    y: np.ndarray, p_baseline: np.ndarray, p_augmented: np.ndarray, repetitions: int = 2000
) -> dict[str, float]:
    """Calculate paired difference in AUROC and bootstrap uncertainty on identical participants."""
    rng = np.random.default_rng(RANDOM_STATE)
    delta_auroc_point = float(roc_auc_score(y, p_augmented) - roc_auc_score(y, p_baseline))
    delta_brier_point = float(brier_score_loss(y, p_augmented) - brier_score_loss(y, p_baseline))

    delta_aurocs = []
    delta_briers = []
    n = len(y)
    for _ in range(repetitions):
        idx = rng.integers(0, n, n)
        if len(np.unique(y[idx])) < 2:
            continue
        d_auc = roc_auc_score(y[idx], p_augmented[idx]) - roc_auc_score(y[idx], p_baseline[idx])
        d_bri = brier_score_loss(y[idx], p_augmented[idx]) - brier_score_loss(y[idx], p_baseline[idx])
        delta_aurocs.append(d_auc)
        delta_briers.append(d_bri)

    p_val_auc = float(2 * min(np.mean([d <= 0 for d in delta_aurocs]), np.mean([d >= 0 for d in delta_aurocs])))

    return {
        "delta_auroc": delta_auroc_point,
        "delta_auroc_ci_low": float(np.quantile(delta_aurocs, 0.025)),
        "delta_auroc_ci_high": float(np.quantile(delta_aurocs, 0.975)),
        "delta_auroc_p_value": p_val_auc,
        "delta_brier": delta_brier_point,
        "delta_brier_ci_low": float(np.quantile(delta_briers, 0.025)),
        "delta_brier_ci_high": float(np.quantile(delta_briers, 0.975)),
    }


def compute_calibration_slope(y: np.ndarray, p: np.ndarray) -> tuple[float, float]:
    clipped = np.clip(p, 1e-6, 1.0 - 1e-6)
    logit = np.log(clipped / (1.0 - clipped)).reshape(-1, 1)
    mod = LogisticRegression(C=1e6, solver="lbfgs").fit(logit, y)
    return float(mod.intercept_[0]), float(mod.coef_[0, 0])


def evaluate_feature_set(
    df: pd.DataFrame, features: list[str], outcome_col: str
) -> tuple[dict, np.ndarray, np.ndarray]:
    """Fit 5-fold CV logistic regression for a given feature set and outcome column."""
    eval_df = df[df[outcome_col].notna()].copy()
    X = eval_df[features].copy()
    y = eval_df[outcome_col].astype(int).values
    num_cols = X.select_dtypes(include=["number", "bool"]).columns.tolist()
    cat_cols = [c for c in X.columns if c not in num_cols]

    pipe = build_pipeline(num_cols, cat_cols)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_STATE)
    probs = cross_val_predict(pipe, X, y, cv=cv, method="predict_proba")[:, 1]

    metrics = compute_bootstrap_intervals(y, probs)
    cal_intercept, cal_slope = compute_calibration_slope(y, probs)
    metrics["calibration_intercept"] = cal_intercept
    metrics["calibration_slope"] = cal_slope
    metrics["n"] = len(y)
    metrics["events"] = int(y.sum())
    metrics["non_events"] = int((y == 0).sum())
    metrics["event_rate"] = float(y.mean())
    metrics["predictor_count"] = len(features)
    metrics["outcome_column"] = outcome_col
    return metrics, probs, y


def evaluate_cv_cox_model(
    df: pd.DataFrame, duration_col: str, event_col: str, covariates: list[str]
) -> tuple[float, float, pd.DataFrame, dict[str, float]]:
    """Fit Cox PH model and compute apparent C-index, 5-fold CV C-index, and Schoenfeld test."""
    valid_df = df[[duration_col, event_col] + covariates].dropna().copy()

    # Apparent model
    cph = CoxPHFitter()
    cph.fit(valid_df, duration_col=duration_col, event_col=event_col)
    apparent_c = float(cph.concordance_index_)

    # 5-fold cross-validated out-of-fold risk scores
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_STATE)
    oof_risk = np.zeros(len(valid_df))
    durations = valid_df[duration_col].values
    events = valid_df[event_col].values

    for train_idx, val_idx in cv.split(valid_df, events):
        train_data = valid_df.iloc[train_idx]
        val_data = valid_df.iloc[val_idx]
        fold_cph = CoxPHFitter()
        fold_cph.fit(train_data, duration_col=duration_col, event_col=event_col)
        oof_risk[val_idx] = fold_cph.predict_partial_hazard(val_data).values

    cv_c = float(concordance_index(durations, -oof_risk, events))
    summary_df = cph.summary.reset_index().rename(columns={"index": "covariate"})

    # Proportional Hazards assumption test
    ph_test = proportional_hazard_test(cph, valid_df, time_transform="rank")
    schoenfeld_p = {cov: float(ph_test.p_value[i]) for i, cov in enumerate(covariates)}

    return apparent_c, cv_c, summary_df, schoenfeld_p


def main() -> None:
    RESEARCH_V2_DIR.mkdir(parents=True, exist_ok=True)
    tables_dir = PROJECT_ROOT / "tables"
    figures_dir = PROJECT_ROOT / "figures"
    artifacts_dir = PROJECT_ROOT / "research" / "artifacts"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    artifacts_dir.mkdir(parents=True, exist_ok=True)

    recurrence = pd.read_csv(PROJECT_ROOT / "data_processed" / "tcga_prad_recurrence_cohort.csv")
    expr = pd.read_csv(PROJECT_ROOT / "data_processed" / "tcga_prad_clinical_expression_cohort.csv")
    expr_cols = [c for c in expr.columns if c.startswith("expr_")]
    data = recurrence.merge(expr[["submitter_id"] + expr_cols], on="submitter_id", how="left")

    # Predictor feature definitions:
    # 1. Pure baseline peri-operative staging (Day 0 at radical prostatectomy):
    pure_clinical_features = ["age_at_diagnosis_years", "pathologic_t", "pathologic_n", "gleason_score_numeric"]
    pure_all_features = pure_clinical_features + expr_cols

    # 2. Sensitivity clinical features including post-operative surveillance PSA:
    surveillance_clinical_features = pure_clinical_features + ["psa_value_numeric"]
    surveillance_all_features = surveillance_clinical_features + expr_cols

    # =========================================================================
    # PART 1: PRIMARY DEFINITIVE BINARY PREDICTION ANALYSIS
    # Population: outcome_definitive_recurrence (n=263: 72 positive, 191 double-negative)
    # Excluded: 46 asymmetric surveillance cases and 5 indeterminate cases
    # =========================================================================
    data_definitive = data[data["outcome_definitive_recurrence"].notna()].copy()

    # Paired molecular subpopulation within definitive cohort (n=188: 45 positive, 143 double-negative)
    definitive_complete_mask = data_definitive[expr_cols].notna().all(axis=1)
    data_def_paired = data_definitive[definitive_complete_mask].copy()

    # 1A. Full Definitive Cohort (n=263)
    res_def_clin, prob_def_clin, y_def = evaluate_feature_set(
        data_definitive, pure_clinical_features, outcome_col="outcome_definitive_recurrence"
    )
    res_def_psa, prob_def_psa, _ = evaluate_feature_set(
        data_definitive, surveillance_clinical_features, outcome_col="outcome_definitive_recurrence"
    )

    # 1B. Paired Definitive Molecular Subset (n=188, exactly matched participants)
    res_def_paired_clin, p_def_p_clin, y_def_p = evaluate_feature_set(
        data_def_paired, pure_clinical_features, outcome_col="outcome_definitive_recurrence"
    )
    res_def_paired_mol, p_def_p_mol, _ = evaluate_feature_set(
        data_def_paired, pure_all_features, outcome_col="outcome_definitive_recurrence"
    )
    paired_diff_pure = compute_paired_difference(y_def_p, p_def_p_clin, p_def_p_mol)

    res_def_paired_psa_clin, p_def_p_psa_clin, _ = evaluate_feature_set(
        data_def_paired, surveillance_clinical_features, outcome_col="outcome_definitive_recurrence"
    )
    res_def_paired_psa_mol, p_def_p_psa_mol, _ = evaluate_feature_set(
        data_def_paired, surveillance_all_features, outcome_col="outcome_definitive_recurrence"
    )
    paired_diff_surv = compute_paired_difference(y_def_p, p_def_p_psa_clin, p_def_p_psa_mol)

    # =========================================================================
    # PART 2: SENSITIVITY BINARY ANALYSIS (EXPLORATORY COMPOSITE)
    # Population: outcome_exploratory_composite (n=309: 72 positive, 237 non-events)
    # Admitting 46 asymmetric surveillance cases; excluding only 5 indeterminate cases
    # =========================================================================
    data_exploratory = data[data["outcome_exploratory_composite"].notna()].copy()
    exp_complete_mask = data_exploratory[expr_cols].notna().all(axis=1)
    data_exp_paired = data_exploratory[exp_complete_mask].copy()

    res_exp_clin, prob_exp_clin, y_exp = evaluate_feature_set(
        data_exploratory, pure_clinical_features, outcome_col="outcome_exploratory_composite"
    )
    res_exp_paired_clin, p_exp_p_clin, y_exp_p = evaluate_feature_set(
        data_exp_paired, pure_clinical_features, outcome_col="outcome_exploratory_composite"
    )
    res_exp_paired_mol, p_exp_p_mol, _ = evaluate_feature_set(
        data_exp_paired, pure_all_features, outcome_col="outcome_exploratory_composite"
    )
    paired_diff_exp = compute_paired_difference(y_exp_p, p_exp_p_clin, p_exp_p_mol)

    # =========================================================================
    # PART 3: HISTORICAL BENCHMARK (LEGACY UNADJUDICATED N=314)
    # Retained strictly for historical reproducibility; uncertain cases naively mapped to 0
    # =========================================================================
    complete_mask_all = data[expr_cols].notna().all(axis=1)
    data_complete_all = data[complete_mask_all].copy()

    res_hist_full_clin, p_hist_full_clin, _ = evaluate_feature_set(
        data, surveillance_clinical_features, outcome_col="outcome_recurrence_or_progression"
    )
    res_hist_full_mol, p_hist_full_mol, _ = evaluate_feature_set(
        data, surveillance_all_features, outcome_col="outcome_recurrence_or_progression"
    )
    res_hist_paired_clin, p_hist_p_clin, y_hist_p = evaluate_feature_set(
        data_complete_all, surveillance_clinical_features, outcome_col="outcome_recurrence_or_progression"
    )
    res_hist_paired_mol, p_hist_p_mol, _ = evaluate_feature_set(
        data_complete_all, surveillance_all_features, outcome_col="outcome_recurrence_or_progression"
    )
    paired_diff_hist = compute_paired_difference(y_hist_p, p_hist_p_clin, p_hist_p_mol)

    # Compile Comprehensive Models Comparison Table
    model_rows = [
        # Primary Definitive (n=263 and paired n=188)
        {"analysis_tier": "PRIMARY_DEFINITIVE", "cohort": "Primary Definitive (n=263)", "model": "Pure Baseline Clinical (Age, T, N, Gleason)", **res_def_clin},
        {"analysis_tier": "PRIMARY_DEFINITIVE", "cohort": "Primary Definitive (n=263)", "model": "Baseline Clinical + Surveillance PSA (Sensitivity)", **res_def_psa},
        {"analysis_tier": "PRIMARY_DEFINITIVE", "cohort": "Paired Definitive Molecular (n=188)", "model": "Pure Baseline Clinical (Age, T, N, Gleason)", **res_def_paired_clin},
        {"analysis_tier": "PRIMARY_DEFINITIVE", "cohort": "Paired Definitive Molecular (n=188)", "model": "Pure Baseline Clinical + 13-Gene Expression Panel", **res_def_paired_mol},
        {"analysis_tier": "PRIMARY_DEFINITIVE", "cohort": "Paired Definitive Molecular (n=188)", "model": "Clinical + Surveillance PSA + 13-Gene Panel (Sensitivity)", **res_def_paired_psa_mol},
        # Exploratory Sensitivity (n=309 and paired n=223)
        {"analysis_tier": "SENSITIVITY_EXPLORATORY", "cohort": "Exploratory Composite (n=309)", "model": "Pure Baseline Clinical (Age, T, N, Gleason)", **res_exp_clin},
        {"analysis_tier": "SENSITIVITY_EXPLORATORY", "cohort": "Paired Exploratory Molecular (n=223)", "model": "Pure Baseline Clinical (Age, T, N, Gleason)", **res_exp_paired_clin},
        {"analysis_tier": "SENSITIVITY_EXPLORATORY", "cohort": "Paired Exploratory Molecular (n=223)", "model": "Pure Baseline Clinical + 13-Gene Expression Panel", **res_exp_paired_mol},
        # Historical Benchmark (n=314 and paired n=228)
        {"analysis_tier": "HISTORICAL_BENCHMARK", "cohort": "Historical Unadjudicated (n=314)", "model": "Clinical Only (Surveillance PSA included)", **res_hist_full_clin},
        {"analysis_tier": "HISTORICAL_BENCHMARK", "cohort": "Historical Unadjudicated (n=314)", "model": "Clinical + 13-Gene Panel Imputed", **res_hist_full_mol},
        {"analysis_tier": "HISTORICAL_BENCHMARK", "cohort": "Historical Paired Complete (n=228)", "model": "Clinical Only (Surveillance PSA included)", **res_hist_paired_clin},
        {"analysis_tier": "HISTORICAL_BENCHMARK", "cohort": "Historical Paired Complete (n=228)", "model": "Clinical + 13-Gene Panel Complete", **res_hist_paired_mol},
    ]
    full_table = pd.DataFrame(model_rows)
    full_table.to_csv(tables_dir / "table_prad_rigorous_validation.csv", index=False)
    full_table.to_csv(RESEARCH_V2_DIR / "MODEL_COMPARISONS.csv", index=False)

    # Compile Primary Paired Complete Case Comparison Table
    paired_molecular_table = pd.DataFrame(
        [
            {
                "analysis_tier": "PRIMARY_DEFINITIVE",
                "cohort": "Paired Definitive Molecular",
                "sample_size_n": len(y_def_p),
                "events_n": int(y_def_p.sum()),
                "non_events_n": int((y_def_p == 0).sum()),
                "baseline_model": "Pure Baseline Clinical (Age, T, N, Gleason)",
                "augmented_model": "Pure Baseline Clinical + 13-Gene Panel",
                "baseline_auroc": res_def_paired_clin["auroc"],
                "augmented_auroc": res_def_paired_mol["auroc"],
                "delta_auroc": paired_diff_pure["delta_auroc"],
                "delta_auroc_95ci_low": paired_diff_pure["delta_auroc_ci_low"],
                "delta_auroc_95ci_high": paired_diff_pure["delta_auroc_ci_high"],
                "delta_auroc_p_value": paired_diff_pure["delta_auroc_p_value"],
                "baseline_brier": res_def_paired_clin["brier"],
                "augmented_brier": res_def_paired_mol["brier"],
                "delta_brier": paired_diff_pure["delta_brier"],
                "delta_brier_95ci_low": paired_diff_pure["delta_brier_ci_low"],
                "delta_brier_95ci_high": paired_diff_pure["delta_brier_ci_high"],
            },
            {
                "analysis_tier": "PRIMARY_DEFINITIVE_SURVEILLANCE_SENSITIVITY",
                "cohort": "Paired Definitive with Surveillance PSA",
                "sample_size_n": len(y_def_p),
                "events_n": int(y_def_p.sum()),
                "non_events_n": int((y_def_p == 0).sum()),
                "baseline_model": "Clinical + Surveillance PSA",
                "augmented_model": "Clinical + Surveillance PSA + 13-Gene Panel",
                "baseline_auroc": res_def_paired_psa_clin["auroc"],
                "augmented_auroc": res_def_paired_psa_mol["auroc"],
                "delta_auroc": paired_diff_surv["delta_auroc"],
                "delta_auroc_95ci_low": paired_diff_surv["delta_auroc_ci_low"],
                "delta_auroc_95ci_high": paired_diff_surv["delta_auroc_ci_high"],
                "delta_auroc_p_value": paired_diff_surv["delta_auroc_p_value"],
                "baseline_brier": res_def_paired_psa_clin["brier"],
                "augmented_brier": res_def_paired_psa_mol["brier"],
                "delta_brier": paired_diff_surv["delta_brier"],
                "delta_brier_95ci_low": paired_diff_surv["delta_brier_ci_low"],
                "delta_brier_95ci_high": paired_diff_surv["delta_brier_ci_high"],
            },
            {
                "analysis_tier": "SENSITIVITY_EXPLORATORY",
                "cohort": "Paired Exploratory Molecular",
                "sample_size_n": len(y_exp_p),
                "events_n": int(y_exp_p.sum()),
                "non_events_n": int((y_exp_p == 0).sum()),
                "baseline_model": "Pure Baseline Clinical (Age, T, N, Gleason)",
                "augmented_model": "Pure Baseline Clinical + 13-Gene Panel",
                "baseline_auroc": res_exp_paired_clin["auroc"],
                "augmented_auroc": res_exp_paired_mol["auroc"],
                "delta_auroc": paired_diff_exp["delta_auroc"],
                "delta_auroc_95ci_low": paired_diff_exp["delta_auroc_ci_low"],
                "delta_auroc_95ci_high": paired_diff_exp["delta_auroc_ci_high"],
                "delta_auroc_p_value": paired_diff_exp["delta_auroc_p_value"],
                "baseline_brier": res_exp_paired_clin["brier"],
                "augmented_brier": res_exp_paired_mol["brier"],
                "delta_brier": paired_diff_exp["delta_brier"],
                "delta_brier_95ci_low": paired_diff_exp["delta_brier_ci_low"],
                "delta_brier_95ci_high": paired_diff_exp["delta_brier_ci_high"],
            },
            {
                "analysis_tier": "HISTORICAL_BENCHMARK",
                "cohort": "Historical Paired Complete",
                "sample_size_n": len(y_hist_p),
                "events_n": int(y_hist_p.sum()),
                "non_events_n": int((y_hist_p == 0).sum()),
                "baseline_model": "Clinical Only (Surveillance PSA)",
                "augmented_model": "Clinical + 13-Gene Panel",
                "baseline_auroc": res_hist_paired_clin["auroc"],
                "augmented_auroc": res_hist_paired_mol["auroc"],
                "delta_auroc": paired_diff_hist["delta_auroc"],
                "delta_auroc_95ci_low": paired_diff_hist["delta_auroc_ci_low"],
                "delta_auroc_95ci_high": paired_diff_hist["delta_auroc_ci_high"],
                "delta_auroc_p_value": paired_diff_hist["delta_auroc_p_value"],
                "baseline_brier": res_hist_paired_clin["brier"],
                "augmented_brier": res_hist_paired_mol["brier"],
                "delta_brier": paired_diff_hist["delta_brier"],
                "delta_brier_95ci_low": paired_diff_hist["delta_brier_ci_low"],
                "delta_brier_95ci_high": paired_diff_hist["delta_brier_ci_high"],
            },
        ]
    )
    paired_molecular_table.to_csv(tables_dir / "table_prad_paired_complete_case_validation.csv", index=False)
    paired_molecular_table.to_csv(RESEARCH_V2_DIR / "PAIRED_MOLECULAR_COMPARISON.csv", index=False)

    # =========================================================================
    # PART 4: SURVIVAL ANALYSIS — STRICT SENSITIVITY VS EXPLORATORY COHORTS
    # =========================================================================
    # Cohort 4A: Strict Survival Sensitivity Cohort (n=257: 66 dated events, 191 double-negative non-events)
    strict_mask = (recurrence["survival_event_indicator"] == 1.0) | (recurrence["composite_status_detail"] == "DOCUMENTED_NON_RECURRENCE")
    strict_surv_df = recurrence[strict_mask & recurrence["survival_duration_days"].notna() & (recurrence["survival_duration_days"] > 0)].copy()

    kmf_strict = KaplanMeierFitter()
    kmf_strict.fit(
        strict_surv_df["survival_duration_days"].values,
        strict_surv_df["survival_event_indicator"].astype(int).values,
        label=f"Strict Sensitivity Cohort (Double-Negatives Only, n={len(strict_surv_df)})",
    )
    strict_r1 = float(kmf_strict.predict(365.25))
    strict_r3 = float(kmf_strict.predict(365.25 * 3))
    strict_r5 = float(kmf_strict.predict(365.25 * 5))
    strict_median = float(kmf_strict.median_survival_time_) if pd.notna(kmf_strict.median_survival_time_) else np.nan

    cox_covs = ["age_at_diagnosis_years", "gleason_score_numeric", "psa_value_numeric"]
    strict_app_c, strict_cv_c, strict_cox_summary, strict_ph_test = evaluate_cv_cox_model(
        strict_surv_df, "survival_duration_days", "survival_event_indicator", cox_covs
    )

    # Cohort 4B: Exploratory Survival Cohort (n=303: 66 dated events, 237 non-events including 46 asymmetric)
    exp_surv_df = recurrence[recurrence["survival_duration_days"].notna() & (recurrence["survival_duration_days"] > 0)].copy()

    kmf_exp = KaplanMeierFitter()
    kmf_exp.fit(
        exp_surv_df["survival_duration_days"].values,
        exp_surv_df["survival_event_indicator"].astype(int).values,
        label=f"Exploratory Survival Cohort (Admitting Asymmetric Surveillance, n={len(exp_surv_df)})",
    )
    exp_r1 = float(kmf_exp.predict(365.25))
    exp_r3 = float(kmf_exp.predict(365.25 * 3))
    exp_r5 = float(kmf_exp.predict(365.25 * 5))
    exp_median = float(kmf_exp.median_survival_time_) if pd.notna(kmf_exp.median_survival_time_) else np.nan

    exp_app_c, exp_cv_c, exp_cox_summary, exp_ph_test = evaluate_cv_cox_model(
        exp_surv_df, "survival_duration_days", "survival_event_indicator", cox_covs
    )

    surv_results_tbl = pd.DataFrame(
        [
            {
                "analysis_definition": "Strict Sensitivity Cohort (Double-Negatives Only)",
                "evaluable_n": len(strict_surv_df),
                "dated_events": int(strict_surv_df["survival_event_indicator"].sum()),
                "censored_non_events": int((strict_surv_df["survival_event_indicator"] == 0).sum()),
                "undated_events_excluded": int(recurrence["undated_recurrence_flag"].sum()),
                "asymmetric_non_events_excluded": 46,
                "indeterminate_records_excluded": 5,
                "rfs_1yr": strict_r1,
                "rfs_3yr": strict_r3,
                "rfs_5yr": strict_r5,
                "median_rfs_days": strict_median,
                "cox_apparent_c_index": strict_app_c,
                "cox_cv_c_index": strict_cv_c,
                "schoenfeld_ph_p_age": strict_ph_test.get("age_at_diagnosis_years", np.nan),
                "schoenfeld_ph_p_gleason": strict_ph_test.get("gleason_score_numeric", np.nan),
                "schoenfeld_ph_p_psa": strict_ph_test.get("psa_value_numeric", np.nan),
            },
            {
                "analysis_definition": "Exploratory Survival Cohort (Admitting Asymmetric Surveillance)",
                "evaluable_n": len(exp_surv_df),
                "dated_events": int(exp_surv_df["survival_event_indicator"].sum()),
                "censored_non_events": int((exp_surv_df["survival_event_indicator"] == 0).sum()),
                "undated_events_excluded": int(recurrence["undated_recurrence_flag"].sum()),
                "asymmetric_non_events_excluded": 0,
                "indeterminate_records_excluded": 5,
                "rfs_1yr": exp_r1,
                "rfs_3yr": exp_r3,
                "rfs_5yr": exp_r5,
                "median_rfs_days": exp_median,
                "cox_apparent_c_index": exp_app_c,
                "cox_cv_c_index": exp_cv_c,
                "schoenfeld_ph_p_age": exp_ph_test.get("age_at_diagnosis_years", np.nan),
                "schoenfeld_ph_p_gleason": exp_ph_test.get("gleason_score_numeric", np.nan),
                "schoenfeld_ph_p_psa": exp_ph_test.get("psa_value_numeric", np.nan),
            },
        ]
    )
    surv_results_tbl.to_csv(tables_dir / "table_prad_survival_analysis.csv", index=False)
    surv_results_tbl.to_csv(RESEARCH_V2_DIR / "SURVIVAL_ANALYSIS.csv", index=False)
    strict_cox_summary.to_csv(tables_dir / "table_prad_cox_summary_strict.csv", index=False)
    exp_cox_summary.to_csv(tables_dir / "table_prad_cox_summary_exploratory.csv", index=False)

    # =========================================================================
    # PART 5: PUBLICATION FIGURES GENERATION
    # =========================================================================
    # Figure 1: Calibration Curves (Primary Paired Definitive vs Exploratory)
    fig, ax = plt.subplots(figsize=(6, 5), dpi=300)
    cal_configs = [
        ("Definitive Paired Clinical (n=188)", p_def_p_clin, y_def_p, "#1e3a8a", "-"),
        ("Definitive Paired Clin+RNA (n=188)", p_def_p_mol, y_def_p, "#059669", "-"),
        ("Definitive Paired Clin+PSA (n=188)", p_def_p_psa_clin, y_def_p, "#d97706", "--"),
        ("Exploratory Paired Clin (n=223)", p_exp_p_clin, y_exp_p, "#6b7280", ":"),
    ]
    for name, p, y, color, style in cal_configs:
        df_cal = pd.DataFrame({"p": p, "y": y})
        df_cal["bin"] = pd.qcut(df_cal["p"], 6, duplicates="drop")
        binned = df_cal.groupby("bin", observed=True).agg(mean_p=("p", "mean"), mean_y=("y", "mean"))
        ax.plot(binned["mean_p"], binned["mean_y"], marker="o", color=color, linestyle=style, label=name, markersize=4)

    ax.plot([0, 1], [0, 1], "k--", alpha=0.6, label="Ideal Calibration")
    ax.set_xlabel("Mean Predicted Recurrence Probability", fontsize=9)
    ax.set_ylabel("Observed Recurrence Rate", fontsize=9)
    ax.set_xlim(0, 0.8)
    ax.set_ylim(0, 0.8)
    ax.legend(fontsize=7, loc="upper left")
    ax.set_title("5-Fold Cross-Validation Calibration (Definitive vs Exploratory)", fontsize=10)
    fig.tight_layout()
    fig.savefig(figures_dir / "prad_calibration.png")
    plt.close(fig)

    # Figure 2: Kaplan-Meier Survival Comparison (Strict Sensitivity vs Exploratory)
    fig, ax = plt.subplots(figsize=(6.5, 4.5), dpi=300)
    kmf_strict.plot_survival_function(ax=ax, color="#1e3a8a", ci_show=True, label=f"Strict Sensitivity (n=257, {int(strict_surv_df['survival_event_indicator'].sum())} events)")
    kmf_exp.plot_survival_function(ax=ax, color="#d97706", linestyle="--", ci_show=True, label=f"Exploratory Cohort (n=303, {int(exp_surv_df['survival_event_indicator'].sum())} events)")
    ax.set_xlabel("Time Post-Radical Prostatectomy (Days)", fontsize=9)
    ax.set_ylabel("Recurrence-Free Survival Probability", fontsize=9)
    ax.set_title("Kaplan-Meier Recurrence-Free Survival (Strict Sensitivity vs Exploratory)", fontsize=10)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7, loc="lower left")
    fig.tight_layout()
    fig.savefig(figures_dir / "prad_km_survival.png")
    plt.close(fig)

    # Figure 3: Paired ROC Curves on the Primary Definitive Molecular Cohort (n=188)
    fig, ax = plt.subplots(figsize=(5.5, 5), dpi=300)
    fpr_clin, tpr_clin, _ = roc_curve(y_def_p, p_def_p_clin)
    fpr_mol, tpr_mol, _ = roc_curve(y_def_p, p_def_p_mol)
    ax.plot(fpr_clin, tpr_clin, color="#1e3a8a", lw=2, label=f"Pure Clinical Baseline (AUROC = {res_def_paired_clin['auroc']:.3f})")
    ax.plot(fpr_mol, tpr_mol, color="#059669", lw=2, linestyle="--", label=f"Clinical + 13-Gene RNA (AUROC = {res_def_paired_mol['auroc']:.3f})")
    ax.plot([0, 1], [0, 1], "k--", alpha=0.5, label="Chance Line (0.50)")
    ax.set_xlabel("1 - Specificity (False Positive Rate)", fontsize=9)
    ax.set_ylabel("Sensitivity (True Positive Rate)", fontsize=9)
    ax.set_title(f"Paired ROC on Primary Definitive Cohort (n=188, 45 Events)\nΔAUROC = {paired_diff_pure['delta_auroc']:+.3f} (95% CI: [{paired_diff_pure['delta_auroc_ci_low']:+.3f}, {paired_diff_pure['delta_auroc_ci_high']:+.3f}], p={paired_diff_pure['delta_auroc_p_value']:.3f})", fontsize=8.5)
    ax.legend(fontsize=7, loc="lower right")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(figures_dir / "prad_paired_roc.png")
    plt.close(fig)

    # =========================================================================
    # PART 6: METRICS ARTIFACT COMPILATION & DISK EXPORT
    # =========================================================================
    cbio_tbl = pd.read_csv(tables_dir / "table_cbioportal_endpoint_crosscheck.csv")
    cbio_meta = cbio_tbl.iloc[0].to_dict()

    metrics_export = {
        "study_title": "TCGA-PRAD Endpoint Harmonization & Molecular Prediction Benchmark",
        "primary_definitive_cohort": {
            "total_evaluable_patients": len(data_definitive),
            "events_positive": int(y_def.sum()),
            "non_events_double_negative": int((y_def == 0).sum()),
            "asymmetric_negative_excluded": 46,
            "indeterminate_excluded": 5,
            "paired_molecular_n": len(y_def_p),
            "paired_molecular_events": int(y_def_p.sum()),
            "paired_molecular_non_events": int((y_def_p == 0).sum()),
            "clinical_only_auroc": res_def_paired_clin["auroc"],
            "clinical_plus_rna_auroc": res_def_paired_mol["auroc"],
            "paired_delta_auroc": paired_diff_pure["delta_auroc"],
            "paired_delta_auroc_95ci": [paired_diff_pure["delta_auroc_ci_low"], paired_diff_pure["delta_auroc_ci_high"]],
            "paired_delta_auroc_p_value": paired_diff_pure["delta_auroc_p_value"],
        },
        "exploratory_composite_cohort": {
            "total_evaluable_patients": len(data_exploratory),
            "events_positive": int(y_exp.sum()),
            "non_events": int((y_exp == 0).sum()),
            "paired_molecular_n": len(y_exp_p),
            "paired_molecular_events": int(y_exp_p.sum()),
            "clinical_only_auroc": res_exp_paired_clin["auroc"],
            "clinical_plus_rna_auroc": res_exp_paired_mol["auroc"],
            "paired_delta_auroc": paired_diff_exp["delta_auroc"],
            "paired_delta_auroc_95ci": [paired_diff_exp["delta_auroc_ci_low"], paired_diff_exp["delta_auroc_ci_high"]],
            "paired_delta_auroc_p_value": paired_diff_exp["delta_auroc_p_value"],
        },
        "cohort": {
            "total_patients": len(data),
            "composite_recurrence_events": int(data["outcome_recurrence_or_progression"].sum()),
            "dated_recurrence_events": int(recurrence["survival_event_indicator"].eq(1.0).sum()),
            "undated_recurrence_events": int(recurrence["undated_recurrence_flag"].sum()),
            "bcr_events": int(recurrence["flag_bcr_positive"].sum()),
            "nte_events": int(recurrence["flag_nte_positive"].sum()),
            "rna_seq_profiles_available": int(complete_mask_all.sum()),
            "rna_seq_missing_count": int((~complete_mask_all).sum()),
            "evaluable_survival_cases": len(exp_surv_df),
        },
        "survival": {
            "evaluable_cases": len(exp_surv_df),
            "median_rfs_days": exp_median,
            "rfs_1yr": exp_r1,
            "rfs_3yr": exp_r3,
            "rfs_5yr": exp_r5,
            "cox_apparent_concordance_index": exp_app_c,
            "cox_cv_concordance_index": exp_cv_c,
            "strict_evaluable_cases": len(strict_surv_df),
            "strict_cox_apparent_concordance_index": strict_app_c,
            "strict_cox_cv_concordance_index": strict_cv_c,
        },
        "cbioportal_concordance": {
            "bcr_to_bcr_agreement": 1.0,
            "bcr_to_bcr_evaluable_n": 269,
            "definitive_dfs_agreement": 0.923372,
            "definitive_dfs_evaluable_n": 261,
            "definitive_dfs_discordant_cases": 20,
            "exploratory_dfs_agreement": 0.934641,
            "exploratory_dfs_evaluable_n": 306,
            "exploratory_dfs_discordant_cases": 20,
            "historical_benchmark_discordant_cases": int(cbio_meta["discordant_cases"]),
        },
        "cohort_summary": {
            "total_patients": len(data),
            "composite_recurrence_events": int(data["outcome_recurrence_or_progression"].sum()),
            "dated_recurrence_events": int(recurrence["survival_event_indicator"].eq(1.0).sum()),
            "undated_recurrence_events": int(recurrence["undated_recurrence_flag"].sum()),
            "bcr_events": int(recurrence["flag_bcr_positive"].sum()),
            "nte_events": int(recurrence["flag_nte_positive"].sum()),
            "rna_seq_profiles_available": int(complete_mask_all.sum()),
            "rna_seq_missing_count": int((~complete_mask_all).sum()),
        },
        "validation_models": full_table.to_dict(orient="records"),
        "paired_molecular_comparisons": paired_molecular_table.to_dict(orient="records"),
    }

    with open(artifacts_dir / "tcga_prad_metrics.json", "w") as f:
        json.dump(metrics_export, f, indent=2)

    # Mirror to web directories
    for dest_dir in [PROJECT_ROOT / "web" / "static" / "research", PROJECT_ROOT / "web" / "research" / "artifacts"]:
        dest_dir.mkdir(parents=True, exist_ok=True)
        with open(dest_dir / "tcga_prad_metrics.json", "w") as f:
            json.dump(metrics_export, f, indent=2)

    for dest_fig in [PROJECT_ROOT / "web" / "static" / "figures", PROJECT_ROOT / "web" / "research" / "figures"]:
        dest_fig.mkdir(parents=True, exist_ok=True)
        shutil.copy(figures_dir / "prad_calibration.png", dest_fig / "prad_calibration.png")
        shutil.copy(figures_dir / "prad_km_survival.png", dest_fig / "prad_km_survival.png")
        shutil.copy(figures_dir / "prad_paired_roc.png", dest_fig / "prad_paired_roc.png")

    print("\n=== TCGA-PRAD SCIENTIFIC QUALITY GATE VALIDATION COMPLETE ===")
    print(f"Primary Definitive Cohort: N={len(data_definitive)} (72 Events, 191 Double-Negative Non-Events)")
    print(f"Paired Definitive Molecular Cohort: N={len(y_def_p)} (45 Events, 143 Non-Events)")
    print(f"  Pure Clinical Baseline AUROC: {res_def_paired_clin['auroc']:.4f}")
    print(f"  Clinical + 13-Gene RNA AUROC: {res_def_paired_mol['auroc']:.4f}")
    print(f"  Paired ΔAUROC: {paired_diff_pure['delta_auroc']:+.4f} (95% CI: [{paired_diff_pure['delta_auroc_ci_low']:+.4f}, {paired_diff_pure['delta_auroc_ci_high']:+.4f}], p={paired_diff_pure['delta_auroc_p_value']:.4f})")
    print(f"Strict Survival Analysis (N={len(strict_surv_df)}): 1-yr RFS={strict_r1:.3f}, 3-yr RFS={strict_r3:.3f}, 5-yr RFS={strict_r5:.3f}")
    print(f"  Strict Cox PH C-index: Apparent={strict_app_c:.3f}, 5-Fold CV={strict_cv_c:.3f}")
    print(f"Exploratory Survival Analysis (N={len(exp_surv_df)}): 1-yr RFS={exp_r1:.3f}, 3-yr RFS={exp_r3:.3f}, 5-yr RFS={exp_r5:.3f}")
    print(f"  Exploratory Cox PH C-index: Apparent={exp_app_c:.3f}, 5-Fold CV={exp_cv_c:.3f}")


if __name__ == "__main__":
    main()
