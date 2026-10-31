from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.utils import concordance_index
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
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


def compute_calibration_slope(y: np.ndarray, p: np.ndarray) -> tuple[float, float]:
    clipped = np.clip(p, 1e-6, 1.0 - 1e-6)
    logit = np.log(clipped / (1.0 - clipped)).reshape(-1, 1)
    mod = LogisticRegression(C=1e6, solver="lbfgs").fit(logit, y)
    return float(mod.intercept_[0]), float(mod.coef_[0, 0])


def evaluate_feature_set(df: pd.DataFrame, features: list[str], outcome_col: str = "outcome_recurrence_or_progression") -> tuple[dict, np.ndarray]:
    X = df[features].copy()
    y = df[outcome_col].astype(int).values
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
    metrics["predictor_count"] = len(features)
    return metrics, probs


def evaluate_cv_cox_model(df: pd.DataFrame, duration_col: str, event_col: str, covariates: list[str]) -> tuple[float, float, pd.DataFrame]:
    """Fit Cox PH model and compute both apparent and 5-fold cross-validated C-index."""
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
    return apparent_c, cv_c, summary_df


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

    clinical_features = ["age_at_diagnosis_years", "pathologic_t", "pathologic_n", "gleason_score_numeric", "psa_value_numeric"]
    all_features = clinical_features + expr_cols

    # Subset paired complete cases (patients with non-null RNA data)
    complete_mask = data[expr_cols].notna().all(axis=1)
    data_complete = data[complete_mask].copy()

    # Subset high confidence cases (excluding 5 indeterminate missing cases and 46 asymmetric negatives)
    high_conf_mask = data["high_confidence_outcome_flag"] == 1
    data_high_conf = data[high_conf_mask].copy()

    # 1. Full cohort evaluation (n=314)
    res_clin_full, prob_clin_full = evaluate_feature_set(data, clinical_features)
    res_mol_full, prob_mol_full = evaluate_feature_set(data, all_features)

    # 2. Paired complete-case evaluation (n=228)
    res_clin_comp, prob_clin_comp = evaluate_feature_set(data_complete, clinical_features)
    res_mol_comp, prob_mol_comp = evaluate_feature_set(data_complete, all_features)

    # 3. High-confidence outcome sensitivity cohort (n=263)
    res_clin_high, prob_clin_high = evaluate_feature_set(data_high_conf, clinical_features)
    res_mol_high, prob_mol_high = evaluate_feature_set(data_high_conf, all_features)

    full_table = pd.DataFrame(
        [
            {"cohort": "Full (n=314)", "model": "clinical_only", **res_clin_full},
            {"cohort": "Full (n=314)", "model": "clinical_expression_imputed", **res_mol_full},
            {"cohort": "Paired Complete (n=228)", "model": "clinical_only_complete_cases", **res_clin_comp},
            {"cohort": "Paired Complete (n=228)", "model": "clinical_expression_complete_cases", **res_mol_comp},
            {"cohort": "High Confidence (n=263)", "model": "clinical_only_high_conf", **res_clin_high},
            {"cohort": "High Confidence (n=263)", "model": "clinical_expression_high_conf", **res_mol_high},
        ]
    )
    full_table.to_csv(tables_dir / "table_prad_rigorous_validation.csv", index=False)
    full_table.to_csv(RESEARCH_V2_DIR / "MODEL_COMPARISONS.csv", index=False)

    paired_table = pd.DataFrame(
        [
            {"cohort": "Paired Complete (n=228)", "model": "clinical_only_complete_cases", **res_clin_comp},
            {"cohort": "Paired Complete (n=228)", "model": "clinical_expression_complete_cases", **res_mol_comp},
        ]
    )
    paired_table.to_csv(tables_dir / "table_prad_paired_complete_case_validation.csv", index=False)

    # 4. Longitudinal Survival Analysis (Primary Defensible Population)
    # Valid continuous-time cases: 66 dated recurrence events + 242 censored non-events (total n=308)
    surv_df = recurrence[recurrence["survival_duration_days"].notna() & (recurrence["survival_duration_days"] > 0)].copy()
    surv_duration = surv_df["survival_duration_days"].values
    surv_event = surv_df["survival_event_indicator"].astype(int).values

    kmf = KaplanMeierFitter()
    kmf.fit(surv_duration, surv_event, label="Recurrence-Free Survival (Dated Events, n=308)")

    rfs_1yr = float(kmf.predict(365.25))
    rfs_3yr = float(kmf.predict(365.25 * 3))
    rfs_5yr = float(kmf.predict(365.25 * 5))
    median_rfs = float(kmf.median_survival_time_) if pd.notna(kmf.median_survival_time_) else np.nan

    # Cox Proportional Hazards Model
    cox_covs = ["age_at_diagnosis_years", "gleason_score_numeric", "psa_value_numeric"]
    apparent_c, cv_c, cox_summary = evaluate_cv_cox_model(surv_df, "survival_duration_days", "survival_event_indicator", cox_covs)

    cox_summary.to_csv(tables_dir / "table_prad_survival_analysis.csv", index=False)
    
    surv_results_tbl = pd.DataFrame(
        [
            {
                "analysis": "Primary Defensible Survival (Dated Events Only)",
                "evaluable_n": len(surv_df),
                "dated_events": int(surv_event.sum()),
                "undated_events_excluded": int(recurrence["undated_recurrence_flag"].sum()),
                "censored_cases": int((surv_event == 0).sum()),
                "rfs_1yr": rfs_1yr,
                "rfs_3yr": rfs_3yr,
                "rfs_5yr": rfs_5yr,
                "median_rfs_days": median_rfs,
                "cox_apparent_c_index": apparent_c,
                "cox_cv_c_index": cv_c,
            }
        ]
    )
    surv_results_tbl.to_csv(RESEARCH_V2_DIR / "SURVIVAL_ANALYSIS.csv", index=False)

    # 5. Generate Publication Figures
    # Calibration Plot
    fig, ax = plt.subplots(figsize=(6, 5), dpi=300)
    for name, p, y in [
        ("Clinical Only (n=314)", prob_clin_full, data["outcome_recurrence_or_progression"].values),
        ("Clinical+Expr Imputed (n=314)", prob_mol_full, data["outcome_recurrence_or_progression"].values),
        ("Clinical Only Paired (n=228)", prob_clin_comp, data_complete["outcome_recurrence_or_progression"].values),
        ("Clinical+Expr Paired (n=228)", prob_mol_comp, data_complete["outcome_recurrence_or_progression"].values),
    ]:
        df_cal = pd.DataFrame({"p": p, "y": y})
        df_cal["bin"] = pd.qcut(df_cal["p"], 8, duplicates="drop")
        binned = df_cal.groupby("bin", observed=True).agg(mean_p=("p", "mean"), mean_y=("y", "mean"))
        ax.plot(binned["mean_p"], binned["mean_y"], marker="o", label=name)
    ax.plot([0, 1], [0, 1], "k--", label="Ideal Calibration")
    ax.set_xlabel("Mean Predicted Probability")
    ax.set_ylabel("Observed Event Rate")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend(fontsize=7, loc="upper left")
    ax.set_title("5-Fold Cross-Validation Calibration Curves", fontsize=10)
    fig.tight_layout()
    fig.savefig(figures_dir / "prad_calibration.png")
    plt.close(fig)

    # Kaplan-Meier Curve
    fig, ax = plt.subplots(figsize=(6, 4.5), dpi=300)
    kmf.plot_survival_function(ax=ax, color="#1e3a8a", ci_show=True)
    ax.set_xlabel("Time post-radical prostatectomy (Days)")
    ax.set_ylabel("Recurrence-Free Survival Probability")
    ax.set_title("TCGA-PRAD Kaplan-Meier Recurrence-Free Survival (V2 Audited)", fontsize=10)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(figures_dir / "prad_km_survival.png")
    plt.close(fig)

    # 6. Load cBioPortal Crosscheck Dynamically (NO HARDCODING)
    cbio_tbl = pd.read_csv(tables_dir / "table_cbioportal_endpoint_crosscheck.csv")
    cbio_meta = cbio_tbl.iloc[0].to_dict()

    metrics_export = {
        "study_title": "TCGA-PRAD Endpoint Harmonization & Molecular Prediction Benchmark",
        "cohort": {
            "total_patients": len(data),
            "composite_recurrence_events": int(data["outcome_recurrence_or_progression"].sum()),
            "dated_recurrence_events": int(surv_event.sum()),
            "undated_recurrence_events": int(recurrence["undated_recurrence_flag"].sum()),
            "bcr_events": int(recurrence["flag_bcr_positive"].sum()),
            "nte_events": int(recurrence["flag_nte_positive"].sum()),
            "rna_seq_profiles_available": int(complete_mask.sum()),
            "rna_seq_missing_count": int((~complete_mask).sum()),
            "rna_seq_missing_fraction": float((~complete_mask).mean()),
            "evaluable_survival_cases": len(surv_df),
        },
        "survival": {
            "median_rfs_days": median_rfs,
            "rfs_1yr": rfs_1yr,
            "rfs_3yr": rfs_3yr,
            "rfs_5yr": rfs_5yr,
            "cox_apparent_concordance_index": apparent_c,
            "cox_cv_concordance_index": cv_c,
        },
        "validation_models": full_table.to_dict(orient="records"),
        "cbioportal_crosscheck": {
            "evaluable_matched": int(cbio_meta["evaluable_matched"]),
            "agreement_fraction": float(cbio_meta["agreement_fraction"]),
            "concordant_non_events": int(cbio_meta["concordant_non_events"]),
            "concordant_events": int(cbio_meta["concordant_events"]),
            "discordant_cases": int(cbio_meta["discordant_cases"]),
            "local_only_events": int(cbio_meta["local_only_events"]),
            "cbio_only_events": int(cbio_meta["cbio_only_events"]),
            "discordance_reason": "Local composite captured documented clinical New Tumor Events / Distant Progression without recorded biochemical PSA elevation",
        },
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
        import shutil
        shutil.copy(figures_dir / "prad_calibration.png", dest_fig / "prad_calibration.png")
        shutil.copy(figures_dir / "prad_km_survival.png", dest_fig / "prad_km_survival.png")

    print("\n=== TCGA-PRAD V2 RIGOROUS VALIDATION COMPLETE ===")
    print(full_table[["cohort", "model", "n", "events", "auroc", "auprc", "calibration_slope"]].to_string(index=False))
    print(f"\nSurvival Analysis (n={len(surv_df)}): 1-yr RFS={rfs_1yr:.3f}, 3-yr RFS={rfs_3yr:.3f}, 5-yr RFS={rfs_5yr:.3f}")
    print(f"Cox PH Concordance: Apparent={apparent_c:.3f}, 5-Fold CV={cv_c:.3f}")
    print(f"Artifact written to {artifacts_dir / 'tcga_prad_metrics.json'}")


if __name__ == "__main__":
    main()
