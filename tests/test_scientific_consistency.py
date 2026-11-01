from __future__ import annotations

import json
import re
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]


def test_cohort_denominator_consistency():
    """Verify that cohort denominators are consistent across all artifacts, tables, and web docs."""
    flow = pd.read_csv(REPO_ROOT / "research" / "v3" / "COHORT_FLOW.csv")
    metrics_path = REPO_ROOT / "research" / "artifacts" / "tcga_prad_metrics.json"
    assert metrics_path.exists()
    with open(metrics_path) as f:
        metrics = json.load(f)

    # 1. Total Radical Prostatectomy Patients
    total_flow = int(flow.loc[flow["step"] == "Total TCGA-PRAD Radical Prostatectomy Patients", "count"].iloc[0])
    assert total_flow == 314
    assert metrics["cohort"]["total_patients"] == 314

    # 2. Primary Definitive Cohort
    def_flow = int(flow.loc[flow["step"].str.contains("Primary Definitive Prostate-Specific Cohort"), "count"].iloc[0])
    assert def_flow == 261
    assert metrics["primary_definitive_cohort"]["total_evaluable_patients"] == 261
    assert metrics["primary_definitive_cohort"]["events_positive"] == 70
    assert metrics["primary_definitive_cohort"]["non_events_double_negative"] == 191

    # 3. Paired Molecular Cohort
    paired_tbl = pd.read_csv(REPO_ROOT / "tables" / "table_prad_paired_complete_case_validation.csv")
    def_paired_row = paired_tbl[paired_tbl["analysis_tier"] == "PRIMARY_DEFINITIVE"].iloc[0]
    assert def_paired_row["sample_size_n"] == 187
    assert def_paired_row["events_n"] == 44
    assert def_paired_row["non_events_n"] == 143
    assert metrics["primary_definitive_cohort"]["paired_molecular_n"] == 187

    # 4. Survival Cohorts (Strict Primary + Permissive Sensitivity)
    surv_tbl = pd.read_csv(REPO_ROOT / "tables" / "table_prad_survival_analysis.csv")
    strict_row = surv_tbl[surv_tbl["analysis_definition"].str.contains("Strict")].iloc[0]
    exp_row = surv_tbl[surv_tbl["analysis_definition"].str.contains("Exploratory")].iloc[0]
    assert strict_row["evaluable_n"] == 248
    assert strict_row["dated_events"] == 57
    assert strict_row["censored_non_events"] == 191
    assert exp_row["evaluable_n"] == 294
    assert exp_row["dated_events"] == 57
    assert exp_row["censored_non_events"] == 237

    # Permissive sensitivity check
    perm_strict = surv_tbl[surv_tbl["survival_tier"] == "PERMISSIVE_SENSITIVITY_DOUBLE_NEG"].iloc[0]
    assert perm_strict["evaluable_n"] == 255
    assert perm_strict["dated_events"] == 64

    # 5. Excluded Cases (20 strict total: 13 undated prostate recurrence + 2 second primaries + 5 indeterminate)
    bias_tbl = pd.read_csv(REPO_ROOT / "research" / "v3.2" / "SELECTION_BIAS_AUDIT.csv")
    size_row = bias_tbl[bias_tbl["characteristic"] == "Cohort Size (N)"].iloc[0]
    assert size_row["excluded_cohort_n20"] == 20


def test_auroc_and_discrimination_consistency():
    """Verify that AUROCs and paired differences match across tables, JSON artifacts, and abstracts."""
    paired_tbl = pd.read_csv(REPO_ROOT / "tables" / "table_prad_paired_complete_case_validation.csv")
    def_row = paired_tbl[paired_tbl["analysis_tier"] == "PRIMARY_DEFINITIVE"].iloc[0]

    with open(REPO_ROOT / "research" / "artifacts" / "tcga_prad_metrics.json") as f:
        metrics = json.load(f)

    # Pure Clinical AUROC: 0.7104
    assert abs(def_row["baseline_auroc"] - 0.710426) < 1e-4
    assert abs(metrics["primary_definitive_cohort"]["clinical_only_auroc"] - 0.710426) < 1e-4

    # Clin+RNA AUROC: 0.7028
    assert abs(def_row["augmented_auroc"] - 0.702797) < 1e-4
    assert abs(metrics["primary_definitive_cohort"]["clinical_plus_rna_auroc"] - 0.702797) < 1e-4

    # Paired Delta AUROC: -0.0076
    assert abs(def_row["delta_auroc"] - (-0.007629)) < 1e-4
    assert abs(metrics["primary_definitive_cohort"]["paired_delta_auroc"] - (-0.007629)) < 1e-4

    # Paired Delta Brier: -0.0001
    assert abs(def_row["delta_brier"] - (-0.000125)) < 1e-4
    assert def_row["delta_brier_95ci_low"] < 0.0 < def_row["delta_brier_95ci_high"]

    # Delta bootstrap CI spans zero (no statistically significant difference)
    assert def_row["delta_auroc_95ci_low"] < 0.0 < def_row["delta_auroc_95ci_high"]
    if "delta_auroc_p_value" in def_row and pd.notna(def_row["delta_auroc_p_value"]):
        assert def_row["delta_auroc_p_value"] > 0.05

    # Cross-verify with rigorous validation table
    rigorous_tbl = pd.read_csv(REPO_ROOT / "tables" / "table_prad_rigorous_validation.csv")
    rig_clin = rigorous_tbl[(rigorous_tbl["cohort"] == "Paired Definitive Molecular (n=187)") & (rigorous_tbl["model"].str.startswith("Pure Baseline Clinical ("))].iloc[0]
    rig_mol = rigorous_tbl[(rigorous_tbl["cohort"] == "Paired Definitive Molecular (n=187)") & (rigorous_tbl["model"].str.contains("13-Gene Expression Panel"))].iloc[0]

    assert abs(rig_clin["auroc"] - def_row["baseline_auroc"]) < 1e-6
    assert abs(rig_mol["auroc"] - def_row["augmented_auroc"]) < 1e-6
    assert abs(rig_clin["brier"] - def_row["baseline_brier"]) < 1e-6
    assert abs(rig_mol["brier"] - def_row["augmented_brier"]) < 1e-6

    # Verify abstract text matches within 0.001
    abstract_text = (REPO_ROOT / "conference" / "oncology" / "V2_ABSTRACT.md").read_text()
    assert f"{def_row['baseline_auroc']:.3f}" in abstract_text
    assert f"{def_row['augmented_auroc']:.3f}" in abstract_text
    assert f"{def_row['delta_auroc']:+.3f}" in abstract_text or f"{def_row['delta_auroc']:.3f}" in abstract_text
    assert "0.1673" in abstract_text
    assert "0.1672" in abstract_text


def test_unsupported_significance_wording():
    """Verify that no file claims statistically significant molecular improvement in TCGA-PRAD."""
    files_to_check = [
        REPO_ROOT / "CLAIM_LEDGER.csv",
        REPO_ROOT / "conference" / "oncology" / "V2_ABSTRACT.md",
        REPO_ROOT / "web" / "index.html",
    ]
    for filepath in files_to_check:
        if not filepath.exists():
            continue
        content = filepath.read_text().lower()
        # Phrases indicating an unsupported claim of molecular superiority
        assert "statistically significant improvement from rna" not in content
        assert "expression panel significantly improved" not in content
        assert "superior molecular discrimination" not in content


def test_invalidated_claims_status():
    """Verify that earlier or unadjudicated claims are properly labeled as INVALIDATED or BENCHMARK."""
    claim_tbl = pd.read_csv(REPO_ROOT / "CLAIM_LEDGER.csv")
    hist_claims = claim_tbl[claim_tbl["claim_id"] == "CLM-PRAD-13"]
    assert len(hist_claims) == 1
    text_combined = (str(hist_claims.iloc[0]["comparator_metric"]) + " " + str(hist_claims.iloc[0]["notes"])).lower()
    assert "benchmark" in text_combined


def test_cbioportal_crosscheck_consistency():
    """Verify explicit endpoint concordance definitions match across tables and JSON artifacts."""
    multi_tbl = pd.read_csv(REPO_ROOT / "tables" / "table_cbioportal_endpoint_crosscheck_multidef.csv")
    bcr_row = multi_tbl[multi_tbl["comparison_definition"].str.contains("BCR-to-BCR")].iloc[0]
    dfs_row = multi_tbl[multi_tbl["comparison_definition"].str.contains("Definitive Recurrence")].iloc[0]

    assert bcr_row["agreement_fraction"] == 1.0
    assert bcr_row["discordant_cases"] == 0

    assert abs(dfs_row["agreement_fraction"] - 0.922780) < 1e-4
    assert dfs_row["discordant_cases"] == 20
    assert dfs_row["local_only_events"] == 20
    assert dfs_row["cbio_only_events"] == 0


def test_predictor_timing_integrity():
    """Verify that XML psa_value is classified as post-operative surveillance rather than baseline."""
    timing_tbl = pd.read_csv(REPO_ROOT / "research" / "v3" / "PREDICTOR_TIMING_AUDIT.csv")
    psa_row = timing_tbl[timing_tbl["predictor"] == "psa_value_numeric"].iloc[0]
    assert "post-operative surveillance" in psa_row["clinical_phase"].lower()
    assert "511 days" in psa_row["reference_origin_days"].lower()
