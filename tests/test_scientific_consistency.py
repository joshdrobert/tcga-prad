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
    flow = pd.read_csv(REPO_ROOT / "research" / "v2" / "COHORT_FLOW.csv")
    metrics_path = REPO_ROOT / "research" / "artifacts" / "tcga_prad_metrics.json"
    assert metrics_path.exists()
    with open(metrics_path) as f:
        metrics = json.load(f)

    # 1. Total Radical Prostatectomy Patients
    total_flow = int(flow.loc[flow["step"] == "Total TCGA-PRAD Radical Prostatectomy Patients", "count"].iloc[0])
    assert total_flow == 314
    assert metrics["cohort"]["total_patients"] == 314

    # 2. Primary Definitive Cohort
    def_flow = int(flow.loc[flow["step"] == "Primary Definitive Binary Cohort (Positives + Double Negatives)", "count"].iloc[0])
    assert def_flow == 263
    assert metrics["primary_definitive_cohort"]["total_evaluable_patients"] == 263
    assert metrics["primary_definitive_cohort"]["events_positive"] == 72
    assert metrics["primary_definitive_cohort"]["non_events_double_negative"] == 191

    # 3. Paired Molecular Cohort
    paired_tbl = pd.read_csv(REPO_ROOT / "tables" / "table_prad_paired_complete_case_validation.csv")
    def_paired_row = paired_tbl[paired_tbl["analysis_tier"] == "PRIMARY_DEFINITIVE"].iloc[0]
    assert def_paired_row["sample_size_n"] == 188
    assert def_paired_row["events_n"] == 45
    assert def_paired_row["non_events_n"] == 143
    assert metrics["primary_definitive_cohort"]["paired_molecular_n"] == 188

    # 4. Survival Cohorts
    surv_tbl = pd.read_csv(REPO_ROOT / "tables" / "table_prad_survival_analysis.csv")
    strict_row = surv_tbl[surv_tbl["analysis_definition"].str.contains("Strict")].iloc[0]
    exp_row = surv_tbl[surv_tbl["analysis_definition"].str.contains("Exploratory")].iloc[0]
    assert strict_row["evaluable_n"] == 257
    assert strict_row["dated_events"] == 66
    assert strict_row["censored_non_events"] == 191
    assert exp_row["evaluable_n"] == 303
    assert exp_row["dated_events"] == 66
    assert exp_row["censored_non_events"] == 237

    # 5. Excluded Cases
    bias_tbl = pd.read_csv(REPO_ROOT / "research" / "v2" / "SELECTION_BIAS_AUDIT.csv")
    size_row = bias_tbl[bias_tbl["characteristic"] == "Cohort Size (N)"].iloc[0]
    assert size_row["excluded_cohort_n11"] == 11


def test_auroc_and_discrimination_consistency():
    """Verify that AUROCs and paired differences match across tables, JSON artifacts, and abstracts."""
    paired_tbl = pd.read_csv(REPO_ROOT / "tables" / "table_prad_paired_complete_case_validation.csv")
    def_row = paired_tbl[paired_tbl["analysis_tier"] == "PRIMARY_DEFINITIVE"].iloc[0]

    with open(REPO_ROOT / "research" / "artifacts" / "tcga_prad_metrics.json") as f:
        metrics = json.load(f)

    # Pure Clinical AUROC: 0.6931
    assert abs(def_row["baseline_auroc"] - 0.693084) < 1e-4
    assert abs(metrics["primary_definitive_cohort"]["clinical_only_auroc"] - 0.693084) < 1e-4

    # Clin+RNA AUROC: 0.6476
    assert abs(def_row["augmented_auroc"] - 0.647552) < 1e-4
    assert abs(metrics["primary_definitive_cohort"]["clinical_plus_rna_auroc"] - 0.647552) < 1e-4

    # Paired Delta AUROC: -0.0455
    assert abs(def_row["delta_auroc"] - (-0.045532)) < 1e-4
    assert abs(metrics["primary_definitive_cohort"]["paired_delta_auroc"] - (-0.045532)) < 1e-4

    # Delta p-value: not statistically significant (p > 0.05)
    assert def_row["delta_auroc_p_value"] > 0.05
    assert def_row["delta_auroc_95ci_low"] < 0.0 < def_row["delta_auroc_95ci_high"]

    # Verify abstract text matches within 0.001
    abstract_text = (REPO_ROOT / "conference" / "oncology" / "V2_ABSTRACT.md").read_text()
    assert "0.693" in abstract_text
    assert "0.648" in abstract_text
    assert "-0.046" in abstract_text


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

    assert abs(dfs_row["agreement_fraction"] - 0.923372) < 1e-4
    assert dfs_row["discordant_cases"] == 20
    assert dfs_row["local_only_events"] == 20
    assert dfs_row["cbio_only_events"] == 0


def test_predictor_timing_integrity():
    """Verify that XML psa_value is classified as post-operative surveillance rather than baseline."""
    timing_tbl = pd.read_csv(REPO_ROOT / "research" / "v2" / "PREDICTOR_TIMING_AUDIT.csv")
    psa_row = timing_tbl[timing_tbl["predictor"] == "psa_value_numeric"].iloc[0]
    assert psa_row["clinical_phase"] == "Post-operative Surveillance"
    assert psa_row["timing_relative_to_surgery_days"] > 100.0  # Median 511 days
