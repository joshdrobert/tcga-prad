from __future__ import annotations

import json
from pathlib import Path
import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]


def test_no_postoperative_psa_in_baseline_models():
    """Ensure postoperative surveillance PSA is strictly excluded from baseline Cox survival models."""
    cox_strict_file = REPO_ROOT / "tables" / "table_prad_cox_summary_strict.csv"
    assert cox_strict_file.exists(), "Strict Cox summary table must exist"
    cox_summary = pd.read_csv(cox_strict_file)
    covariates = cox_summary["covariate"].tolist()

    assert "psa_value_numeric" not in covariates, "Postoperative PSA must NOT be included in baseline Cox model"
    assert "psa_value" not in covariates, "Postoperative PSA must NOT be included in baseline Cox model"
    assert "age_at_diagnosis_years" in covariates, "Age must be in baseline Cox model"
    assert "gleason_score_numeric" in covariates, "Gleason score must be in baseline Cox model"

    # Verify timing audit classification
    timing_file = REPO_ROOT / "research" / "v3" / "PREDICTOR_TIMING_AUDIT.csv"
    assert timing_file.exists()
    timing_df = pd.read_csv(timing_file)
    psa_row = timing_df[timing_df["predictor"] == "psa_value_numeric"].iloc[0]
    assert psa_row["clinical_phase"] == "Post-operative Surveillance"
    assert "strictly excluded from baseline cox models" in psa_row["timing_interpretation"].lower()


def test_time_origin_indexed_to_initial_pathologic_diagnosis():
    """Verify that survival analysis explicitly indexes time to initial pathologic diagnosis (biopsy Day 0)."""
    rec_file = REPO_ROOT / "data_processed" / "tcga_prad_recurrence_cohort.csv"
    assert rec_file.exists()
    rec_df = pd.read_csv(rec_file)

    assert "duration_days_from_diagnosis" in rec_df.columns, "Explicit duration_days_from_diagnosis column must exist"
    assert "survival_duration_days" in rec_df.columns

    surv_file = REPO_ROOT / "tables" / "table_prad_survival_analysis.csv"
    assert surv_file.exists()
    surv_df = pd.read_csv(surv_file)
    for _, row in surv_df.iterrows():
        assert row["time_origin"] == "Days from Initial Pathologic Diagnosis"


def test_second_primary_cancers_not_treated_as_prostate_recurrence():
    """Verify that non-prostatic second primary tumors (TCGA-CH-5791, TCGA-HC-7080) are separated from prostate recurrence."""
    rec_file = REPO_ROOT / "data_processed" / "tcga_prad_recurrence_cohort.csv"
    rec_df = pd.read_csv(rec_file)

    second_primaries = ["TCGA-CH-5791", "TCGA-HC-7080"]
    sub_df = rec_df[rec_df["submitter_id"].isin(second_primaries)]

    assert len(sub_df) == 2, "Both second primary cases must be present"
    for _, row in sub_df.iterrows():
        assert row["composite_status_detail"] == "SECOND_PRIMARY_MALIGNANCY", (
            f"Patient {row['submitter_id']} must be adjudicated as SECOND_PRIMARY_MALIGNANCY"
        )
        assert pd.isna(row["outcome_prostate_specific_recurrence"]) or row["outcome_prostate_specific_recurrence"] == 0.0, (
            f"Patient {row['submitter_id']} must NOT have outcome_prostate_specific_recurrence == 1.0"
        )
        assert pd.isna(row["outcome_definitive_recurrence"]) or row["outcome_definitive_recurrence"] == 0.0, (
            f"Patient {row['submitter_id']} must NOT have outcome_definitive_recurrence == 1.0"
        )
        assert pd.isna(row["survival_event_indicator"]) or row["survival_event_indicator"] == 0.0, (
            f"Patient {row['submitter_id']} must NOT have survival_event_indicator == 1.0"
        )

    # Verify dedicated artifact exists
    sec_file = REPO_ROOT / "research" / "v3" / "SECOND_PRIMARY_CASES.csv"
    assert sec_file.exists()
    sec_df = pd.read_csv(sec_file)
    assert len(sec_df) == 2
    assert set(sec_df["submitter_id"]) == set(second_primaries)


def test_no_overclaiming_proportional_hazards_or_molecular_significance():
    """Verify that Schoenfeld test is reported accurately and no unsupported molecular claims exist."""
    files_to_check = [
        REPO_ROOT / "CLAIM_LEDGER.csv",
        REPO_ROOT / "tables" / "table_prad_survival_analysis.csv",
        REPO_ROOT / "research" / "artifacts" / "tcga_prad_metrics.json",
        REPO_ROOT / "conference" / "oncology" / "V2_ABSTRACT.md",
        REPO_ROOT / "web" / "index.html",
    ]

    for filepath in files_to_check:
        if not filepath.exists():
            continue
        content = filepath.read_text().lower()

        # Disallowed over-claiming language for Schoenfeld test:
        assert "confirmed proportional hazards" not in content, f"{filepath} must not over-claim confirmed proportional hazards"
        assert "confirms proportional hazards" not in content, f"{filepath} must not over-claim confirms proportional hazards"
        assert "proven proportional hazards" not in content, f"{filepath} must not over-claim proven proportional hazards"

        # Disallowed molecular superiority claims:
        assert "statistically significant improvement from rna" not in content
        assert "expression panel significantly improved" not in content
        assert "superior molecular discrimination" not in content
