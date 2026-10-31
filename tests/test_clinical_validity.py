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
    assert "post-operative surveillance" in psa_row["clinical_phase"].lower()
    assert "strictly excluded from baseline cox" in psa_row["immortal_time_and_landmark_risk"].lower()


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


def test_strict_json_rfc8259_compliance():
    """Verify metrics JSON strictly adheres to RFC 8259 (no Infinity or NaN)."""
    metrics_file = REPO_ROOT / "research" / "artifacts" / "tcga_prad_metrics.json"
    assert metrics_file.exists(), "Metrics JSON artifact must exist"
    raw_text = metrics_file.read_text()

    assert "Infinity" not in raw_text, "Strict RFC 8259 JSON cannot contain 'Infinity'"
    assert "NaN" not in raw_text, "Strict RFC 8259 JSON cannot contain 'NaN'"

    def fail_on_nonstandard(x):
        raise ValueError(f"Nonstandard JSON token detected: {x}")

    metrics = json.loads(raw_text, parse_constant=fail_on_nonstandard)
    assert metrics["survival"]["median_rfs_days"] is None, "Median RFS days must be null (None)"
    assert metrics["survival"]["median_rfs_status"] == "not_reached", "Median RFS status must be 'not_reached'"
    assert metrics["survival"]["strict_median_rfs_days"] is None
    assert metrics["survival"]["strict_median_rfs_status"] == "not_reached"


def test_event_level_provenance_pairing():
    """Verify EVENT_LEVEL_PROVENANCE.csv unrolls and pairs event types, dates, and sources without decoupling."""
    prov_file = REPO_ROOT / "research" / "v3" / "EVENT_LEVEL_PROVENANCE.csv"
    assert prov_file.exists(), "EVENT_LEVEL_PROVENANCE.csv must exist in research/v3"
    prov_df = pd.read_csv(prov_file)

    required_cols = [
        "submitter_id",
        "source_file",
        "event_category",
        "event_type",
        "event_days_from_diagnosis",
        "prostate_specific",
        "adjudication_class",
        "adjudication_rationale",
    ]
    for col in required_cols:
        assert col in prov_df.columns, f"Column '{col}' must be in provenance table"

    # Verify second primaries are non-prostate specific
    sec_primaries = prov_df[prov_df["submitter_id"].isin(["TCGA-CH-5791", "TCGA-HC-7080"])]
    assert len(sec_primaries) >= 2
    for _, row in sec_primaries.iterrows():
        assert str(row["prostate_specific"]).lower() == "false"

    # Verify cohort earliest prostate recurrence matches a prostate_specific == True event in provenance
    rec_file = REPO_ROOT / "data_processed" / "tcga_prad_recurrence_cohort.csv"
    rec_df = pd.read_csv(rec_file)
    prostate_events = rec_df[rec_df["outcome_prostate_specific_recurrence"] == 1.0]

    for _, p_row in prostate_events.iterrows():
        pid = p_row["submitter_id"]
        pt_prov = prov_df[(prov_df["submitter_id"] == pid) & (prov_df["prostate_specific"].astype(str).str.lower() == "true")]
        assert len(pt_prov) > 0, f"Patient {pid} with positive prostate recurrence must have at least one prostate-specific provenance event"


def test_predictor_timing_audit_landmark():
    """Verify PREDICTOR_TIMING_AUDIT.csv correctly states landmark and avoids false surgery offsets."""
    timing_file = REPO_ROOT / "research" / "v3" / "PREDICTOR_TIMING_AUDIT.csv"
    assert timing_file.exists()
    timing_df = pd.read_csv(timing_file)

    # Ensure no false 'timing_relative_to_surgery_days = 0.0'
    assert "timing_relative_to_surgery_days" not in timing_df.columns, (
        "PREDICTOR_TIMING_AUDIT.csv must not declare timing_relative_to_surgery_days=0.0 when surgery date is unmeasured"
    )

    # Verify initial diagnosis day is 0.0 for baseline staging
    age_row = timing_df[timing_df["predictor"] == "age_at_diagnosis_years"].iloc[0]
    assert "0.0" in str(age_row["reference_origin_days"])

    gleason_row = timing_df[timing_df["predictor"] == "gleason_score_numeric"].iloc[0]
    assert "post-resection" in gleason_row["immortal_time_and_landmark_risk"].lower()
    assert "unmeasured" in gleason_row["reference_origin_days"].lower()
