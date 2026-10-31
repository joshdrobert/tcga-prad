from __future__ import annotations

import json
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

from src.build_recurrence_endpoint import (
    parse_numeric_tokens,
    categorize_bcr_status,
    categorize_nte_status,
    categorize_ht_progression,
)


def test_parse_numeric_tokens_single():
    assert parse_numeric_tokens("120") == 120.0
    assert parse_numeric_tokens("45.5") == 45.5


def test_parse_numeric_tokens_multi_pipe():
    # Multi-encounter longitudinal follow-up tokens
    val = "120|365|730"
    assert parse_numeric_tokens(val, reducer=max) == 730.0
    assert parse_numeric_tokens(val, reducer=min) == 120.0


def test_parse_numeric_tokens_invalid_and_missing():
    assert np.isnan(parse_numeric_tokens(np.nan))
    assert np.isnan(parse_numeric_tokens("[Not Available]"))
    assert np.isnan(parse_numeric_tokens(""))
    # Mixed valid and invalid tokens
    assert parse_numeric_tokens("[Not Available]|450|None", reducer=max) == 450.0


def test_endpoint_categorization_semantics():
    # BCR
    assert categorize_bcr_status("YES") == "DOCUMENTED_POSITIVE"
    assert categorize_bcr_status("NO") == "DOCUMENTED_NEGATIVE"
    assert categorize_bcr_status(np.nan) == "UNKNOWN_MISSING"
    assert categorize_bcr_status("[Unknown]") == "UNKNOWN_MISSING"

    # NTE
    assert categorize_nte_status("YES") == "DOCUMENTED_POSITIVE"
    assert categorize_nte_status("NO|YES") == "DOCUMENTED_POSITIVE_MULTI_ENCOUNTER"
    assert categorize_nte_status("YES|NO") == "DOCUMENTED_POSITIVE_MULTI_ENCOUNTER"
    assert categorize_nte_status("NO") == "DOCUMENTED_NEGATIVE"
    assert categorize_nte_status(np.nan) == "UNKNOWN_MISSING"

    # HT progression
    assert categorize_ht_progression("YES") == "DOCUMENTED_POSITIVE"
    assert categorize_ht_progression("NO") == "DOCUMENTED_NEGATIVE"
    assert categorize_ht_progression(np.nan) == "NOT_APPLICABLE_OR_UNRECORDED"


def test_processed_cohort_invariants():
    repo_root = Path(__file__).resolve().parents[1]
    rec_file = repo_root / "data_processed" / "tcga_prad_recurrence_cohort.csv"
    assert rec_file.exists(), "Processed recurrence cohort must exist"

    df = pd.read_csv(rec_file)
    assert len(df) == 314, f"Expected 314 patients, found {len(df)}"
    assert df["submitter_id"].is_unique, "Patient IDs must be strictly unique"

    # Outcome check
    events = int(df["outcome_recurrence_or_progression"].sum())
    assert events == 72, f"Expected 72 recurrence events, found {events}"

    # Exact timing check
    assert "survival_duration_days" in df.columns
    assert "survival_event_indicator" in df.columns
    assert "undated_recurrence_flag" in df.columns

    # Verify exactly 6 undated events
    undated_count = int(df["undated_recurrence_flag"].sum())
    assert undated_count == 6, f"Expected exactly 6 undated recurrence cases, got {undated_count}"

    # Continuous-time evaluable cohort must be 308
    evaluable = df[df["survival_duration_days"].notna() & (df["survival_duration_days"] > 0)]
    assert len(evaluable) == 308, f"Expected 308 evaluable survival cases, got {len(evaluable)}"
    assert int(evaluable["survival_event_indicator"].sum()) == 66, "Expected 66 dated recurrence events"


def test_cbioportal_concordance_invariants():
    repo_root = Path(__file__).resolve().parents[1]
    crosscheck_file = repo_root / "tables" / "table_cbioportal_endpoint_crosscheck.csv"
    assert crosscheck_file.exists()
    cc = pd.read_csv(crosscheck_file)
    row = cc.iloc[0]
    assert row["evaluable_matched"] == 269
    assert row["cbio_bcr_events"] == 38
    assert row["discordant_cases"] == 31
    assert row["local_only_events"] == 31
    assert row["cbio_only_events"] == 0
    assert abs(row["agreement_fraction"] - 0.884758) < 1e-4


def test_metrics_json_integrity():
    repo_root = Path(__file__).resolve().parents[1]
    artifact = repo_root / "research" / "artifacts" / "tcga_prad_metrics.json"
    assert artifact.exists()
    with open(artifact) as f:
        data = json.load(f)
    assert data["cohort"]["total_patients"] == 314
    assert data["cohort"]["composite_recurrence_events"] == 72
    assert data["cohort"]["dated_recurrence_events"] == 66
    assert data["cohort"]["undated_recurrence_events"] == 6
    assert data["cohort"]["rna_seq_profiles_available"] == 228
    assert data["cohort"]["rna_seq_missing_count"] == 86
    assert data["survival"]["evaluable_cases"] if "evaluable_cases" in data["survival"] else data["cohort"]["evaluable_survival_cases"] == 308
    assert data["survival"]["cox_apparent_concordance_index"] > 0.65
    assert data["survival"]["cox_cv_concordance_index"] > 0.65
