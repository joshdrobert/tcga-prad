from __future__ import annotations

import json
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

from src.build_recurrence_endpoint import parse_numeric_tokens, parse_boolean_flag


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


def test_parse_boolean_flag_tokens():
    assert parse_boolean_flag("YES") == 1
    assert parse_boolean_flag("NO") == 0
    assert parse_boolean_flag("NO|YES") == 1
    assert parse_boolean_flag("YES|NO") == 1
    assert parse_boolean_flag("BIOCHEMICAL") == 1
    assert parse_boolean_flag(np.nan) == 0
    assert parse_boolean_flag("[Unknown]") == 0


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

    # Timing check
    assert "time_to_recurrence_days" in df.columns
    valid_timing = df["time_to_recurrence_days"].dropna()
    assert len(valid_timing) >= 300, "At least 300 patients must have non-null follow-up/event days"
    assert (valid_timing > 0).all(), "All event/censor times must be positive"


def test_cbioportal_concordance_invariants():
    repo_root = Path(__file__).resolve().parents[1]
    crosscheck_file = repo_root / "tables" / "table_cbioportal_endpoint_crosscheck.csv"
    assert crosscheck_file.exists()
    cc = pd.read_csv(crosscheck_file)
    row = cc.iloc[0]
    assert row["matched_patients"] == 269
    assert row["cbio_bcr_events"] == 38
    # Exact discordances = 31
    evaluable = 269
    agreement = row["agreement_fraction"]
    concordant = round(evaluable * agreement)
    discordant = evaluable - concordant
    assert discordant == 31, f"Expected exactly 31 discordant cases, got {discordant}"


def test_metrics_json_integrity():
    repo_root = Path(__file__).resolve().parents[1]
    artifact = repo_root / "research" / "artifacts" / "tcga_prad_metrics.json"
    assert artifact.exists()
    with open(artifact) as f:
        data = json.load(f)
    assert data["cohort"]["total_patients"] == 314
    assert data["cohort"]["rna_seq_profiles_available"] == 228
    assert data["cohort"]["rna_seq_missing_count"] == 86
    assert data["survival"]["cox_concordance_index"] > 0.65
