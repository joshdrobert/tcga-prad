from __future__ import annotations

import re
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESEARCH_V2_DIR = PROJECT_ROOT / "research" / "v2"


def parse_numeric_tokens(val: object, reducer: Callable[[list[float]], float] = max) -> float:
    """Parse pipe-separated XML tokens into a single numeric value using a reducer (min/max)."""
    if pd.isna(val):
        return np.nan
    tokens = str(val).split("|")
    nums = []
    for token in tokens:
        cleaned = token.strip()
        try:
            nums.append(float(cleaned))
        except ValueError:
            continue
    return reducer(nums) if nums else np.nan


def categorize_bcr_status(val: object) -> str:
    """Classify Biochemical Recurrence status into explicit clinical categories."""
    if pd.isna(val):
        return "UNKNOWN_MISSING"
    text = str(val).upper().strip()
    tokens = [t.strip() for t in text.split("|")]
    has_pos = any(t in {"YES", "TRUE", "1", "BIOCHEMICAL"} for t in tokens)
    has_neg = any(t in {"NO", "FALSE", "0"} for t in tokens)
    if has_pos and has_neg:
        return "CONFLICTING"
    if has_pos:
        return "DOCUMENTED_POSITIVE"
    if has_neg:
        return "DOCUMENTED_NEGATIVE"
    return "UNKNOWN_MISSING"


def categorize_nte_status(val: object) -> str:
    """Classify New Tumor Event status into explicit clinical categories."""
    if pd.isna(val):
        return "UNKNOWN_MISSING"
    text = str(val).upper().strip()
    tokens = [t.strip() for t in text.split("|")]
    has_pos = any(t in {"YES", "TRUE", "1", "DISTANT", "METASTASIS", "PROGRESSION"} for t in tokens)
    has_neg = any(t in {"NO", "FALSE", "0"} for t in tokens)
    # Multi-encounter YES|NO or NO|YES reflects longitudinal event emergence
    if has_pos and has_neg:
        return "DOCUMENTED_POSITIVE_MULTI_ENCOUNTER"
    if has_pos:
        return "DOCUMENTED_POSITIVE"
    if has_neg:
        return "DOCUMENTED_NEGATIVE"
    return "UNKNOWN_MISSING"


def categorize_ht_progression(val: object) -> str:
    """Classify Post-Hormone Therapy progression status."""
    if pd.isna(val):
        return "NOT_APPLICABLE_OR_UNRECORDED"
    text = str(val).upper().strip()
    tokens = [t.strip() for t in text.split("|")]
    has_pos = any(t in {"YES", "TRUE", "1", "PROGRESSION"} for t in tokens)
    has_neg = any(t in {"NO", "FALSE", "0"} for t in tokens)
    if has_pos:
        return "DOCUMENTED_POSITIVE"
    if has_neg:
        return "DOCUMENTED_NEGATIVE"
    return "NOT_APPLICABLE_OR_UNRECORDED"


def main() -> None:
    RESEARCH_V2_DIR.mkdir(parents=True, exist_ok=True)
    base = pd.read_csv(PROJECT_ROOT / "data_processed" / "tcga_prad_analysis_cohort.csv")
    flat = pd.read_csv(PROJECT_ROOT / "data_processed" / "tcga_prad_clinical_xml_flat.csv")

    keep = [
        "submitter_id",
        "psa_value",
        "gleason_score",
        "biochemical_recurrence",
        "new_tumor_event_after_initial_treatment",
        "tumor_progression_post_ht",
        "days_to_first_biochemical_recurrence",
        "days_to_new_tumor_event_after_initial_treatment",
        "days_to_last_followup",
        "days_to_death",
        "supplement_kind",
        "source_file",
    ]
    # Prioritize clinical supplement XML over OMF (other malignancy form)
    flat_sorted = flat[[c for c in keep if c in flat.columns]].sort_values(
        by=["submitter_id", "supplement_kind"], ascending=[True, True]
    )
    clinical = flat_sorted.drop_duplicates("submitter_id", keep="first").copy()
    cohort = base.merge(clinical, on="submitter_id", how="left")

    # Clean longitudinal follow-up and event days
    cohort["days_to_last_followup_clean"] = cohort["days_to_last_followup_y"].apply(lambda v: parse_numeric_tokens(v, max))
    cohort["days_to_last_followup_clean"] = cohort["days_to_last_followup_clean"].fillna(
        cohort["days_to_last_followup_x"].apply(lambda v: parse_numeric_tokens(v, max))
    )
    # Check overall survival / days to death for non-recurrence censoring
    cohort["days_to_death_clean"] = cohort["days_to_death_y"].apply(lambda v: parse_numeric_tokens(v, max))
    cohort["days_to_death_clean"] = cohort["days_to_death_clean"].fillna(
        cohort["days_to_death_x"].apply(lambda v: parse_numeric_tokens(v, max))
    )

    cohort["days_to_bcr_clean"] = cohort["days_to_first_biochemical_recurrence"].apply(lambda v: parse_numeric_tokens(v, min))
    cohort["days_to_nte_clean"] = cohort["days_to_new_tumor_event_after_initial_treatment"].apply(lambda v: parse_numeric_tokens(v, min))

    # Explicit clinical status adjudication
    cohort["bcr_status"] = cohort["biochemical_recurrence"].apply(categorize_bcr_status)
    cohort["nte_status"] = cohort["new_tumor_event_after_initial_treatment"].apply(categorize_nte_status)
    cohort["ht_progression_status"] = cohort["tumor_progression_post_ht"].apply(categorize_ht_progression)

    # Boolean binary indicators (strictly distinguishing positive from negative vs missing)
    cohort["flag_bcr_positive"] = cohort["bcr_status"].eq("DOCUMENTED_POSITIVE").astype(int)
    cohort["flag_nte_positive"] = cohort["nte_status"].isin(["DOCUMENTED_POSITIVE", "DOCUMENTED_POSITIVE_MULTI_ENCOUNTER"]).astype(int)
    cohort["flag_ht_progression_positive"] = cohort["ht_progression_status"].eq("DOCUMENTED_POSITIVE").astype(int)

    # Composite Recurrence Categorization
    def adjudicate_composite(row: pd.Series) -> str:
        if row["flag_bcr_positive"] == 1 or row["flag_nte_positive"] == 1 or row["flag_ht_progression_positive"] == 1:
            return "DOCUMENTED_RECURRENCE"
        if row["bcr_status"] == "UNKNOWN_MISSING" and row["nte_status"] == "UNKNOWN_MISSING":
            return "INDETERMINATE_ALL_MISSING"
        if row["bcr_status"] == "DOCUMENTED_NEGATIVE" and row["nte_status"] == "DOCUMENTED_NEGATIVE":
            return "DOCUMENTED_NON_RECURRENCE"
        if row["bcr_status"] == "DOCUMENTED_NEGATIVE" and row["nte_status"] == "UNKNOWN_MISSING":
            return "NEGATIVE_BCR_MISSING_NTE"
        if row["bcr_status"] == "UNKNOWN_MISSING" and row["nte_status"] == "DOCUMENTED_NEGATIVE":
            return "MISSING_BCR_NEGATIVE_NTE"
        return "INDETERMINATE_OTHER"

    cohort["composite_status_detail"] = cohort.apply(adjudicate_composite, axis=1)

    # Binary outcome flag for modeling:
    # 1: Documented recurrence/progression (n=72)
    # 0: Documented recurrence-free (n=242 in full cohort, with sensitivity tracking for indeterminate cases)
    cohort["outcome_recurrence_or_progression"] = cohort["composite_status_detail"].eq("DOCUMENTED_RECURRENCE").astype(int)

    # Quality flag for binary analysis:
    # High confidence: documented positive (72) OR both BCR and NTE documented negative (191) = 263 patients
    cohort["high_confidence_outcome_flag"] = cohort["composite_status_detail"].isin(
        ["DOCUMENTED_RECURRENCE", "DOCUMENTED_NON_RECURRENCE"]
    ).astype(int)

    # Clean numeric clinical predictors
    cohort["psa_value_numeric"] = pd.to_numeric(cohort.get("psa_value"), errors="coerce")
    gleason_source = cohort.get("gleason_score_y", cohort.get("gleason_score_x", cohort.get("gleason_score")))
    cohort["gleason_score_numeric"] = pd.to_numeric(gleason_source, errors="coerce")

    # Longitudinal Event Timing Construction
    # NO FALLBACK TO LAST FOLLOWUP TIME FOR POSITIVE EVENTS WITHOUT DOCUMENTED DATES
    exact_event_times = []
    exact_censoring_times = []
    survival_analysis_times = []
    survival_event_flags = []
    undated_recurrence_flags = []

    for _, row in cohort.iterrows():
        is_event = row["outcome_recurrence_or_progression"] == 1
        t_event_candidates = []
        if pd.notna(row["days_to_bcr_clean"]) and row["days_to_bcr_clean"] > 0:
            t_event_candidates.append(row["days_to_bcr_clean"])
        if pd.notna(row["days_to_nte_clean"]) and row["days_to_nte_clean"] > 0:
            t_event_candidates.append(row["days_to_nte_clean"])

        t_censor_candidates = []
        if pd.notna(row["days_to_last_followup_clean"]) and row["days_to_last_followup_clean"] > 0:
            t_censor_candidates.append(row["days_to_last_followup_clean"])
        elif pd.notna(row["days_to_death_clean"]) and row["days_to_death_clean"] > 0:
            t_censor_candidates.append(row["days_to_death_clean"])

        t_event = min(t_event_candidates) if t_event_candidates else np.nan
        t_censor = max(t_censor_candidates) if t_censor_candidates else np.nan

        exact_event_times.append(t_event)
        exact_censoring_times.append(t_censor)

        if is_event:
            if pd.notna(t_event):
                # Valid dated event
                survival_analysis_times.append(t_event)
                survival_event_flags.append(1)
                undated_recurrence_flags.append(0)
            else:
                # Documented recurrence but NO exact date recorded:
                # DO NOT ASSIGN LAST FOLLOWUP TIME AS EVENT DATE!
                # Mark as NaN in exact-time survival analysis to prevent falsification.
                survival_analysis_times.append(np.nan)
                survival_event_flags.append(np.nan)
                undated_recurrence_flags.append(1)
        else:
            if pd.notna(t_censor):
                survival_analysis_times.append(t_censor)
                survival_event_flags.append(0)
                undated_recurrence_flags.append(0)
            else:
                survival_analysis_times.append(np.nan)
                survival_event_flags.append(np.nan)
                undated_recurrence_flags.append(0)

    cohort["days_to_recurrence_exact"] = exact_event_times
    cohort["days_to_censoring_exact"] = exact_censoring_times
    cohort["survival_duration_days"] = survival_analysis_times
    cohort["survival_event_indicator"] = survival_event_flags
    cohort["undated_recurrence_flag"] = undated_recurrence_flags

    # Historical compatibility column name (only populated for valid dated cases or non-events)
    cohort["time_to_recurrence_days"] = cohort["survival_duration_days"]

    # Save processed recurrence cohort
    out = PROJECT_ROOT / "data_processed" / "tcga_prad_recurrence_cohort.csv"
    cohort.to_csv(out, index=False)

    # Save V2 Longitudinal Event Audit CSV
    audit_cols = [
        "submitter_id",
        "composite_status_detail",
        "outcome_recurrence_or_progression",
        "bcr_status",
        "nte_status",
        "ht_progression_status",
        "days_to_bcr_clean",
        "days_to_nte_clean",
        "days_to_last_followup_clean",
        "days_to_death_clean",
        "days_to_recurrence_exact",
        "days_to_censoring_exact",
        "survival_duration_days",
        "survival_event_indicator",
        "undated_recurrence_flag",
    ]
    cohort[audit_cols].to_csv(RESEARCH_V2_DIR / "LONGITUDINAL_EVENT_AUDIT.csv", index=False)

    # Audit the 6 undated recurrence cases specifically
    undated = cohort[cohort["undated_recurrence_flag"] == 1][audit_cols]
    undated.to_csv(RESEARCH_V2_DIR / "UNDATED_RECURRENCE_CASES.csv", index=False)

    # Save V2 Cohort Flow Breakdown
    flow = pd.DataFrame(
        [
            {"step": "Total TCGA-PRAD Radical Prostatectomy Patients", "count": len(cohort)},
            {"step": "Documented Recurrence Events (Composite)", "count": int(cohort["outcome_recurrence_or_progression"].sum())},
            {"step": "  - Dated Recurrence Events (Exact Timing Available)", "count": int((cohort["survival_event_indicator"] == 1).sum())},
            {"step": "  - Undated Recurrence Events (Timing Missing in XML)", "count": int(cohort["undated_recurrence_flag"].sum())},
            {"step": "Documented Recurrence-Free (Non-Events)", "count": int((cohort["outcome_recurrence_or_progression"] == 0).sum())},
            {"step": "  - High Confidence (BCR=NO & NTE=NO)", "count": int((cohort["composite_status_detail"] == "DOCUMENTED_NON_RECURRENCE").sum())},
            {"step": "  - Asymmetric Negative (BCR=NO, NTE missing)", "count": int((cohort["composite_status_detail"] == "NEGATIVE_BCR_MISSING_NTE").sum())},
            {"step": "  - Asymmetric Negative (NTE=NO, BCR missing)", "count": int((cohort["composite_status_detail"] == "MISSING_BCR_NEGATIVE_NTE").sum())},
            {"step": "  - Indeterminate (BCR missing & NTE missing)", "count": int((cohort["composite_status_detail"] == "INDETERMINATE_ALL_MISSING").sum())},
            {"step": "Evaluable Continuous-Time Survival Cohort (Dated events + Follow-up)", "count": int(cohort["survival_duration_days"].notna().sum())},
        ]
    )
    flow.to_csv(RESEARCH_V2_DIR / "COHORT_FLOW.csv", index=False)

    # Missingness Table
    missing_tbl = pd.DataFrame(
        [
            {"field": "psa_value", "missing_count": int(cohort["psa_value_numeric"].isna().sum()), "pct": float(cohort["psa_value_numeric"].isna().mean() * 100)},
            {"field": "gleason_score", "missing_count": int(cohort["gleason_score_numeric"].isna().sum()), "pct": float(cohort["gleason_score_numeric"].isna().mean() * 100)},
            {"field": "biochemical_recurrence", "missing_count": int(cohort["biochemical_recurrence"].isna().sum()), "pct": float(cohort["biochemical_recurrence"].isna().mean() * 100)},
            {"field": "new_tumor_event", "missing_count": int(cohort["new_tumor_event_after_initial_treatment"].isna().sum()), "pct": float(cohort["new_tumor_event_after_initial_treatment"].isna().mean() * 100)},
            {"field": "days_to_last_followup", "missing_count": int(cohort["days_to_last_followup_clean"].isna().sum()), "pct": float(cohort["days_to_last_followup_clean"].isna().mean() * 100)},
            {"field": "event_recurrence_exact_date", "missing_count": int(cohort[cohort['outcome_recurrence_or_progression'] == 1]['days_to_recurrence_exact'].isna().sum()), "pct": float(cohort[cohort['outcome_recurrence_or_progression'] == 1]['days_to_recurrence_exact'].isna().mean() * 100)},
        ]
    )
    missing_tbl.to_csv(RESEARCH_V2_DIR / "MISSINGNESS_TABLE.csv", index=False)

    # Legacy table for backward compatibility
    counts = pd.DataFrame(
        [
            {"endpoint": "Composite Recurrence/Progression (Total)", "events": int(cohort["outcome_recurrence_or_progression"].sum()), "total": len(cohort)},
            {"endpoint": "Dated Recurrence Events (Exact Timing)", "events": int((cohort["survival_event_indicator"] == 1).sum()), "total": len(cohort)},
            {"endpoint": "Undated Recurrence Events (Excluded from exact time)", "events": int(cohort["undated_recurrence_flag"].sum()), "total": len(cohort)},
            {"endpoint": "Biochemical Recurrence (BCR)", "events": int(cohort["flag_bcr_positive"].sum()), "total": len(cohort)},
            {"endpoint": "New Tumor Event (NTE)", "events": int(cohort["flag_nte_positive"].sum()), "total": len(cohort)},
            {"endpoint": "Post-HT Progression", "events": int(cohort["flag_ht_progression_positive"].sum()), "total": len(cohort)},
            {"endpoint": "Evaluable RFS (with valid timing)", "events": int((cohort["survival_event_indicator"] == 1).sum()), "total": int(cohort["survival_duration_days"].notna().sum())},
        ]
    )
    counts.to_csv(PROJECT_ROOT / "tables" / "table_recurrence_outcome_counts.csv", index=False)

    print(f"=== TCGA-PRAD V2 ENDPOINT RECONSTRUCTION COMPLETE ===")
    print(flow.to_string(index=False))


if __name__ == "__main__":
    main()
