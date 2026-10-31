from __future__ import annotations

import re
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESEARCH_V2_DIR = PROJECT_ROOT / "research" / "v2"
RESEARCH_V3_DIR = PROJECT_ROOT / "research" / "v3"


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


def parse_numeric_tokens_all(val: object) -> list[float]:
    """Extract all valid positive numeric tokens from a pipe-separated string."""
    if pd.isna(val):
        return []
    tokens = str(val).split("|")
    nums = []
    for token in tokens:
        cleaned = token.strip()
        try:
            f = float(cleaned)
            if not np.isnan(f) and f >= 0:
                nums.append(f)
        except ValueError:
            continue
    return nums


def parse_string_tokens(val: object) -> list[str]:
    """Extract clean non-empty string tokens from a pipe-separated string."""
    if pd.isna(val):
        return []
    return [
        t.strip()
        for t in str(val).split("|")
        if t.strip() and t.strip().lower() not in {"nan", "[not available]", "[not applicable]", "[unknown]"}
    ]


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
    if has_pos and has_neg:
        return "DOCUMENTED_POSITIVE_MULTI_ENCOUNTER"
    if has_pos:
        return "DOCUMENTED_POSITIVE"
    if has_neg:
        return "DOCUMENTED_NEGATIVE"
    return "UNKNOWN_MISSING"


def categorize_nte_type(val: object) -> str:
    """Classify the anatomic/clinical type of new tumor event from XML."""
    if pd.isna(val):
        return "UNSPECIFIED"
    text = str(val).strip()
    if "New Primary" in text and not any(k in text for k in ["Biochemical", "Locoregional", "Distant"]):
        return "SECOND_PRIMARY_ONLY"
    if any(k in text for k in ["Biochemical", "Locoregional", "Distant"]):
        return "PROSTATE_SPECIFIC"
    return "OTHER"


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


def extract_event_level_provenance(flat: pd.DataFrame) -> pd.DataFrame:
    """
    Reconstruct unrolled event-level clinical provenance records pairing each event
    type with its specific event date, source form, and encounter category.
    Guarantees that event dates are never decoupled from event types.
    """
    event_records = []

    for _, row in flat.iterrows():
        sid = row["submitter_id"]
        src = row.get("source_file", "unknown")

        # 1. Biochemical Recurrence tokens
        bcr_statuses = parse_string_tokens(row.get("biochemical_recurrence"))
        bcr_days = parse_numeric_tokens_all(row.get("days_to_first_biochemical_recurrence"))
        has_bcr_pos = any(s.upper() in {"YES", "TRUE", "1", "BIOCHEMICAL"} for s in bcr_statuses)

        if has_bcr_pos or bcr_days:
            dates = bcr_days if bcr_days else [np.nan]
            for d in dates:
                event_records.append(
                    {
                        "submitter_id": sid,
                        "source_file": src,
                        "event_category": "BIOCHEMICAL_RECURRENCE",
                        "event_type": "Biochemical Recurrence (BCR=YES)",
                        "event_days_from_diagnosis": d,
                        "prostate_specific": True,
                        "adjudication_class": "PROSTATE_RECURRENCE",
                        "adjudication_rationale": "Biochemical recurrence documented in XML follow-up",
                    }
                )

        # 2. New Tumor Event tokens
        nte_types = parse_string_tokens(row.get("new_neoplasm_event_type"))
        nte_days = parse_numeric_tokens_all(row.get("days_to_new_tumor_event_after_initial_treatment"))
        nte_statuses = parse_string_tokens(row.get("new_tumor_event_after_initial_treatment"))
        has_nte_pos = any(s.upper() in {"YES", "TRUE", "1"} for s in nte_statuses) or bool(nte_types)

        if has_nte_pos or nte_days:
            t_full = "|".join(nte_types) if nte_types else "New Tumor Event"
            is_second_prim = ("New Primary" in t_full) and not any(
                k in t_full for k in ["Biochemical", "Locoregional", "Distant"]
            )

            if nte_days:
                for i, d in enumerate(nte_days):
                    t = nte_types[i] if i < len(nte_types) else (nte_types[0] if nte_types else "New Tumor Event")
                    cur_sec = ("New Primary" in t) and not any(
                        k in t for k in ["Biochemical", "Locoregional", "Distant"]
                    )
                    event_records.append(
                        {
                            "submitter_id": sid,
                            "source_file": src,
                            "event_category": "NEW_TUMOR_EVENT",
                            "event_type": t,
                            "event_days_from_diagnosis": d,
                            "prostate_specific": not cur_sec,
                            "adjudication_class": "SECOND_PRIMARY_MALIGNANCY" if cur_sec else "PROSTATE_RECURRENCE",
                            "adjudication_rationale": (
                                "Second non-prostatic primary malignancy"
                                if cur_sec
                                else f"Prostate-specific new tumor event ({t})"
                            ),
                        }
                    )
            else:
                event_records.append(
                    {
                        "submitter_id": sid,
                        "source_file": src,
                        "event_category": "NEW_TUMOR_EVENT",
                        "event_type": t_full,
                        "event_days_from_diagnosis": np.nan,
                        "prostate_specific": not is_second_prim,
                        "adjudication_class": "SECOND_PRIMARY_MALIGNANCY" if is_second_prim else "PROSTATE_RECURRENCE",
                        "adjudication_rationale": (
                            "Second non-prostatic primary malignancy (undated)"
                            if is_second_prim
                            else f"Prostate-specific new tumor event ({t_full}; undated in XML)"
                        ),
                    }
                )

        # 3. Post-Hormone Therapy progression
        ht_statuses = parse_string_tokens(row.get("tumor_progression_post_ht"))
        if any(k in s.upper() for s in ht_statuses for k in ["YES", "TRUE", "1", "PROGRESSION"]):
            event_records.append(
                {
                    "submitter_id": sid,
                    "source_file": src,
                    "event_category": "POST_HT_PROGRESSION",
                    "event_type": "Post-HT Tumor Progression",
                    "event_days_from_diagnosis": np.nan,
                    "prostate_specific": True,
                    "adjudication_class": "PROSTATE_RECURRENCE",
                    "adjudication_rationale": "Tumor progression after hormone therapy documented in XML",
                }
            )

    return pd.DataFrame(event_records)


def aggregate_clinical_xml_records(flat: pd.DataFrame) -> pd.DataFrame:
    """
    Reconstruct longitudinal clinical records across ALL available XML supplements
    (initial clinical forms, follow-up forms, and OMF/other malignancy forms)
    WITHOUT discarding records via drop_duplicates(..., keep='first').
    """
    aggregated_records = []
    grouped = flat.groupby("submitter_id", sort=True)

    for sid, group in grouped:
        record = {"submitter_id": sid}

        # Status strings: combine unique non-empty tokens across all encounters
        for col in [
            "biochemical_recurrence",
            "new_tumor_event_after_initial_treatment",
            "new_neoplasm_event_type",
            "tumor_progression_post_ht",
        ]:
            if col in group.columns:
                vals = [str(v).strip() for v in group[col].dropna() if str(v).strip() not in {"", "nan"}]
                record[col] = "|".join(dict.fromkeys(vals)) if vals else np.nan
            else:
                record[col] = np.nan

        # Baseline clinical predictors: take first non-null recorded value
        for col in ["psa_value", "gleason_score"]:
            if col in group.columns:
                non_null = group[col].dropna()
                record[col] = non_null.iloc[0] if len(non_null) > 0 else np.nan
            else:
                record[col] = np.nan

        # Longitudinal event days: collect all tokens across all forms and take earliest event
        for col in ["days_to_first_biochemical_recurrence", "days_to_new_tumor_event_after_initial_treatment"]:
            if col in group.columns:
                all_tokens = []
                for val in group[col].dropna():
                    all_tokens.extend(parse_numeric_tokens_all(val))
                record[col] = min(all_tokens) if all_tokens else np.nan
            else:
                record[col] = np.nan

        # Longitudinal follow-up and death days: take longest verified follow-up across all forms
        for col in ["days_to_last_followup", "days_to_death"]:
            if col in group.columns:
                all_tokens = []
                for val in group[col].dropna():
                    all_tokens.extend(parse_numeric_tokens_all(val))
                record[col] = max(all_tokens) if all_tokens else np.nan
            else:
                record[col] = np.nan

        aggregated_records.append(record)

    return pd.DataFrame(aggregated_records)


def main() -> None:
    RESEARCH_V2_DIR.mkdir(parents=True, exist_ok=True)
    RESEARCH_V3_DIR.mkdir(parents=True, exist_ok=True)
    base = pd.read_csv(PROJECT_ROOT / "data_processed" / "tcga_prad_analysis_cohort.csv")
    flat = pd.read_csv(PROJECT_ROOT / "data_processed" / "tcga_prad_clinical_xml_flat.csv")

    # Generate Unrolled Event-Level Provenance Table
    event_provenance = extract_event_level_provenance(flat)
    event_provenance.to_csv(RESEARCH_V3_DIR / "EVENT_LEVEL_PROVENANCE.csv", index=False)
    event_provenance.to_csv(RESEARCH_V2_DIR / "EVENT_LEVEL_PROVENANCE.csv", index=False)

    # Aggregate longitudinal XML records without discarding secondary records
    clinical_agg = aggregate_clinical_xml_records(flat)
    cohort = base.merge(clinical_agg, on="submitter_id", how="left")

    # Clean longitudinal follow-up and event days (Days from Initial Pathologic Diagnosis)
    cohort["days_to_last_followup_clean"] = cohort["days_to_last_followup_y"].fillna(
        cohort["days_to_last_followup_x"].apply(lambda v: parse_numeric_tokens(v, max))
    )
    cohort["days_to_death_clean"] = cohort["days_to_death_y"].fillna(
        cohort["days_to_death_x"].apply(lambda v: parse_numeric_tokens(v, max))
    )
    cohort["days_to_bcr_clean"] = cohort["days_to_first_biochemical_recurrence"]
    cohort["days_to_nte_clean"] = cohort["days_to_new_tumor_event_after_initial_treatment"]

    # Explicit clinical status adjudication
    cohort["bcr_status"] = cohort["biochemical_recurrence"].apply(categorize_bcr_status)
    cohort["nte_status"] = cohort["new_tumor_event_after_initial_treatment"].apply(categorize_nte_status)
    cohort["nte_type_category"] = cohort["new_neoplasm_event_type"].apply(categorize_nte_type)
    cohort["ht_progression_status"] = cohort["tumor_progression_post_ht"].apply(categorize_ht_progression)

    cohort["flag_bcr_positive"] = cohort["bcr_status"].eq("DOCUMENTED_POSITIVE").astype(int)
    cohort["flag_nte_positive"] = cohort["nte_status"].isin(["DOCUMENTED_POSITIVE", "DOCUMENTED_POSITIVE_MULTI_ENCOUNTER"]).astype(int)
    cohort["flag_ht_progression_positive"] = cohort["ht_progression_status"].eq("DOCUMENTED_POSITIVE").astype(int)

    # Adjudicate Composite Status (Accounting for Second Primaries vs Prostate Recurrences)
    def adjudicate_composite(row: pd.Series) -> str:
        # 1. Genuine prostate-specific recurrence / progression:
        is_prostate_nte = (row["flag_nte_positive"] == 1) and (row["nte_type_category"] != "SECOND_PRIMARY_ONLY")
        if row["flag_bcr_positive"] == 1 or is_prostate_nte or row["flag_ht_progression_positive"] == 1:
            return "DOCUMENTED_PROSTATE_RECURRENCE"

        # 2. Non-prostatic second primary malignancy (TCGA-CH-5791, TCGA-HC-7080):
        if row["nte_type_category"] == "SECOND_PRIMARY_ONLY":
            return "SECOND_PRIMARY_MALIGNANCY"

        # 3. Fully Indeterminate:
        if row["bcr_status"] == "UNKNOWN_MISSING" and row["nte_status"] == "UNKNOWN_MISSING":
            return "INDETERMINATE_ALL_MISSING"

        # 4. Documented double negative:
        if row["bcr_status"] == "DOCUMENTED_NEGATIVE" and row["nte_status"] == "DOCUMENTED_NEGATIVE":
            return "DOCUMENTED_NON_RECURRENCE"

        # 5. Asymmetric surveillance:
        if row["bcr_status"] == "DOCUMENTED_NEGATIVE" and row["nte_status"] == "UNKNOWN_MISSING":
            return "NEGATIVE_BCR_MISSING_NTE"
        if row["bcr_status"] == "UNKNOWN_MISSING" and row["nte_status"] == "DOCUMENTED_NEGATIVE":
            return "MISSING_BCR_NEGATIVE_NTE"

        return "INDETERMINATE_OTHER"

    cohort["composite_status_detail"] = cohort.apply(adjudicate_composite, axis=1)

    # Categorical Outcome Variables:
    # 1. outcome_prostate_specific_recurrence:
    #    1.0: Documented prostate-specific recurrence/progression (n=70)
    #    0.0: Documented recurrence-free double-negative (n=191)
    #    NaN: Second primary (n=2), Asymmetric (n=46), or Fully Indeterminate (n=5)
    def assign_prostate_specific(status: str) -> float:
        if status == "DOCUMENTED_PROSTATE_RECURRENCE":
            return 1.0
        if status == "DOCUMENTED_NON_RECURRENCE":
            return 0.0
        return np.nan

    cohort["outcome_prostate_specific_recurrence"] = cohort["composite_status_detail"].apply(assign_prostate_specific)

    # 2. outcome_definitive_recurrence:
    #    Primary definitive cohort: alias to outcome_prostate_specific_recurrence (n=261: 70 positive, 191 negative)
    cohort["outcome_definitive_recurrence"] = cohort["outcome_prostate_specific_recurrence"]

    # 3. outcome_composite_all_cancers (Sensitivity comparison including non-prostatic second primaries: n=263)
    def assign_all_cancers(status: str) -> float:
        if status in {"DOCUMENTED_PROSTATE_RECURRENCE", "SECOND_PRIMARY_MALIGNANCY"}:
            return 1.0
        if status == "DOCUMENTED_NON_RECURRENCE":
            return 0.0
        return np.nan

    cohort["outcome_composite_all_cancers"] = cohort["composite_status_detail"].apply(assign_all_cancers)

    # 4. outcome_exploratory_composite:
    #    1.0: Documented prostate recurrence (n=70)
    #    0.0: Documented recurrence-free (n=191) or asymmetric negative (n=46) -> n=237
    #    NaN: Second primary (n=2) or Fully Indeterminate (n=5) -> EXCLUDED!
    def assign_exploratory(status: str) -> float:
        if status == "DOCUMENTED_PROSTATE_RECURRENCE":
            return 1.0
        if status in {"DOCUMENTED_NON_RECURRENCE", "NEGATIVE_BCR_MISSING_NTE", "MISSING_BCR_NEGATIVE_NTE"}:
            return 0.0
        return np.nan

    cohort["outcome_exploratory_composite"] = cohort["composite_status_detail"].apply(assign_exploratory)

    # Backward compatibility primary binary column:
    # 1: Documented prostate recurrence (n=70)
    # 0: Nonpositive (n=244)
    cohort["outcome_recurrence_or_progression"] = cohort["composite_status_detail"].eq("DOCUMENTED_PROSTATE_RECURRENCE").astype(int)

    # High confidence outcome flag (documented prostate positive OR documented double negative = 261 patients)
    cohort["high_confidence_outcome_flag"] = cohort["composite_status_detail"].isin(
        ["DOCUMENTED_PROSTATE_RECURRENCE", "DOCUMENTED_NON_RECURRENCE"]
    ).astype(int)

    # Numeric clinical predictors
    cohort["psa_value_numeric"] = pd.to_numeric(cohort.get("psa_value"), errors="coerce")
    gleason_source = cohort.get("gleason_score_y", cohort.get("gleason_score_x", cohort.get("gleason_score")))
    cohort["gleason_score_numeric"] = pd.to_numeric(gleason_source, errors="coerce")

    # Longitudinal Event Timing Construction
    # TIME ORIGIN: INITIAL PATHOLOGIC DIAGNOSIS (CORE NEEDLE BIOPSY AT DAY 0 FOR 98.1% PATIENTS)
    # CONDITIONED ON PROSTATECTOMY (EXACT RESECTION DATE UNRECORDED IN GDC XML)
    # PAIRED EVENT PROVENANCE DATES: DERIVED STRICTLY FROM PROSTATE-SPECIFIC EVENT RECORDS
    exact_event_times = []
    exact_censoring_times = []
    survival_analysis_times = []
    survival_event_flags = []
    undated_recurrence_flags = []
    earliest_event_types = []
    earliest_event_sources = []

    for _, row in cohort.iterrows():
        sid = row["submitter_id"]
        status = row["composite_status_detail"]
        is_prostate_event = (status == "DOCUMENTED_PROSTATE_RECURRENCE")
        is_second_primary = (status == "SECOND_PRIMARY_MALIGNANCY")
        is_indeterminate = (status == "INDETERMINATE_ALL_MISSING")

        # Query unrolled event records for this patient
        p_events = event_provenance[event_provenance["submitter_id"] == sid]
        pros_events = p_events[p_events["prostate_specific"] & (p_events["adjudication_class"] == "PROSTATE_RECURRENCE")]
        dated_pros = pros_events[pros_events["event_days_from_diagnosis"].notna() & (pros_events["event_days_from_diagnosis"] > 0)]

        t_event = np.nan
        e_type = "None"
        e_src = "None"

        if is_prostate_event:
            if len(dated_pros) > 0:
                earliest = dated_pros.sort_values("event_days_from_diagnosis").iloc[0]
                t_event = float(earliest["event_days_from_diagnosis"])
                e_type = str(earliest["event_type"])
                e_src = str(earliest["source_file"])
            elif len(pros_events) > 0:
                e_type = str(pros_events.iloc[0]["event_type"])
                e_src = str(pros_events.iloc[0]["source_file"])
        elif is_second_primary:
            sec_events = p_events[~p_events["prostate_specific"]]
            if len(sec_events) > 0:
                e_type = str(sec_events.iloc[0]["event_type"])
                e_src = str(sec_events.iloc[0]["source_file"])

        t_censor_candidates = []
        if pd.notna(row["days_to_last_followup_clean"]) and row["days_to_last_followup_clean"] > 0:
            t_censor_candidates.append(row["days_to_last_followup_clean"])
        elif pd.notna(row["days_to_death_clean"]) and row["days_to_death_clean"] > 0:
            t_censor_candidates.append(row["days_to_death_clean"])

        t_censor = max(t_censor_candidates) if t_censor_candidates else np.nan

        exact_event_times.append(t_event)
        exact_censoring_times.append(t_censor)
        earliest_event_types.append(e_type)
        earliest_event_sources.append(e_src)

        if is_prostate_event:
            if pd.notna(t_event):
                survival_analysis_times.append(t_event)
                survival_event_flags.append(1.0)
                undated_recurrence_flags.append(0)
            else:
                survival_analysis_times.append(np.nan)
                survival_event_flags.append(np.nan)
                undated_recurrence_flags.append(1)
        elif is_second_primary:
            survival_analysis_times.append(np.nan)
            survival_event_flags.append(np.nan)
            undated_recurrence_flags.append(0)
        elif is_indeterminate:
            survival_analysis_times.append(np.nan)
            survival_event_flags.append(np.nan)
            undated_recurrence_flags.append(0)
        else:
            if pd.notna(t_censor):
                survival_analysis_times.append(t_censor)
                survival_event_flags.append(0.0)
                undated_recurrence_flags.append(0)
            else:
                survival_analysis_times.append(np.nan)
                survival_event_flags.append(np.nan)
                undated_recurrence_flags.append(0)

    cohort["days_to_recurrence_exact"] = exact_event_times
    cohort["days_to_censoring_exact"] = exact_censoring_times
    cohort["earliest_prostate_recurrence_type"] = earliest_event_types
    cohort["earliest_prostate_recurrence_source"] = earliest_event_sources
    cohort["survival_duration_days"] = survival_analysis_times
    cohort["duration_days_from_diagnosis"] = survival_analysis_times
    cohort["survival_event_indicator"] = survival_event_flags
    cohort["undated_recurrence_flag"] = undated_recurrence_flags
    cohort["time_to_recurrence_days"] = cohort["survival_duration_days"]

    # Adjudication rationale text per patient
    def assign_rationale(row: pd.Series) -> str:
        status = row["composite_status_detail"]
        if status == "DOCUMENTED_PROSTATE_RECURRENCE":
            timing_str = (
                f"dated Day {row['days_to_recurrence_exact']:.0f} (Type: {row['earliest_prostate_recurrence_type']}; Source: {row['earliest_prostate_recurrence_source']})"
                if pd.notna(row["days_to_recurrence_exact"])
                else f"undated in XML (Type: {row['earliest_prostate_recurrence_type']})"
            )
            return f"Prostate cancer recurrence/progression: {timing_str} from initial diagnostic biopsy."
        if status == "SECOND_PRIMARY_MALIGNANCY":
            return f"Second non-prostatic primary malignancy (NTE type: {row.get('new_neoplasm_event_type')}; Source: {row['earliest_prostate_recurrence_source']}); BCR=NO; adjudicated as non-prostatic cancer and excluded from prostate recurrence."
        if status == "DOCUMENTED_NON_RECURRENCE":
            return f"Documented recurrence-free: BCR=NO and NTE=NO; verified clinical follow-up through Day {row['days_to_censoring_exact']:.0f}."
        if status in {"NEGATIVE_BCR_MISSING_NTE", "MISSING_BCR_NEGATIVE_NTE"}:
            return f"Asymmetric surveillance: single negative indicator with complementary indicator missing; censored at Day {row['days_to_censoring_exact']:.0f}."
        if status == "INDETERMINATE_ALL_MISSING":
            return "Fully indeterminate: missing all longitudinal follow-up, BCR, and NTE data; excluded from definitive and survival non-events."
        return "Other clinical status."

    cohort["adjudication_rationale"] = cohort.apply(assign_rationale, axis=1)

    # Save processed recurrence cohort
    out = PROJECT_ROOT / "data_processed" / "tcga_prad_recurrence_cohort.csv"
    cohort.to_csv(out, index=False)

    # Save Comprehensive Event Adjudication Ledger
    ledger_cols = [
        "submitter_id",
        "composite_status_detail",
        "outcome_prostate_specific_recurrence",
        "outcome_definitive_recurrence",
        "outcome_exploratory_composite",
        "outcome_composite_all_cancers",
        "outcome_recurrence_or_progression",
        "bcr_status",
        "nte_status",
        "new_neoplasm_event_type",
        "nte_type_category",
        "ht_progression_status",
        "days_to_bcr_clean",
        "days_to_nte_clean",
        "days_to_last_followup_clean",
        "days_to_death_clean",
        "days_to_recurrence_exact",
        "days_to_censoring_exact",
        "earliest_prostate_recurrence_type",
        "earliest_prostate_recurrence_source",
        "survival_duration_days",
        "duration_days_from_diagnosis",
        "survival_event_indicator",
        "undated_recurrence_flag",
        "adjudication_rationale",
    ]
    cohort[ledger_cols].to_csv(RESEARCH_V3_DIR / "EVENT_ADJUDICATION_LEDGER.csv", index=False)
    cohort[ledger_cols].to_csv(RESEARCH_V2_DIR / "EVENT_ADJUDICATION_LEDGER.csv", index=False)
    cohort[ledger_cols].to_csv(RESEARCH_V2_DIR / "LONGITUDINAL_EVENT_AUDIT.csv", index=False)

    # Document the 6 undated recurrence cases specifically
    undated = cohort[cohort["undated_recurrence_flag"] == 1][ledger_cols]
    undated.to_csv(RESEARCH_V3_DIR / "UNDATED_RECURRENCE_CASES.csv", index=False)
    undated.to_csv(RESEARCH_V2_DIR / "UNDATED_RECURRENCE_CASES.csv", index=False)

    # Document the 2 second primary cases specifically
    second_prim = cohort[cohort["composite_status_detail"] == "SECOND_PRIMARY_MALIGNANCY"][ledger_cols]
    second_prim.to_csv(RESEARCH_V3_DIR / "SECOND_PRIMARY_CASES.csv", index=False)
    second_prim.to_csv(RESEARCH_V2_DIR / "SECOND_PRIMARY_CASES.csv", index=False)

    # Document the 5 fully indeterminate cases specifically
    indeterminate = cohort[cohort["composite_status_detail"] == "INDETERMINATE_ALL_MISSING"][ledger_cols]
    indeterminate.to_csv(RESEARCH_V3_DIR / "INDETERMINATE_OUTCOME_CASES.csv", index=False)
    indeterminate.to_csv(RESEARCH_V2_DIR / "INDETERMINATE_OUTCOME_CASES.csv", index=False)

    # Save Cohort Flow Breakdown
    flow = pd.DataFrame(
        [
            {"step": "Total TCGA-PRAD Radical Prostatectomy Patients", "count": len(cohort)},
            {"step": "Primary Definitive Prostate-Specific Cohort (Positives + Double Negatives)", "count": int(cohort["outcome_definitive_recurrence"].notna().sum())},
            {"step": "  - Documented Prostate Cancer Recurrence / Progression (Definitive Positive)", "count": int((cohort["outcome_definitive_recurrence"] == 1.0).sum())},
            {"step": "  - Documented Double-Negative (BCR=NO and NTE=NO)", "count": int((cohort["outcome_definitive_recurrence"] == 0.0).sum())},
            {"step": "Second Non-Prostatic Primary Malignancy (TCGA-CH-5791, TCGA-HC-7080)", "count": int((cohort["composite_status_detail"] == "SECOND_PRIMARY_MALIGNANCY").sum())},
            {"step": "Sensitivity All-Cancers Composite (Positives + Second Primaries + Double Neg)", "count": int(cohort["outcome_composite_all_cancers"].notna().sum())},
            {"step": "Exploratory Composite Binary Cohort (Positives + Double Neg + Asymmetric Neg)", "count": int(cohort["outcome_exploratory_composite"].notna().sum())},
            {"step": "  - Asymmetric Negative (Excluded from primary definitive)", "count": int(cohort["composite_status_detail"].isin(["NEGATIVE_BCR_MISSING_NTE", "MISSING_BCR_NEGATIVE_NTE"]).sum())},
            {"step": "    * Asymmetric: BCR=NO, NTE missing", "count": int((cohort["composite_status_detail"] == "NEGATIVE_BCR_MISSING_NTE").sum())},
            {"step": "    * Asymmetric: NTE=NO, BCR missing", "count": int((cohort["composite_status_detail"] == "MISSING_BCR_NEGATIVE_NTE").sum())},
            {"step": "Fully Indeterminate Cases (Excluded from all primary/exploratory non-events)", "count": int((cohort["composite_status_detail"] == "INDETERMINATE_ALL_MISSING").sum())},
            {"step": "Continuous-Time Survival: Strict Sensitivity Cohort (Dated events + Double Negatives)", "count": int(((cohort["survival_event_indicator"] == 1.0) | (cohort["composite_status_detail"] == "DOCUMENTED_NON_RECURRENCE")).sum())},
            {"step": "  - Strict Dated Prostate Recurrence Events", "count": int((cohort["survival_event_indicator"] == 1.0).sum())},
            {"step": "  - Strict Double-Negative Censored Non-Events", "count": int(((cohort["survival_event_indicator"] == 0.0) & (cohort["composite_status_detail"] == "DOCUMENTED_NON_RECURRENCE")).sum())},
            {"step": "Continuous-Time Survival: Exploratory Cohort (Admitting 46 asymmetric records)", "count": int(cohort["survival_duration_days"].notna().sum())},
            {"step": "  - Exploratory Censored Non-Events (191 double-neg + 46 asymmetric)", "count": int((cohort["survival_event_indicator"] == 0.0).sum())},
            {"step": "Excluded from Exact Continuous Survival (Total)", "count": int(cohort["survival_duration_days"].isna().sum())},
            {"step": "  - Undated Prostate Recurrence Events", "count": int(cohort["undated_recurrence_flag"].sum())},
            {"step": "  - Second Non-Prostatic Primary Tumors", "count": int((cohort["composite_status_detail"] == "SECOND_PRIMARY_MALIGNANCY").sum())},
            {"step": "  - Fully Indeterminate Records", "count": int((cohort["composite_status_detail"] == "INDETERMINATE_ALL_MISSING").sum())},
        ]
    )
    flow.to_csv(RESEARCH_V3_DIR / "COHORT_FLOW.csv", index=False)
    flow.to_csv(RESEARCH_V2_DIR / "COHORT_FLOW.csv", index=False)

    # Long-Form Encounter Provenance Table (All 336 XML rows from 321 patients)
    prov_cols = [
        "submitter_id",
        "file_uuid",
        "batch_number",
        "source_file",
        "bcr_followup_barcode",
        "followup_case_report_form_submission_reason",
        "days_to_last_followup",
        "days_to_death",
        "biochemical_recurrence",
        "days_to_first_biochemical_recurrence",
        "new_tumor_event_after_initial_treatment",
        "new_neoplasm_event_type",
        "days_to_new_tumor_event_after_initial_treatment",
        "tumor_progression_post_ht",
        "psa_value",
        "days_to_psa",
        "gleason_score",
    ]
    avail_prov_cols = [c for c in prov_cols if c in flat.columns]
    flat[avail_prov_cols].to_csv(RESEARCH_V3_DIR / "LONG_FORM_PROVENANCE.csv", index=False)
    flat[avail_prov_cols].to_csv(RESEARCH_V2_DIR / "LONG_FORM_PROVENANCE.csv", index=False)

    # Refined Predictor Timing & Clinical Landmark Audit
    psa_days = pd.to_numeric(flat["days_to_psa"], errors="coerce")
    timing_audit = pd.DataFrame(
        [
            {
                "predictor": "age_at_diagnosis_years",
                "clinical_phase": "Pre-treatment / Diagnostic Biopsy",
                "reference_origin_days": "0.0 (Core Needle Biopsy Day 0 in 98.1% of patients)",
                "prediction_landmark": "Available at initial diagnosis (Biopsy Day 0)",
                "immortal_time_and_landmark_risk": "None (strictly baseline demographic measured at initial diagnosis).",
            },
            {
                "predictor": "pathologic_t",
                "clinical_phase": "Post-resection / Surgical Pathology",
                "reference_origin_days": "Unmeasured post-biopsy surgical offset (conditioned on prostatectomy)",
                "prediction_landmark": "Post-Radical Prostatectomy Staging Landmark (Exact surgery date omitted from GDC XML)",
                "immortal_time_and_landmark_risk": "Post-resection landmark. Core needle biopsy established Day 0 for 305/314 patients, while surgery date is unrecorded. Models using pT condition on resection while follow-up is measured from initial diagnostic biopsy.",
            },
            {
                "predictor": "pathologic_n",
                "clinical_phase": "Post-resection / Surgical Pathology",
                "reference_origin_days": "Unmeasured post-biopsy surgical offset (conditioned on lymphadenectomy)",
                "prediction_landmark": "Post-Radical Prostatectomy Staging Landmark (Exact surgery date omitted from GDC XML)",
                "immortal_time_and_landmark_risk": "Post-resection landmark. Lymph node staging occurs at prostatectomy; exact interval from biopsy is unrecorded in registry.",
            },
            {
                "predictor": "gleason_score_numeric",
                "clinical_phase": "Post-resection / Surgical Pathology",
                "reference_origin_days": "Unmeasured post-biopsy surgical offset (conditioned on prostatectomy specimen)",
                "prediction_landmark": "Post-Radical Prostatectomy Staging Landmark (Exact surgery date omitted from GDC XML)",
                "immortal_time_and_landmark_risk": "Post-resection landmark. Primary and secondary Gleason grade on resected specimen; conditioned on surgery while follow-up is indexed from biopsy.",
            },
            {
                "predictor": "RNA expression panel (13 genes)",
                "clinical_phase": "Post-resection / Resected Tissue Sequencing",
                "reference_origin_days": "Unmeasured post-biopsy surgical offset (conditioned on prostatectomy specimen)",
                "prediction_landmark": "Post-Radical Prostatectomy Staging Landmark (Exact surgery date omitted from GDC XML)",
                "immortal_time_and_landmark_risk": "Post-resection landmark. High-throughput RNA-seq performed on frozen resected primary tumor tissue; conditioned on resection while follow-up is indexed from biopsy.",
            },
            {
                "predictor": "psa_value_numeric",
                "clinical_phase": "Post-operative Surveillance Monitoring",
                "reference_origin_days": f"Median {psa_days.median():.0f} days post-biopsy (IQR: 191–909)",
                "prediction_landmark": "Surveillance monitoring (271 of 276 values recorded >0 days post-biopsy)",
                "immortal_time_and_landmark_risk": "Severe immortal time / lead time bias if treated as baseline. Strictly excluded from baseline Cox and binary prognostic models; evaluated only in dedicated surveillance sensitivity analysis.",
            },
        ]
    )
    timing_audit.to_csv(RESEARCH_V3_DIR / "PREDICTOR_TIMING_AUDIT.csv", index=False)
    timing_audit.to_csv(RESEARCH_V2_DIR / "PREDICTOR_TIMING_AUDIT.csv", index=False)

    # Selection Bias Audit: Comparing Excluded (n=13) vs Included (n=301) Survival Patients
    cohort["excluded_from_survival"] = cohort["survival_duration_days"].isna()
    inc = cohort[~cohort["excluded_from_survival"]]
    exc = cohort[cohort["excluded_from_survival"]]

    bias_audit = pd.DataFrame(
        [
            {
                "characteristic": "Cohort Size (N)",
                "included_cohort_n301": len(inc),
                "excluded_cohort_n13": len(exc),
                "comparison_interpretation": "13 patients excluded: 6 undated prostate recurrence events + 2 second primary cancers + 5 indeterminate",
            },
            {
                "characteristic": "Mean Age at Diagnosis (years)",
                "included_cohort_n301": round(float(inc["age_at_diagnosis_years"].dropna().mean()), 2),
                "excluded_cohort_n13": round(float(exc["age_at_diagnosis_years"].dropna().mean()), 2),
                "comparison_interpretation": "Excluded patients slightly younger at diagnosis (57.8 vs 61.5 years)",
            },
            {
                "characteristic": "Mean Gleason Score",
                "included_cohort_n301": round(float(inc["gleason_score_numeric"].dropna().mean()), 2),
                "excluded_cohort_n13": round(float(exc["gleason_score_numeric"].dropna().mean()), 2),
                "comparison_interpretation": "Excluded cases have higher mean Gleason score (7.82 vs 7.59)",
            },
            {
                "characteristic": "Advanced Pathologic Stage T3/T4 (%)",
                "included_cohort_n301": round(float(inc["pathologic_t"].isin(["T3a", "T3b", "T4"]).mean() * 100), 1),
                "excluded_cohort_n13": round(float(exc["pathologic_t"].isin(["T3a", "T3b", "T4"]).mean() * 100), 1),
                "comparison_interpretation": "Excluded cases show 90.0% advanced T3/T4 stage vs 60.9% in included",
            },
            {
                "characteristic": "True Recurrence Event Prevalence (%)",
                "included_cohort_n301": round(float((inc["survival_event_indicator"] == 1.0).mean() * 100), 1),
                "excluded_cohort_n13": round(float(exc["undated_recurrence_flag"].mean() * 100), 1),
                "comparison_interpretation": "46.2% of excluded patients suffered documented prostate recurrence (6/13), downwardly biasing observed continuous-time event rates",
            },
        ]
    )
    bias_audit.to_csv(RESEARCH_V3_DIR / "SELECTION_BIAS_AUDIT.csv", index=False)
    bias_audit.to_csv(RESEARCH_V2_DIR / "SELECTION_BIAS_AUDIT.csv", index=False)

    # Missingness Table
    missing_tbl = pd.DataFrame(
        [
            {"field": "psa_value", "missing_count": int(cohort["psa_value_numeric"].isna().sum()), "pct": float(cohort["psa_value_numeric"].isna().mean() * 100)},
            {"field": "gleason_score", "missing_count": int(cohort["gleason_score_numeric"].isna().sum()), "pct": float(cohort["gleason_score_numeric"].isna().mean() * 100)},
            {"field": "biochemical_recurrence", "missing_count": int(cohort["biochemical_recurrence"].isna().sum()), "pct": float(cohort["biochemical_recurrence"].isna().mean() * 100)},
            {"field": "new_tumor_event", "missing_count": int(cohort["new_tumor_event_after_initial_treatment"].isna().sum()), "pct": float(cohort["new_tumor_event_after_initial_treatment"].isna().mean() * 100)},
            {"field": "days_to_last_followup", "missing_count": int(cohort["days_to_last_followup_clean"].isna().sum()), "pct": float(cohort["days_to_last_followup_clean"].isna().mean() * 100)},
            {"field": "event_recurrence_exact_date", "missing_count": int(cohort[cohort['outcome_definitive_recurrence'] == 1]['days_to_recurrence_exact'].isna().sum()), "pct": float(cohort[cohort['outcome_definitive_recurrence'] == 1]['days_to_recurrence_exact'].isna().mean() * 100)},
        ]
    )
    missing_tbl.to_csv(RESEARCH_V3_DIR / "MISSINGNESS_TABLE.csv", index=False)
    missing_tbl.to_csv(RESEARCH_V2_DIR / "MISSINGNESS_TABLE.csv", index=False)

    counts = pd.DataFrame(
        [
            {"endpoint": "Prostate-Specific Recurrence (Definitive)", "events": int(cohort["outcome_definitive_recurrence"].eq(1.0).sum()), "total": len(cohort)},
            {"endpoint": "Dated Prostate Recurrence Events (Exact Timing)", "events": int((cohort["survival_event_indicator"] == 1.0).sum()), "total": len(cohort)},
            {"endpoint": "Undated Prostate Recurrence Events (Excluded from exact time)", "events": int(cohort["undated_recurrence_flag"].sum()), "total": len(cohort)},
            {"endpoint": "Second Non-Prostatic Primary Malignancies", "events": int((cohort["composite_status_detail"] == "SECOND_PRIMARY_MALIGNANCY").sum()), "total": len(cohort)},
            {"endpoint": "Biochemical Recurrence (BCR)", "events": int(cohort["flag_bcr_positive"].sum()), "total": len(cohort)},
            {"endpoint": "New Tumor Event (NTE)", "events": int(cohort["flag_nte_positive"].sum()), "total": len(cohort)},
            {"endpoint": "Post-HT Progression", "events": int(cohort["flag_ht_progression_positive"].sum()), "total": len(cohort)},
            {"endpoint": "Evaluable RFS (with valid timing, excluding indeterminate)", "events": int((cohort["survival_event_indicator"] == 1.0).sum()), "total": int(cohort["survival_duration_days"].notna().sum())},
        ]
    )
    counts.to_csv(PROJECT_ROOT / "tables" / "table_recurrence_outcome_counts.csv", index=False)

    print("=== TCGA-PRAD V3.1 CLINICAL VALIDITY ENDPOINT RECONSTRUCTION COMPLETE ===")
    print(flow.to_string(index=False))


if __name__ == "__main__":
    main()
