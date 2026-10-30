from __future__ import annotations

import re
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]


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


def parse_boolean_flag(val: object) -> int:
    """Safely parse boolean clinical response from XML, handling multi-encounter strings."""
    if pd.isna(val):
        return 0
    text = str(val).upper()
    tokens = [t.strip() for t in text.split("|")]
    # If any encounter documented a positive recurrence/progression event, mark as 1
    for t in tokens:
        if t in {"YES", "TRUE", "1", "PROGRESSION", "RECURRENCE", "BIOCHEMICAL", "DISTANT"}:
            return 1
    return 0


def main() -> None:
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
    ]
    clinical = flat[[c for c in keep if c in flat.columns]].drop_duplicates("submitter_id", keep="first").copy()
    cohort = base.merge(clinical, on="submitter_id", how="left")

    # Clean longitudinal follow-up and event days
    cohort["days_to_last_followup_clean"] = cohort["days_to_last_followup_y"].apply(lambda v: parse_numeric_tokens(v, max))
    # Fallback to initial followup if missing
    cohort["days_to_last_followup_clean"] = cohort["days_to_last_followup_clean"].fillna(
        cohort["days_to_last_followup_x"].apply(lambda v: parse_numeric_tokens(v, max))
    )

    cohort["days_to_bcr_clean"] = cohort["days_to_first_biochemical_recurrence"].apply(lambda v: parse_numeric_tokens(v, min))
    cohort["days_to_nte_clean"] = cohort["days_to_new_tumor_event_after_initial_treatment"].apply(lambda v: parse_numeric_tokens(v, min))

    # Parse distinct endpoint components
    cohort["flag_bcr"] = cohort["biochemical_recurrence"].apply(parse_boolean_flag)
    cohort["flag_nte"] = cohort["new_tumor_event_after_initial_treatment"].apply(parse_boolean_flag)
    cohort["flag_progression_ht"] = cohort["tumor_progression_post_ht"].apply(parse_boolean_flag)

    # Composite recurrence endpoint (BCR or NTE or post-HT progression)
    cohort["outcome_recurrence_or_progression"] = cohort[["flag_bcr", "flag_nte", "flag_progression_ht"]].max(axis=1)

    # Clean numeric clinical predictors
    cohort["psa_value_numeric"] = pd.to_numeric(cohort.get("psa_value"), errors="coerce")
    gleason_source = cohort.get("gleason_score_y", cohort.get("gleason_score_x", cohort.get("gleason_score")))
    cohort["gleason_score_numeric"] = pd.to_numeric(gleason_source, errors="coerce")

    # Longitudinal survival time calculation
    time_to_event = []
    for _, row in cohort.iterrows():
        is_event = row["outcome_recurrence_or_progression"] == 1
        t_candidates = []
        if is_event:
            if pd.notna(row["days_to_bcr_clean"]) and row["days_to_bcr_clean"] > 0:
                t_candidates.append(row["days_to_bcr_clean"])
            if pd.notna(row["days_to_nte_clean"]) and row["days_to_nte_clean"] > 0:
                t_candidates.append(row["days_to_nte_clean"])
            if not t_candidates and pd.notna(row["days_to_last_followup_clean"]) and row["days_to_last_followup_clean"] > 0:
                t_candidates.append(row["days_to_last_followup_clean"])
        else:
            if pd.notna(row["days_to_last_followup_clean"]) and row["days_to_last_followup_clean"] > 0:
                t_candidates.append(row["days_to_last_followup_clean"])
        time_to_event.append(min(t_candidates) if t_candidates else np.nan)

    cohort["time_to_recurrence_days"] = time_to_event

    out = PROJECT_ROOT / "data_processed" / "tcga_prad_recurrence_cohort.csv"
    cohort.to_csv(out, index=False)

    # Summary table
    counts = pd.DataFrame(
        [
            {"endpoint": "Composite Recurrence/Progression", "events": int(cohort["outcome_recurrence_or_progression"].sum()), "total": len(cohort)},
            {"endpoint": "Biochemical Recurrence (BCR)", "events": int(cohort["flag_bcr"].sum()), "total": len(cohort)},
            {"endpoint": "New Tumor Event (NTE)", "events": int(cohort["flag_nte"].sum()), "total": len(cohort)},
            {"endpoint": "Post-HT Progression", "events": int(cohort["flag_progression_ht"].sum()), "total": len(cohort)},
            {"endpoint": "Evaluable RFS (with timing)", "events": int(cohort.dropna(subset=['time_to_recurrence_days'])['outcome_recurrence_or_progression'].sum()), "total": int(cohort['time_to_recurrence_days'].notna().sum())},
        ]
    )
    counts.to_csv(PROJECT_ROOT / "tables" / "table_recurrence_outcome_counts.csv", index=False)
    print(f"Wrote PRAD recurrence cohort: {len(cohort):,} patients, events={int(cohort['outcome_recurrence_or_progression'].sum())}")
    print(counts.to_string(index=False))


if __name__ == "__main__":
    main()
