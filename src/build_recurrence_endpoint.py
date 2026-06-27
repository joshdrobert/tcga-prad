from __future__ import annotations

from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def yes_no(series: pd.Series) -> pd.Series:
    text = series.fillna("").astype(str).str.upper()
    return text.str.contains("YES|BIOCHEMICAL|DISTANT|RECURRENCE|PROGRESSION").astype(int)


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
    ]
    clinical = flat[[c for c in keep if c in flat.columns]].drop_duplicates("submitter_id", keep="first")
    cohort = base.merge(clinical, on="submitter_id", how="left")
    recurrence_flags = []
    for col in ["biochemical_recurrence", "new_tumor_event_after_initial_treatment", "tumor_progression_post_ht"]:
        if col in cohort.columns:
            recurrence_flags.append(yes_no(cohort[col]))
    if recurrence_flags:
        cohort["outcome_recurrence_or_progression"] = pd.concat(recurrence_flags, axis=1).max(axis=1)
    else:
        cohort["outcome_recurrence_or_progression"] = pd.NA
    cohort["psa_value_numeric"] = pd.to_numeric(cohort.get("psa_value"), errors="coerce")
    gleason_source = cohort.get("gleason_score_y", cohort.get("gleason_score_x", cohort.get("gleason_score")))
    cohort["gleason_score_numeric"] = pd.to_numeric(gleason_source, errors="coerce")
    out = PROJECT_ROOT / "data_processed" / "tcga_prad_recurrence_cohort.csv"
    cohort.to_csv(out, index=False)
    cohort["outcome_recurrence_or_progression"].value_counts(dropna=False).rename_axis("outcome").reset_index(name="n").to_csv(
        PROJECT_ROOT / "tables" / "table_recurrence_outcome_counts.csv", index=False
    )
    print(f"Wrote PRAD recurrence cohort: {len(cohort):,} patients, events={int(cohort['outcome_recurrence_or_progression'].sum())}")


if __name__ == "__main__":
    main()
