from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = PROJECT_ROOT
sys.path.append(str(REPO_ROOT / "_shared" / "02_common_functions"))

from validation_extras import cross_validated_metrics, simple_km_table  # noqa: E402


def main() -> None:
    data = pd.read_csv(PROJECT_ROOT / "data_processed" / "tcga_prad_recurrence_cohort.csv")
    expr = pd.read_csv(PROJECT_ROOT / "data_processed" / "tcga_prad_clinical_expression_cohort.csv")
    data = data.merge(expr[["submitter_id"] + [c for c in expr.columns if c.startswith("expr_")]], on="submitter_id", how="left")
    predictors = ["age_at_diagnosis_years", "pathologic_t", "pathologic_n", "gleason_score_numeric", "psa_value_numeric"] + [c for c in data.columns if c.startswith("expr_")]
    cross_validated_metrics(data, predictors, "outcome_recurrence_or_progression", "prad_recurrence_expression").to_csv(PROJECT_ROOT / "tables" / "table_cross_validation.csv", index=False)
    if "days_to_first_biochemical_recurrence" in data.columns:
        data["recurrence_time_days"] = pd.to_numeric(data["days_to_first_biochemical_recurrence"], errors="coerce").fillna(
            pd.to_numeric(data["days_to_new_tumor_event_after_initial_treatment"], errors="coerce")
        )
        simple_km_table(data, "recurrence_time_days", "outcome_recurrence_or_progression", None, PROJECT_ROOT / "tables" / "table_recurrence_horizon_summary.csv")
    print("Wrote PRAD validation summaries")


if __name__ == "__main__":
    main()

