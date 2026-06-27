from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = PROJECT_ROOT
sys.path.append(str(REPO_ROOT / "_shared" / "02_common_functions"))

from rigorous_validation import validate_feature_sets  # noqa: E402


def main() -> None:
    recurrence = pd.read_csv(PROJECT_ROOT / "data_processed" / "tcga_prad_recurrence_cohort.csv")
    expr = pd.read_csv(PROJECT_ROOT / "data_processed" / "tcga_prad_clinical_expression_cohort.csv")
    data = recurrence.merge(
        expr[["submitter_id"] + [c for c in expr.columns if c.startswith("expr_")]],
        on="submitter_id",
        how="left",
    )
    clinical = ["age_at_diagnosis_years", "pathologic_t", "pathologic_n", "gleason_score_numeric", "psa_value_numeric"]
    expression = [column for column in data.columns if column.startswith("expr_")]
    validate_feature_sets(
        data,
        "outcome_recurrence_or_progression",
        {
            "clinical_only": clinical,
            "clinical_expression": clinical + expression,
        },
        PROJECT_ROOT / "tables",
        PROJECT_ROOT / "figures",
        "prad",
    )
    print("Wrote PRAD rigorous validation tables and figures")


if __name__ == "__main__":
    main()
