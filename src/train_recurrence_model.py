from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = PROJECT_ROOT
sys.path.append(str(REPO_ROOT / "_shared" / "02_common_functions"))

from tabular_risk_model import train_logistic_score  # noqa: E402


def main() -> None:
    enriched = PROJECT_ROOT / "data_processed" / "tcga_prad_clinical_expression_cohort.csv"
    recurrence = PROJECT_ROOT / "data_processed" / "tcga_prad_recurrence_cohort.csv"
    if enriched.exists():
        data = pd.read_csv(recurrence).merge(
            pd.read_csv(enriched)[["submitter_id"] + [c for c in pd.read_csv(enriched, nrows=1).columns if c.startswith("expr_")]],
            on="submitter_id",
            how="left",
        )
    else:
        data = pd.read_csv(recurrence)
    predictors = ["age_at_diagnosis_years", "pathologic_t", "pathologic_n", "gleason_score_numeric", "psa_value_numeric"]
    predictors += [c for c in data.columns if c.startswith("expr_")]
    predictors = [c for c in predictors if c in data.columns]
    train_logistic_score(
        data,
        predictors,
        "outcome_recurrence_or_progression",
        PROJECT_ROOT / "models",
        PROJECT_ROOT / "tables",
        min_events=10,
    )
    print("Trained PRAD recurrence/progression model when event counts allowed")


if __name__ == "__main__":
    main()

