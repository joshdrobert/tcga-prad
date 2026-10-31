from __future__ import annotations

import json
import urllib.request
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESEARCH_V2_DIR = PROJECT_ROOT / "research" / "v2"


def fetch_cbioportal_prad_tcga() -> pd.DataFrame:
    raw_path = PROJECT_ROOT / "data_raw" / "cbioportal_prad_tcga_clinical.csv"
    if raw_path.exists():
        return pd.read_csv(raw_path)

    url = "https://www.cbioportal.org/api/studies/prad_tcga/clinical-data?clinicalDataType=PATIENT"
    try:
        request = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "OpenSpecialtyRiskAtlas/0.1"})
        with urllib.request.urlopen(request, timeout=30) as response:
            data = json.loads(response.read().decode("utf-8"))
        long = pd.DataFrame(data)
        wide = long.pivot_table(index="patientId", columns="clinicalAttributeId", values="value", aggfunc="first").reset_index()
        wide.to_csv(raw_path, index=False)
        return wide
    except Exception as e:
        print(f"Warning: Could not fetch from cBioPortal API ({e}), using local raw file if available.")
        if raw_path.exists():
            return pd.read_csv(raw_path)
        raise


def main() -> None:
    RESEARCH_V2_DIR.mkdir(parents=True, exist_ok=True)
    cbio = fetch_cbioportal_prad_tcga()
    local = pd.read_csv(PROJECT_ROOT / "data_processed" / "tcga_prad_recurrence_cohort.csv")
    local["patientId"] = local["submitter_id"]

    merged = local.merge(cbio, on="patientId", how="left")
    cbio_bcr = merged.get("BIOCHEMICAL_RECURRENCE_INDICATOR")

    # Match indicator: YES vs NO (excluding unassessed/missing in cBio)
    merged["cbio_bcr_positive"] = cbio_bcr.fillna("").astype(str).str.upper().str.contains("YES|Y|TRUE|1").astype(int)

    # Evaluable matched: cases present in cBioPortal with non-null BCR indicator
    evaluable_mask = cbio_bcr.notna()
    comparable = merged.loc[evaluable_mask].copy()

    local_y = comparable["outcome_recurrence_or_progression"].astype(int)
    cbio_y = comparable["cbio_bcr_positive"].astype(int)

    concordant_non_events = int(((local_y == 0) & (cbio_y == 0)).sum())
    concordant_events = int(((local_y == 1) & (cbio_y == 1)).sum())
    local_only_events = int(((local_y == 1) & (cbio_y == 0)).sum())
    cbio_only_events = int(((local_y == 0) & (cbio_y == 1)).sum())
    discordant_cases = local_only_events + cbio_only_events

    agreement_fraction = float((local_y == cbio_y).mean()) if len(comparable) else float("nan")

    summary = pd.DataFrame(
        [
            {
                "local_patients": len(local),
                "cbioportal_matched_records": int(merged["patientId"].isin(cbio["patientId"]).sum()),
                "evaluable_matched": len(comparable),
                "agreement_fraction": agreement_fraction,
                "concordant_non_events": concordant_non_events,
                "concordant_events": concordant_events,
                "discordant_cases": discordant_cases,
                "local_only_events": local_only_events,
                "cbio_only_events": cbio_only_events,
                "local_events_in_comparable": int(local_y.sum()),
                "cbio_bcr_events": int(cbio_y.sum()),
            }
        ]
    )

    merged.to_csv(PROJECT_ROOT / "data_processed" / "tcga_prad_cbioportal_endpoint_crosscheck.csv", index=False)
    summary.to_csv(PROJECT_ROOT / "tables" / "table_cbioportal_endpoint_crosscheck.csv", index=False)
    summary.to_csv(RESEARCH_V2_DIR / "ENDPOINT_CONCORDANCE.csv", index=False)

    # Detailed discordant cases table
    discordant_df = comparable[local_y != cbio_y][
        [
            "submitter_id",
            "outcome_recurrence_or_progression",
            "bcr_status",
            "nte_status",
            "ht_progression_status",
            "BIOCHEMICAL_RECURRENCE_INDICATOR",
            "days_to_last_followup_clean",
            "days_to_recurrence_exact",
        ]
    ].rename(columns={"outcome_recurrence_or_progression": "local_composite_recurrence", "BIOCHEMICAL_RECURRENCE_INDICATOR": "cbio_bcr_indicator"})

    discordant_df.to_csv(PROJECT_ROOT / "tables" / "table_cbioportal_discordant_patients.csv", index=False)
    discordant_df.to_csv(RESEARCH_V2_DIR / "DISCORDANT_CASES.csv", index=False)

    print(f"=== TCGA-PRAD V2 cBioPortal Cross-Check Complete ===")
    print(f"Evaluable Matched Cases: {len(comparable)}")
    print(f"Agreement Fraction: {agreement_fraction:.4f}")
    print(f"Concordant Events: {concordant_events}, Concordant Non-Events: {concordant_non_events}")
    print(f"Discordant Cases: {discordant_cases} (Local Only: {local_only_events}, cBio Only: {cbio_only_events})")


if __name__ == "__main__":
    main()
