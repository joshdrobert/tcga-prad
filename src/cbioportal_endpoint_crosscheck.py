from __future__ import annotations

import json
import urllib.request
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def fetch_cbioportal_prad_tcga() -> pd.DataFrame:
    url = "https://www.cbioportal.org/api/studies/prad_tcga/clinical-data?clinicalDataType=PATIENT"
    request = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "OpenSpecialtyRiskAtlas/0.1"})
    with urllib.request.urlopen(request, timeout=120) as response:
        data = json.loads(response.read().decode("utf-8"))
    long = pd.DataFrame(data)
    wide = long.pivot_table(index="patientId", columns="clinicalAttributeId", values="value", aggfunc="first").reset_index()
    return wide


def main() -> None:
    cbio = fetch_cbioportal_prad_tcga()
    cbio.to_csv(PROJECT_ROOT / "data_raw" / "cbioportal_prad_tcga_clinical.csv", index=False)
    local = pd.read_csv(PROJECT_ROOT / "data_processed" / "tcga_prad_recurrence_cohort.csv")
    local["patientId"] = local["submitter_id"]
    merged = local.merge(cbio, on="patientId", how="left")
    cbio_bcr = merged.get("BIOCHEMICAL_RECURRENCE_INDICATOR")
    merged["cbio_bcr_positive"] = cbio_bcr.fillna("").astype(str).str.upper().str.contains("YES|Y|TRUE|1").astype(int)
    comparable = merged.loc[cbio_bcr.notna()].copy()
    if len(comparable):
        agreement = (comparable["outcome_recurrence_or_progression"].astype(int) == comparable["cbio_bcr_positive"].astype(int)).mean()
    else:
        agreement = float("nan")
    summary = pd.DataFrame(
        [
            {
                "local_patients": len(local),
                "cbioportal_patients": len(cbio),
                "matched_patients": int(merged["BIOCHEMICAL_RECURRENCE_INDICATOR"].notna().sum()) if "BIOCHEMICAL_RECURRENCE_INDICATOR" in merged else 0,
                "agreement_fraction": agreement,
                "local_events_in_comparable": int(comparable["outcome_recurrence_or_progression"].sum()) if len(comparable) else 0,
                "cbio_bcr_events": int(comparable["cbio_bcr_positive"].sum()) if len(comparable) else 0,
            }
        ]
    )
    merged.to_csv(PROJECT_ROOT / "data_processed" / "tcga_prad_cbioportal_endpoint_crosscheck.csv", index=False)
    summary.to_csv(PROJECT_ROOT / "tables" / "table_cbioportal_endpoint_crosscheck.csv", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()

