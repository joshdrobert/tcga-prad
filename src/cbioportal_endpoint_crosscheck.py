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

    # 1. Comparison A: Biochemical Recurrence to Biochemical Recurrence (Strict BCR-to-BCR)
    # Both evaluate biochemical recurrence specifically
    cbio_bcr_raw = merged["BIOCHEMICAL_RECURRENCE_INDICATOR"]
    mask_bcr = cbio_bcr_raw.isin(["YES", "NO"]) & merged["bcr_status"].isin(["DOCUMENTED_POSITIVE", "DOCUMENTED_NEGATIVE"])
    bcr_comp = merged.loc[mask_bcr].copy()
    bcr_local_y = (bcr_comp["bcr_status"] == "DOCUMENTED_POSITIVE").astype(int)
    bcr_cbio_y = (bcr_comp["BIOCHEMICAL_RECURRENCE_INDICATOR"] == "YES").astype(int)

    # 2. Comparison B: Definitive Recurrence vs cBio DFS_STATUS (Disease-Free Survival)
    cbio_dfs_raw = merged["DFS_STATUS"]
    mask_dfs_def = cbio_dfs_raw.isin(["0:DiseaseFree", "1:Recurred/Progressed"]) & merged["outcome_definitive_recurrence"].notna()
    dfs_def_comp = merged.loc[mask_dfs_def].copy()
    dfs_local_def_y = dfs_def_comp["outcome_definitive_recurrence"].astype(int)
    dfs_cbio_def_y = (dfs_def_comp["DFS_STATUS"] == "1:Recurred/Progressed").astype(int)

    # 3. Comparison C: Exploratory Composite Recurrence vs cBio DFS_STATUS
    mask_dfs_exp = cbio_dfs_raw.isin(["0:DiseaseFree", "1:Recurred/Progressed"]) & merged["outcome_exploratory_composite"].notna()
    dfs_exp_comp = merged.loc[mask_dfs_exp].copy()
    dfs_local_exp_y = dfs_exp_comp["outcome_exploratory_composite"].astype(int)
    dfs_cbio_exp_y = (dfs_exp_comp["DFS_STATUS"] == "1:Recurred/Progressed").astype(int)

    # 4. Comparison D: Legacy Composite Recurrence vs cBio BCR (Historical Audit Benchmark)
    mask_legacy = cbio_bcr_raw.notna()
    legacy_comp = merged.loc[mask_legacy].copy()
    legacy_local_y = legacy_comp["outcome_recurrence_or_progression"].astype(int)
    legacy_cbio_y = legacy_comp["BIOCHEMICAL_RECURRENCE_INDICATOR"].fillna("").astype(str).str.upper().str.contains("YES|Y|TRUE|1").astype(int)

    def summarize_comp(name: str, local_y: pd.Series, cbio_y: pd.Series, total_n: int) -> dict:
        concordant_neg = int(((local_y == 0) & (cbio_y == 0)).sum())
        concordant_pos = int(((local_y == 1) & (cbio_y == 1)).sum())
        local_only = int(((local_y == 1) & (cbio_y == 0)).sum())
        cbio_only = int(((local_y == 0) & (cbio_y == 1)).sum())
        discordant = local_only + cbio_only
        agreement = float((local_y == cbio_y).mean()) if len(local_y) else float("nan")
        return {
            "comparison_definition": name,
            "local_patients": total_n,
            "cbioportal_matched_records": int(merged["patientId"].isin(cbio["patientId"]).sum()),
            "evaluable_matched": len(local_y),
            "agreement_fraction": agreement,
            "concordant_non_events": concordant_neg,
            "concordant_events": concordant_pos,
            "discordant_cases": discordant,
            "local_only_events": local_only,
            "cbio_only_events": cbio_only,
            "local_events_in_comparable": int(local_y.sum()),
            "cbio_events_in_comparable": int(cbio_y.sum()),
            "cbio_bcr_events": int(cbio_y.sum()),
        }

    summary_rows = [
        summarize_comp("BCR-to-BCR (Local BCR vs cBio BCR)", bcr_local_y, bcr_cbio_y, len(local)),
        summarize_comp("Primary Definitive Recurrence vs cBio DFS_STATUS", dfs_local_def_y, dfs_cbio_def_y, len(local)),
        summarize_comp("Exploratory Composite Recurrence vs cBio DFS_STATUS", dfs_local_exp_y, dfs_cbio_exp_y, len(local)),
        summarize_comp("Historical Benchmark (Local Composite vs cBio BCR)", legacy_local_y, legacy_cbio_y, len(local)),
    ]
    summary_all = pd.DataFrame(summary_rows)

    # Maintain legacy primary single-row table format for backward compatibility with existing parsers
    # Row 0 is the primary Definitive Recurrence vs cBio DFS_STATUS crosscheck, preserving keys
    summary_primary = pd.DataFrame([summary_rows[1]])
    # Also save the historical row table to retain EXACT backward compatibility for table_cbioportal_endpoint_crosscheck.csv
    # if required by older tests, while providing full multi-endpoint comparisons
    summary_legacy = pd.DataFrame([summary_rows[3]])

    merged.to_csv(PROJECT_ROOT / "data_processed" / "tcga_prad_cbioportal_endpoint_crosscheck.csv", index=False)
    summary_all.to_csv(PROJECT_ROOT / "tables" / "table_cbioportal_endpoint_crosscheck_multidef.csv", index=False)
    summary_all.to_csv(RESEARCH_V2_DIR / "ENDPOINT_CONCORDANCE_MULTIDEF.csv", index=False)
    
    # Keep table_cbioportal_endpoint_crosscheck.csv with historical benchmark + explicit multi-endpoint columns
    summary_legacy.to_csv(PROJECT_ROOT / "tables" / "table_cbioportal_endpoint_crosscheck.csv", index=False)
    summary_all.to_csv(RESEARCH_V2_DIR / "ENDPOINT_CONCORDANCE.csv", index=False)

    # Detailed discordant cases table for Definitive Recurrence vs cBio DFS (20 cases)
    disc_dfs = dfs_def_comp[dfs_local_def_y != dfs_cbio_def_y][
        [
            "submitter_id",
            "outcome_definitive_recurrence",
            "DFS_STATUS",
            "bcr_status",
            "nte_status",
            "ht_progression_status",
            "BIOCHEMICAL_RECURRENCE_INDICATOR",
            "NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT",
            "days_to_last_followup_clean",
            "days_to_recurrence_exact",
        ]
    ].rename(columns={
        "outcome_definitive_recurrence": "local_definitive_recurrence",
        "DFS_STATUS": "cbio_dfs_status",
        "BIOCHEMICAL_RECURRENCE_INDICATOR": "cbio_bcr_indicator",
        "NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT": "cbio_nte_indicator",
    })
    disc_dfs.to_csv(PROJECT_ROOT / "tables" / "table_cbioportal_discordant_patients.csv", index=False)
    disc_dfs.to_csv(RESEARCH_V2_DIR / "DISCORDANT_CASES.csv", index=False)

    print(f"=== TCGA-PRAD V2 cBioPortal Cross-Check Complete ===")
    print(summary_all[["comparison_definition", "evaluable_matched", "agreement_fraction", "discordant_cases"]].to_string(index=False))


if __name__ == "__main__":
    main()
