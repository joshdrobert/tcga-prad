from __future__ import annotations

import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = PROJECT_ROOT
sys.path.append(str(REPO_ROOT / "_shared" / "02_common_functions"))

from gdc_molecular import CANCER_PANELS, build_expression_matrix, download_expression_files, merge_expression  # noqa: E402


def main() -> None:
    manifest = PROJECT_ROOT / "data_raw" / "TCGA-PRAD_gdc_file_manifest.csv"
    expr_dir = PROJECT_ROOT / "data_raw" / "rna_seq_star_counts"
    status_csv = PROJECT_ROOT / "tables" / "table_rna_seq_download_status.csv"
    matrix_csv = PROJECT_ROOT / "data_processed" / "tcga_prad_expression_panel.csv"
    enriched_csv = PROJECT_ROOT / "data_processed" / "tcga_prad_clinical_expression_cohort.csv"
    status = download_expression_files(manifest, expr_dir)
    status.to_csv(status_csv, index=False)
    build_expression_matrix(status_csv, CANCER_PANELS["PRAD"], matrix_csv)
    merged = merge_expression(PROJECT_ROOT / "data_processed" / "tcga_prad_analysis_cohort.csv", matrix_csv, enriched_csv)
    print(f"Wrote TCGA-PRAD expression-enriched cohort: {len(merged):,} patients")


if __name__ == "__main__":
    main()

