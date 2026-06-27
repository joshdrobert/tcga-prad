from __future__ import annotations

import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = PROJECT_ROOT
sys.path.append(str(REPO_ROOT / "_shared" / "02_common_functions"))

from tcga_clinical import build_tcga_analysis_table, download_clinical_xml, parse_xml_directory, write_summary_tables  # noqa: E402


def main() -> None:
    manifest = PROJECT_ROOT / "data_raw" / "TCGA-PRAD_gdc_file_manifest.csv"
    xml_dir = PROJECT_ROOT / "data_raw" / "clinical_xml"
    processed = PROJECT_ROOT / "data_processed"
    tables = PROJECT_ROOT / "tables"
    processed.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)
    status = download_clinical_xml(manifest, xml_dir)
    status.to_csv(tables / "table_clinical_xml_download_status.csv", index=False)
    parsed = parse_xml_directory(xml_dir)
    parsed.to_csv(processed / "tcga_prad_clinical_xml_flat.csv", index=False)
    analysis = build_tcga_analysis_table(parsed, "TCGA-PRAD")
    analysis.to_csv(processed / "tcga_prad_analysis_cohort.csv", index=False)
    write_summary_tables(analysis, tables)
    print(f"Wrote TCGA-PRAD cohort: {len(analysis):,} patients")


if __name__ == "__main__":
    main()
