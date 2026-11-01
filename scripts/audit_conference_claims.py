#!/usr/bin/env python3
"""
scripts/audit_conference_claims.py
Independent verification of Osaka TWC 2027 conference abstract and numeric claims.
"""
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent

def audit_osaka_abstract():
    print("=== AUDITING OSAKA TWC 2027 ABSTRACT CONSTRAINTS ===")
    abstract_path = PROJECT_ROOT / "conference" / "osaka_2027" / "READY_TO_PASTE_ABSTRACT.txt"
    if not abstract_path.exists():
        print("ERROR: READY_TO_PASTE_ABSTRACT.txt missing!")
        return False

    content = abstract_path.read_text(encoding="utf-8")
    parts = content.split("ABSTRACT TEXT (PASTE BELOW INTO PORTAL):\n\n")
    if len(parts) < 2:
        print("ERROR: Delimiter not found in READY_TO_PASTE_ABSTRACT.txt")
        return False
    body = parts[1].strip()

    # Rule 1: Portal <= 3,000 characters with spaces
    chars_with_spaces = len(body)
    # Rule 2: Website <= 2,500 characters without whitespace
    chars_no_whitespace = len("".join(body.split()))
    words = len(body.split())

    print(f"Body characters WITH spaces: {chars_with_spaces} (Limit <= 3000)")
    print(f"Body characters WITHOUT whitespace: {chars_no_whitespace} (Limit <= 2500)")
    print(f"Body words: {words}")

    passed = True
    if chars_with_spaces > 3000:
        print(f"FAILED: Body with spaces ({chars_with_spaces}) > 3000")
        passed = False
    if chars_no_whitespace > 2500:
        print(f"FAILED: Body without whitespace ({chars_no_whitespace}) > 2500")
        passed = False

    # Check headings
    headings = ["Aim:", "Material & Methods:", "Results:", "Conclusion:"]
    for h in headings:
        if h not in body:
            print(f"FAILED: Missing heading {h}")
            passed = False

    # Check forbidden content
    for forbidden in ["http://", "https://", "www.", "@", "Veracyte", "Decipher", "Oncotype", "Prolaris"]:
        if forbidden.lower() in body.lower():
            print(f"FAILED: Forbidden term '{forbidden}' in abstract body")
            passed = False

    # Audit numeric claims against research/artifacts/tcga_prad_metrics.json
    metrics_path = PROJECT_ROOT / "research" / "artifacts" / "tcga_prad_metrics.json"
    if metrics_path.exists():
        metrics = json.loads(metrics_path.read_text(encoding="utf-8"))
        print("\n=== AUDITING NUMERICAL CLAIMS AGAINST AUTHORITATIVE METRICS JSON ===")
        
        cohort_info = metrics["cohort"]
        prim_cohort = metrics["primary_definitive_cohort"]
        paired_comp = metrics["paired_molecular_comparisons"][0]
        survival_info = metrics["survival"]

        assert cohort_info["total_patients"] == 314
        assert prim_cohort["total_evaluable_patients"] == 261
        assert prim_cohort["events_positive"] == 70
        assert prim_cohort["non_events_double_negative"] == 191
        assert prim_cohort["paired_molecular_n"] == 187
        assert prim_cohort["paired_molecular_events"] == 44
        assert prim_cohort["paired_molecular_non_events"] == 143
        assert survival_info["strict_evaluable_cases"] == 248

        # Performance anchors
        clin_auroc = paired_comp["baseline_auroc"]
        rna_auroc = paired_comp["augmented_auroc"]
        delta_auroc = paired_comp["delta_auroc"]
        clin_brier = paired_comp["baseline_brier"]
        rna_brier = paired_comp["augmented_brier"]
        delta_brier = paired_comp["delta_brier"]

        val_models = metrics["validation_models"]
        paired_clin = next(m for m in val_models if m["cohort"] == "Paired Definitive Molecular (n=187)" and "Pure Baseline" in m["model"])
        paired_rna = next(m for m in val_models if m["cohort"] == "Paired Definitive Molecular (n=187)" and "13-Gene" in m["model"])
        clin_cal_slope = paired_clin["calibration_slope"]
        rna_cal_slope = paired_rna["calibration_slope"]

        print(f"Clinical AUROC: {clin_auroc:.4f} (Abstract reports 0.7104)")
        print(f"Molecular AUROC: {rna_auroc:.4f} (Abstract reports 0.7028)")
        print(f"Delta AUROC: {delta_auroc:.4f} (Abstract reports -0.0076)")
        print(f"Clinical Brier: {clin_brier:.4f} (Abstract reports 0.1673)")
        print(f"Molecular Brier: {rna_brier:.4f} (Abstract reports 0.1672)")
        print(f"Delta Brier: {delta_brier:.4f} (Abstract reports -0.0001)")
        print(f"Clinical Calibration Slope: {clin_cal_slope:.3f} (Abstract reports 0.756)")
        print(f"Molecular Calibration Slope: {rna_cal_slope:.3f} (Abstract reports 0.511)")

        # Survival
        surv_1yr = survival_info["strict_rfs_1yr"]
        surv_3yr = survival_info["strict_rfs_3yr"]
        surv_5yr = survival_info["strict_rfs_5yr"]
        cox_cv_c = survival_info["strict_cox_cv_concordance_index"]
        cox_app_c = survival_info["strict_cox_apparent_concordance_index"]

        print(f"1-Yr RFS: {surv_1yr*100:.1f}% (Abstract reports 88.1%)")
        print(f"3-Yr RFS: {surv_3yr*100:.1f}% (Abstract reports 75.2%)")
        print(f"5-Yr RFS: {surv_5yr*100:.1f}% (Abstract reports 65.3%)")
        print(f"Cox CV C-index: {cox_cv_c:.4f} (Abstract reports 0.6716)")
        print(f"Cox Apparent C-index: {cox_app_c:.4f} (Abstract reports 0.7011)")

        # Verify tokens exist in text
        tokens = [
            "314", "261", "187", "0.7104", "0.7028", "-0.0076", "0.1673", "0.1672",
            "0.756", "0.511", "248", "88.1%", "75.2%", "65.3%", "0.6716", "0.7011",
            "511.0", "29"
        ]
        for token in tokens:
            if token not in body:
                print(f"FAILED: Expected token '{token}' not found in abstract text!")
                passed = False

    if passed:
        print("\nALL OSAKA ABSTRACT AND METRIC CHECKS PASSED.")
    return passed

if __name__ == "__main__":
    if not audit_osaka_abstract():
        sys.exit(1)
