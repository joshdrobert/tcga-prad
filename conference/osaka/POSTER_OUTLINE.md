# POSTER & ORAL PRESENTATION OUTLINES (FOR FUTURE REPURPOSED VENUE)

## Poster Structure (90cm x 120cm Vertical Format)
1. **Header Panel**: Title, Author Affiliations, GitHub / Data Provenance QR Code.
2. **Column 1: Clinical Problem & Cohort Architecture**:
   - The challenge of endpoint definition in open oncology databases.
   - TCGA-PRAD inclusion criteria (n=314) and molecular subset (n=228).
   - Cohort attrition diagram illustrating missing RNA sequencing.
3. **Column 2: Endpoint Concordance & Discrepancy Breakdown**:
   - Local XML composite recurrence vs cBioPortal Biochemical Recurrence (BCR).
   - Confusion matrix demonstrating 31 discordances across 269 matched patients.
   - Classification of discordances into clinical new tumor events vs distant metastases.
4. **Column 3: Model Discrimination, Calibration & Failure Modes**:
   - Out-of-fold AUROC comparisons (Clinical vs Molecular vs Paired Complete-Case).
   - Calibration curves showing severe miscalibration (slope 0.547) when missing features are median-imputed.
   - Decision curves illustrating net benefit boundaries.
5. **Footer: Translation Lessons & Open Science Checklist**:
   - Recommendations for clinical machine learning reporting (PROBAST+AI).

## 10-Minute Oral Presentation Outline
- **0:00 - 2:00: The Illusion of Accuracy in Clinical ML**: Why published recurrence models fail when applied to real cohorts.
- **2:00 - 4:00: Endpoint Ground Truth**: Auditing the discrepancy between GDC clinical supplements and cBioPortal curated indicators.
- **4:00 - 6:00: The Peril of Silent Imputation**: How median imputation of 86 unsequenced patients masked true molecular signal.
- **6:00 - 8:00: Paired Complete-Case Re-evaluation**: Transparent findings showing clinical staging remains primary.
- **8:00 - 10:00: Discussion & Open Science Practice**: Mandatory auditing guidelines for open biomedical benchmarks.
