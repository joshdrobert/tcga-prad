# 9TH THERANOSTICS WORLD CONGRESS (OSAKA 2027)
## ABSTRACT STATUS: WITHHELD / NO-GO (AUDIT CONCLUSION)

**Track:** Artificial Intelligence / Prostate Cancer
**Decision Status:** **NO-GO FOR SUBMISSION (November 1, 2026 Deadline)**
**Character Count:** 1,948 characters (excluding spaces) [Limit: 2,500 characters excluding spaces]

---

### Title
Evaluating Endpoint Harmonization, Molecular Missingness, and Imaging Feasibility in Machine Learning Benchmarks for Prostate Cancer Recurrence

### Aim
To rigorously evaluate the clinical and translational validity of public machine learning benchmarks for prostate cancer recurrence, specifically assessing whether surgical genomic cohorts (TCGA-PRAD) can support theranostic risk stratification, and auditing the computational feasibility of automated lesion-level PSMA-PET/CT quality control.

### Material & Methods
We audited the complete execution graph of open TCGA-PRAD clinical XML supplements (n=314) and paired GDC RNA-seq STAR profiles (n=228). Predictor availability was evaluated across age, pathologic TNM staging, serum PSA, Gleason score, and a 13-gene transcriptomic panel. Model discrimination and calibration were evaluated using stratified 5-fold cross-validation with out-of-fold AUROC, Brier score, and logistic calibration slopes. Discrepancies between local composite recurrence endpoints and cBioPortal biochemical recurrence (BCR) indicators were audited. Feasibility for patient-level PSMA-PET/CT segmentation robustness was assessed using TCIA public data (597 studies, 378 patients).

### Results
In TCGA-PRAD, a clinical-only baseline model achieved an out-of-fold AUROC of 0.689 (95% CI 0.620–0.756) and calibration slope of 0.777. Adding the RNA-seq panel to the full cohort degraded calibration to 0.547 due to silent imputation of 86 patients lacking molecular sequencing (27.4% missingness). In complete-case paired evaluation (n=228), molecular features did not provide statistically significant incremental discrimination over standard clinical staging. Concordance auditing against cBioPortal BCR identified 31 discordances among 269 matched patients, arising from local inclusion of distant new tumor events. Evaluation of volumetric PSMA-PET/CT confirmed that valid 3D lesion failure modeling requires dedicated high-performance GPU infrastructure exceeding single-workstation limits.

### Conclusion
Genomic surgical benchmarks cannot serve as proxies for theranostic patient selection. Transparent reporting of missing molecular data, endpoint discordance, and computational limits is mandatory to preserve scientific integrity before clinical deployment.
