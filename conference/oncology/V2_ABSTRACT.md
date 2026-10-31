# ONCOLOGY INFORMATICS RESEARCH ABSTRACT
## FUTURE SUBMISSION (AUA / AMIA 2027)

**Title:** Endpoint Discordance, Predictor Chronology, and Missing Molecular Data in Prostate Cancer Recurrence Benchmarks  
**Author:** Joshua Robert [1]  
**Affiliation:** [1] Texas A&M School of Engineering Medicine, Houston, TX, USA  
**Target Venues:** American Urological Association (AUA 2027) / American Medical Informatics Association (AMIA 2027)  

---

### Abstract Body

#### Background
Public cancer registries combine heterogeneous clinical encounter forms and incomplete molecular profiles, complicating prognostic model evaluation. We audited longitudinal recurrence ascertainment, clinical predictor timing relative to initial diagnosis, and the incremental predictive value of transcriptomic profiling in post-radical prostatectomy patients.

#### Methods
We analyzed 314 patients in the TCGA prostate adenocarcinoma (TCGA-PRAD) cohort across 336 encounter records. We reconstructed a primary definitive prostate cancer recurrence endpoint restricted to documented positive recurrence (n=70) and documented double-negative recurrence-free cases (n=191), separating two non-prostatic second primary malignancies and excluding 46 asymmetric surveillance records and five fully indeterminate cases. We audited clinical predictor chronology, establishing initial pathologic diagnosis (biopsy Day 0) as the clinical time origin and separating baseline staging (age, stage, Gleason sum) from post-operative surveillance PSA. On identical participants with complete transcriptomic profiles (n=187), we evaluated 5-fold cross-validated logistic regression with paired bootstrap differences in AUROC. Continuous-time recurrence-free survival from diagnosis was evaluated on strict sensitivity (n=255) and exploratory (n=301) cohorts using baseline Cox models without post-baseline surveillance PSA leakage.

#### Results
Comparing specific endpoints, biochemical recurrence agreed 100.0% (269/269) with cBioPortal, whereas disease-free survival agreed in 92.3% (239/259), uncovering 20 patients with XML-documented recurrence labeled disease-free in cBioPortal. Predictor timing analysis revealed XML surveillance PSA was measured at a median of 511 days post-diagnosis (271 of 276 values post-baseline). In the primary paired molecular cohort (n=187, 44 events), pure baseline clinical staging achieved an AUROC of 0.710 (95% CI: 0.627–0.787), whereas adding the 13-gene expression panel yielded an AUROC of 0.703 (95% CI: 0.618–0.785; paired ΔAUROC: -0.008, 95% CI: [-0.073, +0.056], p=0.809). In strict continuous-time survival analysis from initial diagnosis (n=255, 64 dated events, 191 double-negative censored), out-of-fold Cox concordance was 0.649 (apparent 0.686; 1-year RFS: 87.2%, 5-year RFS: 62.7%), compared with 0.646 out-of-fold when admitting 46 asymmetric records (n=301). Schoenfeld residual tests did not detect a statistically significant violation of proportional hazards (p > 0.40 for all baseline covariates).

#### Conclusions
Failure to adjudicate endpoint definitions, separate second primary malignancies, or exclude post-baseline surveillance biomarkers distorts prognostic evaluation. The evaluated 13-gene transcriptomic panel provided no statistically significant discrimination gain over baseline clinical staging alone in this matched benchmark.

---

### Osaka Theranostics 2027 Formal Decision: NO-GO
- **Official Program:** Focused exclusively on theranostics, radioligands ($^{68}\text{Ga}$-PSMA-11 / $^{177}\text{Lu}$-PSMA-617), biodistribution dosimetry, and radiomolecular therapeutics.
- **TCGA-PRAD:** 0% theranostic data. Maintaining scientific integrity requires not submitting genomic recurrence prediction to Osaka.
