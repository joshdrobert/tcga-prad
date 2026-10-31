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
We analyzed 314 patients in the TCGA prostate adenocarcinoma (TCGA-PRAD) cohort across 336 encounter records. We reconstructed a primary definitive prostate cancer recurrence endpoint restricted to documented positive recurrence (n=70) and documented double-negative recurrence-free cases (n=191), separating two non-prostatic second primary malignancies and excluding 46 asymmetric surveillance records and five fully indeterminate cases. We audited clinical predictor chronology, establishing initial pathologic diagnosis (needle biopsy Day 0 for 98.1% of patients) as the follow-up time origin while recognizing that surgical pathology (Gleason sum, pT, pN) and transcriptomic profiles are conditioned on radical prostatectomy resection with an unmeasured biopsy-to-surgery interval. Postoperative surveillance PSA was separated from baseline staging. On identical participants with complete transcriptomic profiles (n=187, 44 events), we evaluated 5-fold cross-validated logistic regression with paired bootstrap differences in AUROC and Brier score. Continuous-time recurrence-free survival from diagnosis was evaluated on strict sensitivity (n=255) and exploratory (n=301) cohorts using baseline Cox models without post-baseline surveillance PSA leakage.

#### Results
Comparing specific endpoints, biochemical recurrence agreed 100.0% (269/269) with cBioPortal, whereas disease-free survival agreed in 92.3% (239/259), uncovering 20 patients with XML-documented recurrence labeled disease-free in cBioPortal. Predictor timing analysis revealed XML surveillance PSA was measured at a median of 511 days post-diagnosis (271 of 276 values post-baseline). In the primary paired molecular cohort (n=187, 44 events), pure baseline clinical staging achieved an AUROC of 0.710 (95% CI: 0.621–0.791) and Brier score of 0.1673 (95% CI: 0.1366–0.2018), whereas adding the 13-gene expression panel yielded an AUROC of 0.703 (95% CI: 0.605–0.792) and Brier score of 0.1672 (95% CI: 0.1329–0.2057; paired ΔAUROC: -0.008, 95% CI: [-0.073, +0.056]; paired ΔBrier: -0.0001, 95% CI: [-0.0198, +0.0199]). In continuous-time survival analysis from initial diagnosis, median recurrence-free survival was not reached in either cohort. In the strict sensitivity cohort (n=255, 64 dated events, 191 double-negative censored), out-of-fold Cox concordance was 0.649 (apparent 0.686; 1-year RFS: 87.2%, 5-year RFS: 62.7%), compared with 0.646 out-of-fold when admitting 46 asymmetric records (n=301). Schoenfeld residual tests did not detect a statistically significant violation of proportional hazards (p > 0.40 for all baseline covariates).

#### Conclusions
Failure to adjudicate endpoint definitions, separate second primary malignancies, or exclude post-baseline surveillance biomarkers distorts prognostic evaluation. The evaluated 13-gene transcriptomic panel provided no statistically significant discrimination or calibration gain over baseline clinical staging alone in this matched benchmark.

---

### Osaka Theranostics 2027 Formal Decision: NO-GO
- **Official Program:** Focused exclusively on theranostics, radioligands ($^{68}\text{Ga}$-PSMA-11 / $^{177}\text{Lu}$-PSMA-617), biodistribution dosimetry, and radiomolecular therapeutics.
- **TCGA-PRAD:** 0% theranostic data. Maintaining scientific integrity requires not submitting genomic recurrence prediction to Osaka.
