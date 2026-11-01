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
We analyzed 314 patients in the TCGA prostate adenocarcinoma (TCGA-PRAD) cohort across 336 encounter records. We reconstructed a primary definitive prostate cancer recurrence endpoint restricted to documented positive recurrence (n=70) and documented double-negative recurrence-free cases (n=191), separating two non-prostatic second primary malignancies and excluding 46 asymmetric surveillance records and five fully indeterminate cases. We audited clinical predictor chronology, establishing initial pathologic diagnosis (core needle biopsy Day 0 for 97.1% of patients [305/314; 98.1% of 311 patients with non-missing diagnostic method]) as the follow-up time origin while recognizing that surgical pathology (Gleason sum, pT, pN) and transcriptomic profiles are conditioned on radical prostatectomy resection with an unmeasured biopsy-to-surgery interval. Postoperative surveillance PSA was separated from baseline staging. On identical participants with complete transcriptomic profiles (n=187, 44 events), we evaluated 5-fold cross-validated logistic regression with paired bootstrap differences in AUROC, Brier score, and calibration slope. Longitudinal time-to-recurrence provenance was audited across 29 patients with status/date discrepancies; continuous-time recurrence-free survival from diagnosis was evaluated on strict primary sensitivity (n=248, 57 dated events) and exploratory (n=294, 57 dated events) cohorts using baseline Cox models without post-baseline surveillance PSA leakage, alongside permissive sensitivity benchmarking (n=255 and n=301, 64 dated events).

#### Results
Comparing specific endpoints, biochemical recurrence agreed 100.0% (269/269) with cBioPortal, whereas disease-free survival agreed in 92.3% (239/259), uncovering 20 patients with XML-documented recurrence labeled disease-free in cBioPortal. Predictor timing analysis revealed XML surveillance PSA was measured at a median of 511 days post-diagnosis (271 of 276 values post-baseline). In the primary paired molecular cohort (n=187, 44 events), pure baseline clinical staging achieved an AUROC of 0.710 (95% CI: 0.621–0.791), Brier score of 0.1673 (95% CI: 0.1366–0.2018), and calibration slope of 0.756, whereas adding the 13-gene expression panel yielded an AUROC of 0.703 (95% CI: 0.605–0.792), Brier score of 0.1672 (95% CI: 0.1329–0.2057), and calibration slope of 0.511 (paired ΔAUROC: -0.008, 95% CI: [-0.073, +0.056]; paired ΔBrier: -0.0001, 95% CI: [-0.0198, +0.0199]). In continuous-time survival analysis from initial diagnosis, median recurrence-free survival was not reached in any cohort. Under strict primary adjudication withholding unverified dates (n=248, 57 dated events, 191 double-negative censored), out-of-fold Cox concordance was 0.672 (apparent 0.701; 1-year RFS: 88.1%, 5-year RFS: 65.3%), compared with 0.670 out-of-fold when admitting 46 asymmetric records (n=294). In permissive sensitivity analysis admitting 7 contradictory BCR dates without independent verification (n=255, 64 events), out-of-fold Cox concordance was 0.649 (1-year RFS: 87.2%, 5-year RFS: 62.7%), and 0.646 admitting asymmetric records (n=301). Schoenfeld residual tests did not detect a statistically significant violation of proportional hazards (p > 0.30 for all baseline covariates).

#### Conclusions
Failure to adjudicate endpoint definitions, separate second primary malignancies, resolve contradictory event provenance, or exclude post-baseline surveillance biomarkers distorts prognostic evaluation. The evaluated 13-gene transcriptomic panel provided no statistically significant discrimination or Brier score improvement over baseline clinical staging alone in this matched benchmark, while compressing calibration slope from 0.756 to 0.511. Strict event provenance resolution removes temporal noise and improves continuous-time survival concordance.

---

### Osaka Theranostics 2027 Formal Decision: NO-GO
- **Official Program:** Focused exclusively on theranostics, radioligands ($^{68}\text{Ga}$-PSMA-11 / $^{177}\text{Lu}$-PSMA-617), biodistribution dosimetry, and radiomolecular therapeutics.
- **TCGA-PRAD:** 0% theranostic data. Maintaining scientific integrity requires not submitting genomic recurrence prediction to Osaka.
