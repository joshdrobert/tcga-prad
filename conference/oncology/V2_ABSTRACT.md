# ONCOLOGY INFORMATICS RESEARCH ABSTRACT
## FUTURE SUBMISSION (AUA / AMIA 2027)

**Title:** Endpoint Discordance, Predictor Chronology, and Missing Molecular Data in Prostate Cancer Recurrence Benchmarks  
**Author:** Joshua Robert [1]  
**Affiliation:** [1] Texas A&M School of Engineering Medicine, Houston, TX, USA  
**Target Venues:** American Urological Association (AUA 2027) / American Medical Informatics Association (AMIA 2027)  

---

### Abstract Body

#### Background
Public cancer registries combine heterogeneous clinical encounter forms and incomplete molecular profiles, complicating prognostic model evaluation. We audited longitudinal recurrence ascertainment, clinical predictor timing, and the incremental predictive value of transcriptomic profiling in post-radical prostatectomy patients.

#### Methods
We analyzed 314 patients in the TCGA prostate adenocarcinoma (TCGA-PRAD) cohort across 336 encounter records. We reconstructed a primary definitive binary endpoint restricted to documented positive recurrence (n=72) and documented double-negative recurrence-free cases (n=191), excluding 46 asymmetric surveillance records and five fully indeterminate cases. We audited clinical predictor chronology, separating baseline Day 0 peri-operative pathology (age, stage, Gleason sum) from post-operative surveillance PSA. On identical participants with complete transcriptomic profiles (n=188), we evaluated 5-fold cross-validated logistic regression with paired bootstrap differences in AUROC. Continuous-time recurrence-free survival was evaluated on strict sensitivity (n=257) and exploratory (n=303) cohorts.

#### Results
Comparing specific endpoints, biochemical recurrence agreed 100.0% (269/269) with cBioPortal, whereas disease-free survival agreed in 92.3% (241/261), uncovering 20 patients with XML-documented recurrence or secondary tumor events labeled disease-free in cBioPortal. Predictor timing analysis revealed XML surveillance PSA was measured at a median of 511 days post-surgery (271 of 276 values post-operative). In the primary paired molecular cohort (n=188, 45 events), pure baseline clinical staging achieved an AUROC of 0.693 (95% CI: 0.603–0.776), whereas adding the 13-gene expression panel yielded an AUROC of 0.648 (95% CI: 0.551–0.742; paired ΔAUROC: -0.046, 95% CI: [-0.124, +0.032], p=0.243). In strict continuous-time survival analysis (n=257, 66 dated events, 191 double-negative censored), out-of-fold Cox concordance was 0.647 (1-year RFS: 87.3%, 5-year RFS: 61.9%), compared with 0.693 when admitting 46 asymmetric records (n=303). Schoenfeld tests confirmed proportional hazards (all p > 0.40).

#### Conclusions
Failure to adjudicate endpoint definitions, distinguish baseline from surveillance biomarkers, or exclude indeterminate records distorts prognostic evaluation. The evaluated 13-gene transcriptomic panel provided no statistically significant discrimination gain over baseline clinical staging alone in this matched benchmark.

---

### Osaka Theranostics 2027 Formal Decision: NO-GO
- **Official Program:** Focused exclusively on theranostics, radioligands ($^{68}\text{Ga}$-PSMA-11 / $^{177}\text{Lu}$-PSMA-617), biodistribution dosimetry, and radiomolecular therapeutics.
- **TCGA-PRAD:** 0% theranostic data. Maintaining scientific integrity requires not submitting genomic recurrence prediction to Osaka.
