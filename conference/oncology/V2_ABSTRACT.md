# V2 ONCOLOGY / CLINICAL INFORMATICS ABSTRACT DRAFT

**Target Venues:** American Urological Association (AUA 2027) / AMIA Clinical Informatics Conference 2027  
**Track:** Machine Learning, Clinical Outcomes, and Biomarker Harmonization in Urologic Oncology  

---

### Abstract Title
Harmonizing Longitudinal Recurrence Endpoints and Quantifying Molecular Missingness in the TCGA-PRAD Radical Prostatectomy Cohort

### Authors
Joshua Robert [1]  
[1] Texas A&M School of Engineering Medicine, Houston, TX, USA  

---

### Background
Multi-omic cancer repositories such as The Cancer Genome Atlas (TCGA) prostate adenocarcinoma project (TCGA-PRAD) are widely utilized to benchmark molecular recurrence predictors. However, benchmark transportability is challenged by incomplete clinical follow-up, discordance between curated portals (e.g., cBioPortal) and primary GDC clinical XML supplements, and unmeasured missingness in molecular profiling. We audited longitudinal follow-up ascertainment, harmonized multidimensional recurrence endpoints, and evaluated the incremental prognostic value of bulk RNA-seq against clinical staging.

### Methods
In 314 consecutive patients undergoing radical prostatectomy with curative intent, clinical XML encounter records were reconstructed into an auditable longitudinal event ledger. A composite recurrence endpoint was adjudicated across biochemical recurrence (BCR), new tumor events (NTE), and post-hormone therapy progression. Cross-dataset concordance was evaluated against matched cBioPortal records. Predictive utility was assessed using 5-fold stratified cross-validation for clinical staging alone versus staging plus a 13-gene transcriptomic panel in both full ($n=314$) and paired complete-case ($n=228$, 86 unsequenced patients excluded) cohorts. Recurrence-free survival (RFS) was quantified using Kaplan-Meier and internally validated Cox proportional hazards regression.

### Results
Longitudinal follow-up was successfully ascertained for 311 of 314 patients (99.0%; median follow-up: 908 days). Composite recurrence occurred in 72 patients (22.9%), including 66 dated events and 6 undated clinical recurrences. Across 269 evaluable matched cases in cBioPortal, cross-check revealed 88.5% concordance and exactly 31 discordant cases, all reflecting clinical local recurrences or distant metastases uncaptured by cBioPortal's isolated serum PSA elevation indicator. In paired complete cases ($n=228$), adding transcriptomic features did not improve out-of-fold discrimination over clinical staging (AUROC 0.686 [0.589–0.764] vs 0.689 [0.598–0.765]) and worsened calibration slope from 0.777 to 0.540 under median imputation. On 308 evaluable survival patients, Kaplan-Meier RFS was 88.8% at 1 year, 75.1% at 3 years, and 64.8% at 5 years. Cox proportional hazards regression achieved a 5-fold cross-validated concordance index of 0.688 (apparent: 0.703).

### Conclusion
Secondary public cancer repositories exhibit substantial endpoint discordance and transcriptomic missingness that degrade model calibration when handled naively. Rigorous endpoint harmonization and paired complete-case evaluation demonstrate that standard clinical parameters remain the primary drivers of post-prostatectomy recurrence risk.
