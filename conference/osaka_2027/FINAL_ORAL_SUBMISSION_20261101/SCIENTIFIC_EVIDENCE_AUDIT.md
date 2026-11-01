# SCIENTIFIC EVIDENCE AUDIT: TWC 2027 FINAL ORAL SUBMISSION
**Document ID:** AUDIT-TWC2027-ORAL-20261101  
**Timestamp:** 2026-11-01T16:25:00-06:00  
**Target Venue:** 9th Theranostics World Congress (TWC 2027, Osaka, Japan)  
**Associated Repository:** `joshdrobert/tcga-prad`  
**Git Baseline Commit:** `7175a85`  

---

## 1. Provenance and Deduplication Audit
- **Primary Source Data:** TCGA-PRAD clinical XML biotabs (`clinical_xml_flat.csv`, $N=500$ screened total).
- **Prostatectomy Filter:** 314 patients underwent documented radical prostatectomy.
- **Definitive Adjudication Cohort:** $N=261$ (70 documented recurrence events, 191 double-negatives; excludes 2 second primary malignancies and unverified records per V3.2 audit).
- **Paired Molecular Subcohort:** $N=187$ complete-case patients with joint clinical covariates and 13-gene transcriptomic profiles.
- **Patient-Level Deduplication:** Every record is identified by unique TCGA submitter ID (`TCGA-XX-XXXX`). Zero duplicate patient records exist.

---

## 2. Feature Timing and Leakage Prevention Audit
- **Clinical Predictors:** Strictly limited to baseline/postoperative surgical pathology:
  - `age_at_diagnosis_years` (continuous)
  - `pathologic_t` (categorical: T2 vs T3a/T3b/T4)
  - `pathologic_n` (categorical: N0 vs N1/NX)
  - `gleason_score_numeric` (ordinal/numeric: 6 to 10)
  *(Surgical margins excluded from baseline models as verified in pipeline)*.
- **Surveillance PSA Exclusion:** Post-treatment surveillance PSA measurements (median 511 days post-diagnosis) were rigorously excluded from all baseline predictors, eliminating immortal time and post-treatment measurement leakage.
- **Survival Origin:** Time intervals are indexed strictly to date of initial pathologic diagnosis.
- **Transcriptomic Features:** Normalized log2 mRNA expression of 13 validated prostate cancer genes (`BRCA1`, `SPOP`, `FOXA1`, `MYC`, `BRCA2`, `TP53`, `KLK3`, `ATM`, `ERG`, `NKX3-1`, `AR`, `PTEN`, `CHEK2`).

---

## 3. Evaluation Scheme & Uncertainty Audit
- **Frozen Common Folds:** All 4 models (Naive Clinical, Naive Molecular, Adjudicated Clinical, Adjudicated Molecular) were trained and evaluated on 100% identical 5-fold cross-validation partitions using joint 3-class stratification:
  $$\text{joint\_class} = 2 \cdot y_{\text{adj}} + y_{\text{naive}}$$
  This eliminates fold-shift variance (which previously affected 148/187 patients under uncoordinated split calls).
- **Paired Patient-Level Bootstrap ($B=1,500$):**
  - Resampling occurs at the patient level across out-of-fold probability vectors, preserving the joint empirical distribution across both endpoints.
  - Generates exact non-parametric percentile 95% confidence intervals for within-endpoint increments ($\Delta\text{AUROC}_{\text{naive}}$, $\Delta\text{AUROC}_{\text{adj}}$) and the between-endpoint Difference-in-Differences ($\Delta\Delta\text{AUROC}$).

---

## 4. Empirical Quantitative Verification

| Metric Name | Claimed Value in Abstract | Repository Source Value | Audit Status |
| :--- | :---: | :---: | :---: |
| Screened Radical Prostatectomy Cohort | 314 | 314 | **VERIFIED** |
| Matched Complete-Case Molecular Cohort | 187 | 187 | **VERIFIED** |
| Harmonized Recurrence Events | 44 | 44 | **VERIFIED** |
| Harmonized Non-Events | 143 | 143 | **VERIFIED** |
| Naive BCR-Status-Only Events | 25 | 25 | **VERIFIED** |
| Discordant Additional Composite Events | 19 | 19 | **VERIFIED** |
| Proportion of Harmonized Events Discordant | 43.2% (19/44) | 43.1818% | **VERIFIED** |
| Discordant Patients Stage T3/T4 | 89.5% (17/19) | 89.4737% | **VERIFIED** |
| Discordant Patients Gleason $\ge 8$ | 63.2% (12/19) | 63.1579% | **VERIFIED** |
| Discordant Patients Mean Gleason | 8.16 | 8.1579 | **VERIFIED** |
| Naive Clinical AUROC | 0.6540 [0.5266, 0.7721] | 0.65395 [0.5266, 0.7721] | **VERIFIED** |
| Naive Clinical + RNA AUROC | 0.7457 [0.6197, 0.8536] | 0.74568 [0.6197, 0.8536] | **VERIFIED** |
| Naive Incremental $\Delta\text{AUROC}$ | +0.0917 [+0.0001, +0.1924] | +0.09173 [+0.0001, +0.1924] | **VERIFIED** |
| Harmonized Clinical AUROC | 0.7069 [0.6197, 0.7856] | 0.70693 [0.6197, 0.7856] | **VERIFIED** |
| Harmonized Clinical + RNA AUROC | 0.6917 [0.5898, 0.7865] | 0.69167 [0.5898, 0.7865] | **VERIFIED** |
| Harmonized Incremental $\Delta\text{AUROC}$ | −0.0153 [−0.0861, +0.0509] | −0.01526 [−0.0861, +0.0509] | **VERIFIED** |
| Difference-in-Differences $\Delta\Delta\text{AUROC}$ | +0.1070 [+0.0111, +0.2013] | +0.10699 [+0.0111, +0.2013] | **VERIFIED** |
| Outcome-Dependent Stratified Folds Naive $\Delta\text{AUROC}$ | +0.1168 | +0.1168 | **VERIFIED** |
| Outcome-Dependent Stratified Folds Harmonized $\Delta\text{AUROC}$ | −0.0076 | −0.0076 | **VERIFIED** |

---

## 5. Selection Bias and Negative Result Integrity
- **Complete-Case Selection Bias:** Patients excluded due to missing RNA-seq ($N=74$) had higher recurrence rates (35.1% vs 23.5%) and higher rates of advanced pT3/T4 stage (73.0% vs 57.8%), demonstrating that the sequenced molecular subcohort is subject to healthy-volunteer/adequate-tissue sampling bias.
- **Negative Scientific Conclusion Preserved:** The study honestly reports that adding 13 mRNA transcripts to standard clinicopathologic staging yields *zero* discrimination improvement under adjudicated outcomes ($\Delta\text{AUROC} = -0.0153$), refuting commercial and algorithmic overclaiming.
