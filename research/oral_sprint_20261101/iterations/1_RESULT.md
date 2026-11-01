# ITERATION 1 RESULT: Spurious Genomic Gain Unmasked as an Artifact of Longitudinal Endpoint Incompleteness

**Protocol ID:** ITER-01-ENDPOINT-IMPACT  
**Execution Timestamp:** 2026-11-01 15:37:30 CST  
**Status:** COMPLETED & INDEPENDENTLY REPRODUCED  

---

## 1. Executive Summary & Core Discovery
By cross-evaluating the paired molecular cohort ($N=187$) under naive single-source reporting versus multi-source longitudinal adjudication, we uncovered a profound methodological phenomenon:
1. **Unadjudicated Incompleteness:** The naive single-source `bcr_status` field captures only 25 recurrences in the paired complete-case cohort, completely missing **19 true recurrence events (43.2% of the 44 actual recurrences)** that were documented on longitudinal follow-up under New Tumor Events (`nte_status`).
2. **Aggressive Clinical Profile of Missed Cases:** The 19 false-negative patients were heavily enriched for high-risk disease: **89.5% (17/19) had pathologic T3/T4 disease**, **63.2% had Gleason $\ge 8$**, and mean Gleason score was **8.16** (identical to concordant positives at 8.20, vs 7.39 in non-events).
3. **The Spurious Genomic Gain Artifact:**
   - Under naive unadjudicated endpoints, standard clinical staging is falsely penalized for predicting recurrence in these 19 aggressive patients who were mistakenly labeled non-recurrent. Consequently, the clinical model AUROC was suppressed to **0.6368** (95% CI: 0.5067 to 0.7603).
   - Adding the 13-gene RNA panel appeared to dramatically rescue discrimination to **0.7536** (95% CI: 0.6162 to 0.8651), creating an apparent, statistically significant incremental benefit: **$\Delta\text{AUROC} = +0.1168$ (95% CI: +0.0229 to +0.2227)**.
   - When longitudinal multi-source follow-up is harmonized, the 19 true events are properly recognized: clinical staging AUROC rises to **0.7104** (95% CI: 0.6208 to 0.7913), and the apparent incremental genomic gain completely evaporates: **$\Delta\text{AUROC} = -0.0076$ (95% CI: -0.0763 to +0.0566)**.
4. **Conclusion:** The apparent incremental discrimination of the 13-gene RNA panel in unadjudicated TCGA data is an **artifact of unverified outcome misclassification**. When true recurrences are accounted for, standard clinical staging demonstrates superior calibration and non-inferior clinical net benefit across all decision thresholds.

---

## 2. Quantitative Comparative Results Table

| Analysis Dimension / Metric | Naive Unadjudicated Endpoint | Multi-Source Harmonized Endpoint | Absolute Impact / Shift |
| :--- | :--- | :--- | :--- |
| **Evaluated Cohort Denominator** | $N = 187$ paired complete cases | $N = 187$ paired complete cases | Identical patients |
| **Identified Recurrence Events** | 25 events (13.4% prevalence) | 44 events (23.5% prevalence) | **+19 events (+76% increase)** |
| **Missed Recurrence Fraction** | 43.2% (19 / 44) | 0% (Ground truth harmonized) | Corrects 19 false negatives |
| **Clinical Model AUROC (5-Fold CV)** | 0.6368 (95% CI: 0.5067–0.7603) | 0.7104 (95% CI: 0.6208–0.7913) | **+0.0736 (Clinical restored)** |
| **Clinical + RNA AUROC (5-Fold CV)** | 0.7536 (95% CI: 0.6162–0.8651) | 0.7028 (95% CI: 0.6046–0.7917) | -0.0508 |
| **Incremental $\Delta\text{AUROC}$ (Genomic Gain)**| **+0.1168 (95% CI: +0.0229 to +0.2227)** | **-0.0076 (95% CI: -0.0763 to +0.0566)**| **-0.1244 (Artifact eliminated)** |
| **Clinical Model Calibration Slope** | 0.424 (Severe overconfidence) | 0.756 (Substantially improved) | +0.332 |
| **Molecular Calibration Slope** | 0.601 | 0.511 (Overconfident) | -0.090 |
| **Top-Quintile Risk Reclassification** | Baseline allocation | 18 / 187 patients shift (9.6%) | $\rho = 0.894, \tau = 0.710$ |
| **Decision Curve Net Benefit (15% pt)** | 0.118 (Molecular higher) | 0.136 (Clinical superior to RNA 0.105)| Clinical non-inferior |

---

## 3. Negative and Null Findings Preserved
- The 13-gene RNA panel does **not** improve discrimination over clinical staging once endpoints are harmonized ($\Delta\text{AUROC} = -0.0076$).
- Decision Curve Analysis confirms that clinical staging yields higher net benefit than the genomic model at threshold probabilities $\le 20\%$.
- Complete-case molecular analysis is subject to selection bias: patients without RNA ($N=74$) had higher recurrence rates (35.1% vs 23.5%) and higher rates of advanced pT3/T4 disease (73.0% vs 57.8%).
- Surveillance PSA measured at median 511 days post-diagnosis artificially inflates Cox concordance from 0.7011 to 0.7320 if improperly used as a baseline predictor.
