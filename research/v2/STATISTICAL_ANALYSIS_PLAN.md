# TCGA-PRAD V2 STATISTICAL ANALYSIS PLAN (SAP)

## 1. Study Population & Sample Size
- **Cohort:** 314 consecutive patients diagnosed with prostate adenocarcinoma undergoing radical prostatectomy with curative intent in TCGA-PRAD.
- **Analysis Populations:**
  1. *Full Cohort ($n=314$):* Complete surgical cohort evaluated for primary clinical staging risk discrimination.
  2. *Paired Complete-Case Cohort ($n=228$):* Subset with non-null normalized RNA-seq transcriptomic profiles (86 unsequenced patients excluded).
  3. *High-Confidence Outcome Cohort ($n=263$):* Documented recurrence ($n=72$) vs documented recurrence-free across both BCR and NTE axes ($n=191$), excluding 51 indeterminate or asymmetric follow-up cases.
  4. *Continuous-Time Survival Cohort ($n=308$):* Patients with verified positive event onset dates ($n=66$) or non-zero recurrence-free follow-up ($n=242$), strictly excluding the 6 undated recurrence cases from continuous-time models.

## 2. Statistical Metrics & Validation Design
- **Cross-Validation:** 5-fold stratified cross-validation (`StratifiedKFold`, seed = `20260626`) applied with strict fold-contained preprocessors (imputation, scaling, one-hot encoding).
- **Discrimination:** Area under the ROC Curve (AUROC) and Area under the Precision-Recall Curve (AUPRC), with 95% percentile bootstrap confidence intervals (1,500 iterations).
- **Calibration:** Calibration slope and intercept fitted via logistic calibration on out-of-fold logit probabilities. Brier score loss.
- **Survival Modeling:**
  - Kaplan-Meier estimator for Recurrence-Free Survival at 1, 3, and 5 years post-surgery.
  - Multivariable Cox Proportional Hazards regression for age, Gleason score, and pre-treatment PSA.
  - Reporting both apparent concordance index and 5-fold cross-validated out-of-fold C-index.
