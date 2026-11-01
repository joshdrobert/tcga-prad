# ITERATION 1 PROTOCOL: Impact of Longitudinal Endpoint Adjudication on Model Discrimination, Calibration, and Patient Risk Stratification

**Protocol ID:** ITER-01-ENDPOINT-IMPACT  
**Date & Timestamp:** 2026-11-01 15:33:00 CST  
**Investigators:** Principal Computational Oncology Methodologist & Reproducibility Lead  
**Target Venue:** 9th Theranostics World Congress (TWC 2027, Osaka) — Artificial Intelligence Track  

---

## 1. Research Hypothesis and Clinical Relevance
In localized prostate cancer following radical prostatectomy, accurately identifying patients at high risk of recurrence is critical for guiding salvage therapies and emerging theranostic interventions (such as PSMA-PET/CT restaging or radioligand therapy clinical trials). However, clinical data in public repositories (e.g., TCGA-PRAD) often contain contradictory longitudinal signals between single-encounter biochemical recurrence (`bcr_status`) and multi-encounter new tumor event tracking (`nte_status`).

We hypothesize that:
1. Naive single-source endpoint extraction misclassifies patient recurrence status and creates distorted time-to-event estimates compared to multi-source longitudinal adjudication.
2. Adjudicating these discrepancies alters cross-validated discrimination (AUROC, AUPRC) and calibration (Brier score, calibration slope) for both clinical and genomic models.
3. Patient-level risk ranking is sensitive to endpoint definitions, shifting individual assignment into illustrative high-risk tiers (e.g., the top risk quintile).
4. The incremental discrimination of the 13-gene RNA panel over clinical staging remains null ($\Delta\text{AUROC} \le 0$) regardless of adjudication tier, confirming that genomic panels cannot substitute for rigorous clinical staging or high-resolution imaging.

---

## 2. Literature Gap & Scientific Novelty
Published genomic and machine learning studies in oncology frequently take public clinical tables as ground truth without auditing multi-source XML forms or reconciling status-versus-date contradictions. In TCGA-PRAD, 29 patients possess documented conflicts between BCR status and numeric BCR dates, and 46 patients have asymmetric single-source non-event tracking. Prior reports have never systematically quantified how these specific data curation choices alter model performance, patient risk ranking, and apparent genomic utility on identical, matched patients.

---

## 3. Cohorts, Population Flow, and Time Origin
- **Total Evaluated Cohort:** $N = 314$ radical prostatectomy patients from TCGA-PRAD.
- **Matched Molecular Subcohort:** $N = 187$ complete cases with verified clinical staging and 13-gene RNA-seq expression.
- **Time Origin:** Days from initial pathologic diagnosis (day 0). All predictors (pathologic T stage, lymph node status, surgical margin status, primary/secondary Gleason patterns, age) are verified resection/pathologic features.
- **Endpoint Policies Compared on Identical Patients ($N=187$):**
  1. **Naive Raw-Status Reference:** Outcome defined solely by raw `bcr_status` = YES from initial reporting form (ignoring multi-encounter NTE and conflicts).
  2. **Permissive Adjudication:** Admitting numeric BCR dates and events even where single-form status was coded negative ($N=187$; 50 events, 137 non-events).
  3. **Strict Definitive Adjudication:** Requiring multi-source concordance between BCR and NTE; excluding unverified single-source non-events and non-prostate second primaries ($N=187$; 44 events, 143 non-events).

---

## 4. Models, Predictors, and Cross-Validation
- **Clinical Model:** Logistic regression with L2 regularization using standard post-prostatectomy pathologic features:
  - Pathologic T stage (pT2, pT3a, pT3b/pT4)
  - Regional lymph node status (pN0, pN1)
  - Surgical margin status (Negative, Positive)
  - Primary and secondary Gleason patterns (Gleason grade groups)
  - Patient age at diagnosis
- **Clinical + 13-Gene RNA Model:** Clinical features augmented with standardized expression values of the 13 preselected recurrence-associated genes.
- **Validation Design:** Stratified 5-fold cross-validation with a fixed seed (`seed=42`). Exactly identical fold splits are preserved across all endpoint comparison policies. All preprocessing, imputation, and scaling are strictly fit within training folds to eliminate data leakage.

---

## 5. Primary Effect Metrics & Uncertainty Estimation
- **Discrimination:** Out-of-fold (OOF) AUROC and AUPRC (Average Precision).
- **Calibration:** Brier score, Expected Calibration Error (ECE, 10 equal-frequency bins), and calibration slope.
- **Incremental Genomic Utility:** $\Delta\text{AUROC} = \text{AUROC}_{\text{Clinical+RNA}} - \text{AUROC}_{\text{Clinical}}$ and $\Delta\text{Brier}$.
- **Patient Risk Ranking Stability:** Spearman rank correlation ($\rho$) and Kendall tau ($\tau$) of predicted probabilities between endpoint policies.
- **Top-Quintile Risk Reclassification:** Proportion of patients whose assignment to the top risk quintile (highest 20% predicted risk) shifts between models or endpoint policies.
- **Uncertainty Approach:** 1,500-iteration patient-level clustered bootstrap with 95% percentile confidence intervals.

---

## 6. Stop Criteria & Safeguards
- Zero tolerance for data leakage (no post-treatment PSA, no cross-fold peeking).
- If incremental RNA $\Delta\text{AUROC}$ remains $\le 0$, it must be reported honestly as negative.
- Results must be labeled methodological groundwork for prostate theranostics/AI triage, explicitly disclaiming actual PET or radioligand therapy data.
