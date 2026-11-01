# QUANTITATIVE RESULTS COMPARISON: NAIVE VS HARMONIZED RECURRENCE ENDPOINTS

**Study Cohort:** TCGA-PRAD Radical Prostatectomy ($N=314$)  
**Matched Complete-Case Molecular Analysis:** $N=187$ patients with verified 13-gene RNA sequencing  
**Cross-Validation:** 5-Fold Stratified Cross-Validation (Seed 42, Identical Patient Folds)  
**Bootstrap Uncertainty:** 1,500 patient-level resamples (95% Percentile Confidence Intervals)  

---

## 1. Summary Comparison Table

| Metric / Dimension | Naive Unadjudicated Endpoint (Single-Form BCR) | Multi-Source Harmonized Endpoint (BCR + NTE) | Absolute Difference ($\Delta$) | Interpretation / Significance |
| :--- | :--- | :--- | :--- | :--- |
| **Recurrence Events** | 25 events (13.4%) | 44 events (23.5%) | +19 events (+76%) | 43.2% of recurrences missed naively |
| **Non-Events (Censored)** | 162 non-events | 143 non-events | -19 non-events | False-negative cases corrected |
| **Pure Clinical AUROC** | 0.6368 (95% CI: 0.5067–0.7603) | 0.7104 (95% CI: 0.6208–0.7913) | +0.0736 | Clinical staging restored |
| **Clinical + 13-Gene RNA AUROC** | 0.7536 (95% CI: 0.6162–0.8651) | 0.7028 (95% CI: 0.6046–0.7917) | -0.0508 | Genomic overfit corrected |
| **Incremental $\Delta\text{AUROC}$** | **+0.1168 (95% CI: +0.0229 to +0.2227)** | **-0.0076 (95% CI: -0.0763 to +0.0566)** | **-0.1244** | **Spurious gain eliminated** |
| **Clinical Calibration Slope** | 0.424 | 0.756 | +0.332 | Calibration substantially improved |
| **Molecular Calibration Slope** | 0.601 | 0.511 | -0.090 | Overconfident calibration |
| **Clinical Brier Score** | 0.1179 | 0.1673 | +0.0494 | Higher baseline event prevalence |
| **Molecular Brier Score** | 0.1020 | 0.1672 | +0.0652 | No Brier benefit over clinical |
| **Top-Quintile Risk Movement**| Baseline | 18 / 187 patients shift (9.6%) | 18 patients | Rank $\rho = 0.894, \tau = 0.710$ |

---

## 2. Clinical Profile of the 19 Discordant False-Negative Cases

| Subgroup | Sample Size ($N$) | Pathologic Stage T3/T4 (%) | Gleason Score $\ge 8$ (%) | Mean Gleason Score | Recurrence Follow-Up (Days) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Missed Cases (Adjudicated=1, Naive=0)** | **19** | **89.5% (17 / 19)** | **63.2% (12 / 19)** | **8.16** | Median 925d (Range: 57–2036d) |
| **Concordant Positives (Both=1)** | 25 | 88.0% (22 / 25) | 68.0% (17 / 25) | 8.20 | Median 692d (Range: 85–1634d) |
| **Concordant Negatives (Both=0)** | 143 | 52.4% (75 / 143) | 31.5% (45 / 143) | 7.39 | Median 754d (Range: 32–3265d) |

---

## 3. Decision Curve Analysis (DCA) Net Benefit

| Threshold Probability ($p_t$) | Treat None | Treat All | Pure Clinical Staging | Clinical + 13-Gene RNA | Delta Net Benefit (Molecular - Clinical) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **0.10** | 0.0000 | 0.1503 | **0.1551** | 0.1545 | -0.0006 |
| **0.15** | 0.0000 | 0.1003 | **0.1362** | 0.1051 | -0.0311 (Clinical superior) |
| **0.20** | 0.0000 | 0.0441 | **0.1003** | 0.0976 | -0.0027 (Clinical superior) |
| **0.25** | 0.0000 | -0.0196 | 0.0713 | **0.0856** | +0.0143 |
| **0.30** | 0.0000 | -0.0924 | 0.0535 | **0.0726** | +0.0191 |
| **0.35** | 0.0000 | -0.1765 | **0.0403** | 0.0354 | -0.0049 (Clinical superior) |
