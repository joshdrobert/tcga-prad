# Osaka TWC 2027 Abstract Submission Portal Checklist

**Target Conference:** 9th Theranostics World Congress (Osaka TWC 2027)  
**Portal URL:** https://www.theranostics-world-congress.org/abstracts  
**Submission Deadline:** November 1, 2026, 23:59 CET (**TODAY — 4:59 PM Central**)  
**Submission Status:** Ready for Human Author Verification & Submission

---

## Pre-Submission Verification Summary

- [x] Abstract body strictly $\le$ 2,500 characters excluding spaces (Actual: **2,297 characters without spaces**; 2,425 with title).
- [x] Mandatory section headings present: `Aim:`, `Material & Methods:`, `Results:`, `Conclusion:`.
- [x] Title word count verified: **15 words** in Title Case.
- [x] All numbers verified against TCGA-PRAD V3.2 audited baseline.
- [x] Potential future theranostic relevance clearly stated; no false claim of validated treatment selection.
- [x] Exclusion of surveillance PSA correctly framed as reducing temporal leakage; unmeasured prostatectomy dates acknowledged.
- [x] Statistically accurate panel finding: 13-gene panel did not demonstrate an AUROC improvement in this matched dataset.
- [ ] **Author Disclosures (COI, Funding, IRB Determination) — AUTHOR CONFIRMATION REQUIRED in portal.**

---

## Step-by-Step Portal Submission Guide

### Step 1: Open Portal & Log In
1. Navigate to: [https://www.theranostics-world-congress.org/abstracts](https://www.theranostics-world-congress.org/abstracts)
2. Log into your account or complete author registration.

### Step 2: Track & Presentation Selection
1. **Primary Topic / Subject Area:** Select **Prostate Cancer, including PSMA-Targeted Theranostics** (or **Artificial Intelligence**).
2. **Presentation Type:** Select **Poster Presentation** (or *Oral Presentation* if desired).

### Step 3: Enter Title
- Paste exact title (15 words, Title Case):
  ```text
  Auditing Clinical Provenance, Molecular Prognosis, and Endpoint Harmonization in Machine Learning Benchmarks for Post-Prostatectomy Recurrence
  ```

### Step 4: Author Information & Contact (AUTHOR ACTION REQUIRED)
- **Name:** Joshua Robert
- **Affiliation:** Texas A&M School of Engineering Medicine
- **City / Country:** Houston, TX, USA
- **Email:** `[Enter your verified institutional email]`
- **Presenting / Corresponding Author:** Select Yes

### Step 5: Disclosures & Regulatory Affirmations (AUTHOR ACTION REQUIRED)
- **Conflict of Interest:**
  - Select "No conflict of interest" or enter appropriate declaration if commercial relationships exist.
- **Funding Statement:**
  - Confirm whether study received institutional or grant support.
- **IRB / Ethics Statement:**
  - Affirm secondary analysis of public de-identified data per Texas A&M institutional policies.

### Step 6: Paste Abstract Content by Section
Depending on whether the portal uses separate fields or a single box:

#### If Portal Has Separate Section Fields:
- **Aim:**
  ```text
  To evaluate endpoint integrity, molecular missingness, and prognostic discrimination in machine learning benchmarks for post-prostatectomy prostate cancer recurrence, establishing an audited clinicopathologic baseline to inform future risk stratification prior to potential salvage evaluation.
  ```
- **Material & Methods:**
  ```text
  We conducted a provenance audit of the TCGA-PRAD radical prostatectomy cohort (n=314) linking multi-source clinical XML records with GDC RNA-sequencing profiles. Primary post-prostatectomy recurrence was defined using an audited provenance ledger (n=261; 70 events, 191 double-negatives), excluding indeterminate outcomes, non-primary events, and unconfirmed intervals. Predictor discrimination and calibration were evaluated for standard clinical staging (age, pT, pN, Gleason score) versus clinical plus a candidate 13-gene transcriptomic panel using stratified 5-fold cross-validation in a complete-case paired cohort (n=187; 44 events). Survival was modeled from initial pathologic diagnosis using Cox proportional hazards and Kaplan-Meier estimation (n=248; 57 dated events, 191 censored). Postoperative surveillance PSA (median 511.0 days post-diagnosis) was excluded to reduce temporal leakage; however, exact prostatectomy dates remain unavailable and post-resection pathology is modeled with diagnosis-indexed follow-up.
  ```
- **Results:**
  ```text
  In the paired molecular cohort (n=187), the clinical staging model achieved an out-of-fold AUROC of 0.7104 (95% CI 0.6208–0.7913), Brier score of 0.1673, and calibration slope of 0.756. Adding the 13-gene transcriptomic panel yielded an AUROC of 0.7028 (95% CI 0.6046–0.7917), Brier score of 0.1672, and calibration slope of 0.511 (paired ΔAUROC -0.0076 [95% CI -0.0735 to +0.0559]; paired ΔBrier -0.0001 [95% CI -0.0198 to +0.0199]). Audit of recurrence provenance resolved 29 biochemical recurrence status/date discrepancies across registry supplements. In strict continuous survival analysis (n=248), 1-, 3-, and 5-year recurrence-free survival rates were 88.1%, 75.2%, and 65.3% (median not reached), with a cross-validated Harrell's C-index of 0.6716 (apparent C-index 0.7011).
  ```
- **Conclusion:**
  ```text
  The evaluated 13-gene panel did not demonstrate an AUROC improvement over standard clinicopathologic staging in this matched dataset. While this retrospective study did not evaluate PSMA imaging, radioligand therapy, or treatment response, transparent endpoint harmonization, mitigation of temporal leakage, and complete-case reporting represent essential methodological prerequisites before recurrence benchmarks can inform future patient selection in prostate theranostics.
  ```

#### If Portal Has a Single Text Box:
- Paste the full text directly from [`conference/osaka_2027/READY_TO_PASTE_ABSTRACT.txt`](file:///Users/josh/Research/March2027/tcga-prad/conference/osaka_2027/READY_TO_PASTE_ABSTRACT.txt).

### Step 7: Final Review & Confirmation
1. Verify the character counter in the portal shows $\le$ 2,500 characters.
2. Complete submission before **November 1, 2026, 23:59 CET (4:59 PM Central)**.
3. Save confirmation email and abstract tracking ID.
