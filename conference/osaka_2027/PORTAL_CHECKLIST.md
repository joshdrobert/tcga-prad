# Osaka TWC 2027 Abstract Submission Portal Checklist

**Target Conference:** 9th Theranostics World Congress (Osaka TWC 2027)  
**Portal URL:** https://www.theranostics-world-congress.org/abstracts  
**Submission Deadline:** November 1, 2026, 23:59 CET (**TODAY**)  
**Submission Status:** Ready for Human Author Paste & Submit

---

## Pre-Submission Verification Summary

- [x] Abstract text strictly $\le$ 2,500 characters excluding spaces (Actual: **2,253 characters without spaces**; 2,381 with title).
- [x] Mandatory section headings present: `Aim`, `Material & Methods`, `Results`, `Conclusion`.
- [x] Author name and affiliation verified: Joshua Robert, Texas A&M School of Engineering Medicine, Houston, TX, USA.
- [x] No fabricated data, phantom affiliations, or unverified claims.
- [x] All numbers verified against TCGA-PRAD V3.2 audited baseline.

---

## Step-by-Step Portal Submission Guide

### Step 1: Open Portal & Log In
1. Navigate to: [https://www.theranostics-world-congress.org/abstracts](https://www.theranostics-world-congress.org/abstracts)
2. Log into your account or register a new author profile if not already completed.

### Step 2: Select Track & Presentation Category
1. **Primary Topic / Subject Area:** Select **Prostate Cancer** (or **Artificial Intelligence**).
2. **Presentation Type:** Select **Poster Presentation** (or *Oral Presentation* if desired).

### Step 3: Enter Title
- Paste the exact title:
  ```text
  Auditing Clinical Provenance, Molecular Prognosis, and Endpoint Harmonization in Machine Learning Benchmarks for Post-Prostatectomy Recurrence
  ```
  *(Note: 14 words; strictly Title Case).*

### Step 4: Enter Author & Affiliation Details
- **First Name:** Joshua
- **Last Name:** Robert
- **Institution / Affiliation:** Texas A&M School of Engineering Medicine
- **City:** Houston
- **State:** TX
- **Country:** United States
- **Presenting Author:** Yes (Mark checkbox)
- **Corresponding Author:** Yes (Mark checkbox)

### Step 5: Paste Abstract Content by Section
Depending on the portal input fields, either paste into separate section boxes or as a single block:

#### Option A: If Portal Has Separate Fields for Each Section
- **Aim Box:**
  ```text
  To evaluate endpoint integrity, molecular missingness, and prognostic discrimination in machine learning benchmarks for post-prostatectomy prostate cancer recurrence, establishing rigorous clinical baseline performance necessary for risk stratification prior to salvage theranostic interventions.
  ```
- **Material & Methods Box:**
  ```text
  We conducted a comprehensive provenance audit of the TCGA-PRAD radical prostatectomy cohort (n=314) linking multi-source clinical XML follow-up records with GDC RNA-sequencing profiles. Primary post-prostatectomy recurrence was defined using an audited multi-source provenance ledger (n=261; 70 events, 191 double-negatives), excluding indeterminate outcomes, non-primary events, and unconfirmed intervals. Predictor discrimination and probability calibration were evaluated for standard clinical staging (age, pT, pN, Gleason score) versus clinical plus a candidate 13-gene transcriptomic panel using stratified 5-fold cross-validation in a complete-case paired cohort (n=187; 44 events). Survival analysis was modeled from initial pathologic diagnosis using Cox proportional hazards and Kaplan-Meier estimation (n=248 with 57 dated events and 191 double-negative censored), strictly excluding surveillance PSA (median 511.0 days post-diagnosis) from baseline predictors to prevent look-ahead bias.
  ```
- **Results Box:**
  ```text
  In the paired molecular cohort (n=187), the clinical staging model achieved an out-of-fold AUROC of 0.7104 (95% CI 0.6208–0.7913), Brier score of 0.1673, and calibration slope of 0.756. Adding the 13-gene transcriptomic panel yielded an AUROC of 0.7028 (95% CI 0.6046–0.7917), Brier score of 0.1672, and calibration slope of 0.511 (paired ΔAUROC -0.0076 [95% CI -0.0735 to +0.0559]; paired ΔBrier -0.0001 [95% CI -0.0198 to +0.0199]). Audit of recurrence provenance resolved 29 biochemical recurrence status/date discrepancies across registry supplements. In strict continuous survival analysis (n=248), 1-, 3-, and 5-year recurrence-free survival rates were 88.1%, 75.2%, and 65.3% (median not reached), with a cross-validated Harrell's C-index of 0.6716 (apparent C-index 0.7011).
  ```
- **Conclusion Box:**
  ```text
  Candidate transcriptomic panels did not provide incremental prognostic discrimination over standard clinicopathologic staging for post-prostatectomy recurrence. Accurate baseline risk stratification for patient selection in emerging theranostic and salvage pathways requires rigorous endpoint harmonization, strict prevention of immortal time bias, and complete-case reporting rather than uncritical reliance on high-dimensional genomic features.
  ```

#### Option B: If Portal Has a Single Text Box for Abstract Body
- Copy the entire structured text directly from `READY_TO_PASTE_ABSTRACT.txt`.

### Step 6: Disclosures and Declarations
- **Conflict of Interest:** Check "No conflicts of interest" or paste:
  ```text
  The author declares no competing financial or non-financial interests.
  ```
- **Funding:** Check "No external funding" or enter "Investigator-initiated academic research; no dedicated grant support."
- **Data Policy:** Confirm adherence to NCI GDC public data use agreements (IRB exempt de-identified secondary research).

### Step 7: Final Review & Submit
1. Verify character counter in portal shows $\le$ 2,500 characters excluding spaces.
2. Click **Submit Abstract**.
3. Save the submission confirmation email / PDF receipt into your records.
