# Osaka TWC 2027 — Portal Entry Values & Submission Guide (Draft A-325)

**Conference:** 9th Theranostics World Congress (Osaka TWC 2027)  
**Portal URL:** https://www.theranostics-world-congress.org/abstracts  
**Submission Deadline:** November 1, 2026, 23:59 CET (**TODAY — 16:59 CST**)  
**Draft Tracking Number:** **A-325**  
**Action Required:** Human author must log in, review draft fields, and click **SUBMIT**.

---

## 1. Character & Word Count Compliance Audit

| Metric | Online Portal Rule | Public Website Rule | Measured Value | Margin / Safety Headroom |
| :--- | :--- | :--- | :--- | :--- |
| **Abstract Body Characters WITH Spaces** | $\le$ 3,000 chars | N/A | **2,631 chars** | **PASS (+369 chars headroom)** |
| **Abstract Body Characters WITHOUT Whitespace** | N/A | $\le$ 2,500 chars | **2,297 chars** | **PASS (+203 chars headroom)** |
| **Full Text (with Title) WITHOUT Whitespace** | N/A | $\le$ 2,500 chars | **2,425 chars** | **PASS (+75 chars headroom)** |
| **Word Count** | Informational | Informational | **332 words** | Compliant |
| **Title Word Count** | Informational | Informational | **15 words** | Compliant (Title Case) |
| **Forbidden Content Audit** | No commercial names, URLs, or emails in text | **PASS (0 found)** | Compliant |

---

## 2. Portal Tab-by-Tab Entry Values

### Tab 1: Basic Information
- **Title (Exact, 15 words):**
  ```text
  Auditing Clinical Provenance, Molecular Prognosis, and Endpoint Harmonization in Machine Learning Benchmarks for Post-Prostatectomy Recurrence
  ```
- **Category / Topic Track:**
  - **Recommended Primary:** **Artificial Intelligence (AI)**
  - *Secondary / Alternative:* **Prostate Cancer, including PSMA-Targeted Theranostics**
  - *Reviewer-Fit Rationale:* The research is a rigorous statistical machine-learning benchmark and clinical registry provenance audit. It does **not** contain radiotracer synthesis, PSMA-PET/CT lesion segmentation, or radioligand therapy (RLT) outcomes. Selecting `Artificial Intelligence (AI)` is the most scientifically truthful categorization. If the portal forces disease-first indexing, `Prostate Cancer, including PSMA-Targeted Theranostics` is acceptable because post-prostatectomy recurrence risk stratification directly defines the clinical threshold for restaging and salvage theranostics.
- **Presentation Preference:**
  - Select **Poster Presentation** (or *Oral or Poster* if offered, targeting poster).

---

### Tab 2: Add Authors
- **First Name:** `Joshua`
- **Last Name:** `Robert`
- **Institution / Affiliation:** `Texas A&M School of Engineering Medicine`
- **City:** `Houston`
- **State:** `TX`
- **Country:** `United States`
- **Email:** `[AUTHOR CONFIRMATION REQUIRED — Enter primary institutional email in portal]`
- **Presenting Author:** `Yes` (Check box)
- **Corresponding Author:** `Yes` (Check box)
- **Co-Author Permissions:** None (Single-author investigator study).

---

### Tab 3: Abstract Details & Disclosures (AUTHOR CONFIRMATION REQUIRED)
- **Conflict of Interest (COI):**
  - *Action:* Confirm whether any commercial relationships, consultancies, or competing financial interests exist.
  - *Default text if no conflicts exist:* `"The author declares no competing financial or non-financial interests."`
- **Funding Statement:**
  - *Action:* Confirm whether any departmental, institutional, or grant funding applies.
  - *Default text if unfunded:* `"Investigator-initiated academic study; no dedicated grant support."`
- **Ethical Oversight & IRB Status:**
  - *Action:* Declare that this is a secondary analysis of public, de-identified genomic data (NCI GDC TCGA-PRAD). Note that access to public de-identified datasets does not by itself document an institutional IRB determination; confirm Texas A&M institutional policies.
- **Congress Registration Notice:**
  - *Portal vs. Website Policy Note:* The online portal may display a warning that registration is required to submit, whereas the public website states that presenting author registration/payment is required upon abstract acceptance. **Do not pay registration prematurely unless the portal strictly blocks draft submission without an invoice number.** If required by the portal UI, follow portal payment prompts; otherwise, submit the abstract and register upon acceptance.

---

### Tab 4: Abstract Text

#### Option A: Individual Section Fields
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

#### Option B: Combined Text Box
- Copy directly from [`tcga-prad/conference/osaka_2027/PORTAL_READY_2026-11-01.txt`](file:///Users/josh/Research/March2027/tcga-prad/conference/osaka_2027/PORTAL_READY_2026-11-01.txt).

---

### Tab 5: Summary & Final Submission Checklist (Draft A-325)
1. [ ] Log in to `https://www.theranostics-world-congress.org/abstracts`.
2. [ ] Open existing Draft **A-325**.
3. [ ] Confirm Title, Category (`Artificial Intelligence (AI)` or `Prostate Cancer`), and Author details.
4. [ ] Paste and verify Abstract Text.
5. [ ] Confirm character count counter displays within portal bounds ($\le$ 3,000 characters with spaces, $\le$ 2,500 without spaces).
6. [ ] Complete COI, Funding, and IRB confirmations.
7. [ ] Click **SUBMIT ABSTRACT** (Ensure status transitions from **DRAFT** to **SUBMITTED**).
8. [ ] Verify receipt of confirmation email and store final submission PDF.
