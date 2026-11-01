# 9th Theranostics World Congress (Osaka TWC 2027) — Official Submission Packet

**Conference:** 9th Theranostics World Congress (TWC 2027)  
**Location:** Osaka, Japan  
**Official Portal:** https://www.theranostics-world-congress.org/abstracts  
**Submission Deadline:** November 1, 2026, 23:59 CET (**TODAY — 4:59 PM Central**)  
**Project Repository:** https://github.com/joshdrobert/tcga-prad  
**Live Research Site:** https://joshdrobert.github.io/tcga-prad/  
**Scientific Baseline:** TCGA-PRAD V3.2 Provenance Reconciliation Baseline (`483ab18`)  
**Overall Readiness Classification:** **READY WITH VENUE-FIT RISK (AUTHOR REVIEW REQUIRED)**

---

## 1. Submission Metadata & Author Details

- **Title:** Auditing Clinical Provenance, Molecular Prognosis, and Endpoint Harmonization in Machine Learning Benchmarks for Post-Prostatectomy Recurrence
- **Primary Topic / Track:** Prostate Cancer, including PSMA-Targeted Theranostics *(Alternative: Artificial Intelligence)*
- **Presentation Type Preference:** Poster Presentation *(Oral if invited)*
- **Author Information:**
  - **Presenting / Corresponding Author:** Joshua Robert
  - **Affiliation:** Texas A&M School of Engineering Medicine, Houston, TX, USA
  - **Email:** `[AUTHOR CONFIRMATION REQUIRED]` *(Enter primary institutional email upon login)*
  - **Co-Author Permissions:** None (Single author; verify with home institution if co-investigators need listing)

---

## 2. Character & Word Count Compliance Verification

| Metric | Measured Value | Portal Limit | Status |
| :--- | :--- | :--- | :--- |
| **Title Word Count** | **15 words** | N/A (Title Case recommended) | **PASS** |
| **Title Character Count (excluding spaces)** | 128 characters | N/A | Informational |
| **Abstract Body Characters (excluding spaces)** | **2,297** | $\le$ 2,500 characters | **PASS (Under limit by 203 chars)** |
| **Full Text Characters (excluding spaces, with Title)** | **2,425** | $\le$ 2,500 characters | **PASS (Under limit by 75 chars)** |
| **Abstract Body Word Count** | 332 words | N/A (Character-based venue) | Informational |
| **Mandatory Headings Included** | `Aim`, `Material & Methods`, `Results`, `Conclusion` | Required | **PASS (100% compliant)** |

---

## 3. Abstract Text (Formatted for Direct Copy-Paste)

### Title
Auditing Clinical Provenance, Molecular Prognosis, and Endpoint Harmonization in Machine Learning Benchmarks for Post-Prostatectomy Recurrence

### Aim
To evaluate endpoint integrity, molecular missingness, and prognostic discrimination in machine learning benchmarks for post-prostatectomy prostate cancer recurrence, establishing an audited clinicopathologic baseline to inform future risk stratification prior to potential salvage evaluation.

### Material & Methods
We conducted a provenance audit of the TCGA-PRAD radical prostatectomy cohort (n=314) linking multi-source clinical XML records with GDC RNA-sequencing profiles. Primary post-prostatectomy recurrence was defined using an audited provenance ledger (n=261; 70 events, 191 double-negatives), excluding indeterminate outcomes, non-primary events, and unconfirmed intervals. Predictor discrimination and calibration were evaluated for standard clinical staging (age, pT, pN, Gleason score) versus clinical plus a candidate 13-gene transcriptomic panel using stratified 5-fold cross-validation in a complete-case paired cohort (n=187; 44 events). Survival was modeled from initial pathologic diagnosis using Cox proportional hazards and Kaplan-Meier estimation (n=248; 57 dated events, 191 censored). Postoperative surveillance PSA (median 511.0 days post-diagnosis) was excluded to reduce temporal leakage; however, exact prostatectomy dates remain unavailable and post-resection pathology is modeled with diagnosis-indexed follow-up.

### Results
In the paired molecular cohort (n=187), the clinical staging model achieved an out-of-fold AUROC of 0.7104 (95% CI 0.6208–0.7913), Brier score of 0.1673, and calibration slope of 0.756. Adding the 13-gene transcriptomic panel yielded an AUROC of 0.7028 (95% CI 0.6046–0.7917), Brier score of 0.1672, and calibration slope of 0.511 (paired ΔAUROC -0.0076 [95% CI -0.0735 to +0.0559]; paired ΔBrier -0.0001 [95% CI -0.0198 to +0.0199]). Audit of recurrence provenance resolved 29 biochemical recurrence status/date discrepancies across registry supplements. In strict continuous survival analysis (n=248), 1-, 3-, and 5-year recurrence-free survival rates were 88.1%, 75.2%, and 65.3% (median not reached), with a cross-validated Harrell's C-index of 0.6716 (apparent C-index 0.7011).

### Conclusion
The evaluated 13-gene panel did not demonstrate an AUROC improvement over standard clinicopathologic staging in this matched dataset. While this retrospective study did not evaluate PSMA imaging, radioligand therapy, or treatment response, transparent endpoint harmonization, mitigation of temporal leakage, and complete-case reporting represent essential methodological prerequisites before recurrence benchmarks can inform future patient selection in prostate theranostics.

---

## 4. Mandatory Submission Disclosures & Policy Statements

> [!WARNING]
> **AUTHOR ACTION REQUIRED BEFORE SUBMISSION:**
> The following items cannot be automatically affirmed by repository evidence and require direct author confirmation in the submission portal.

### Financial & Competing Interests Disclosure
- **Status:** **AUTHOR CONFIRMATION REQUIRED**
- **Action:** The author must explicitly verify and declare in the portal whether any commercial relationships, consultancies, or competing financial interests exist.
- *Default text if no conflicts exist:* `"The author declares that they have no competing financial interests or personal relationships that could have appeared to influence the work reported in this abstract."`

### Funding Disclosure
- **Status:** **AUTHOR CONFIRMATION REQUIRED**
- **Action:** The author must confirm whether any departmental, grant, or institutional funding supported this research.
- *Default text if unfunded:* `"This research was conducted as an investigator-initiated academic study with no external grant funding."`

### Ethical Oversight & Institutional IRB Determination
- **Status:** **AUTHOR CONFIRMATION REQUIRED**
- **Action:** In accordance with institutional research compliance, access to public, de-identified secondary data (NCI GDC) does **not** automatically document a formal Institutional Review Board (IRB) determination. The author must verify their home institution's policy (Texas A&M Human Research Protection Program) regarding whether a formal exempt determination letter is on file or required.
- *Portal description of data source:* `"Data were obtained from public, de-identified clinical and genomic repositories hosted by the NCI Genomic Data Commons (TCGA-PRAD). Formal institutional human subjects determinations must follow home institution policies."`

---

## 5. Explicit Venue-Fit Risk Assessment

- **Official Category:** "Prostate Cancer, including PSMA-Targeted Theranostics" is an official topic on the TWC 2027 abstract portal.
- **Scope Alignment:** The study evaluates risk stratification methodology for post-prostatectomy recurrence. In clinical practice, this is the pivotal juncture where clinicians consider restaging with PSMA-PET/CT and salvage radioligand therapy.
- **Primary Venue Risk:** The work is an oncology informatics and clinical methodology study. It did **not** evaluate PSMA tracer kinetics, radioligand therapeutic administration, radiochemistry, or dosimetry. Reviewers focused strictly on nuclear medicine radiotherapeutics may deem the work out of scope.
- **Classification:** **READY WITH VENUE-FIT RISK**.
