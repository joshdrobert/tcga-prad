# 9th Theranostics World Congress (Osaka TWC 2027) — Official Submission Packet

**Conference:** 9th Theranostics World Congress (TWC 2027)  
**Location:** Osaka, Japan  
**Official Portal:** https://www.theranostics-world-congress.org/abstracts  
**Submission Deadline:** November 1, 2026, 23:59 CET (**TODAY — IMMEDIATE PRIORITY**)  
**Project Repository:** https://github.com/joshdrobert/tcga-prad  
**Live Research Site:** https://joshdrobert.github.io/tcga-prad/  
**Scientific Baseline:** TCGA-PRAD V3.2 Provenance Reconciliation Baseline (`483ab18`)  
**Overall Readiness Classification:** **READY WITH VENUE-FIT RISK**

---

## 1. Submission Metadata & Author Details

- **Title:** Auditing Clinical Provenance, Molecular Prognosis, and Endpoint Harmonization in Machine Learning Benchmarks for Post-Prostatectomy Recurrence
- **Primary Topic / Track:** Prostate Cancer (Alternative / Cross-listing: Artificial Intelligence)
- **Presentation Type Preference:** Poster (Oral presentation accepted if invited)
- **Author:**
  - **Name:** Joshua Robert
  - **Affiliation:** Texas A&M School of Engineering Medicine, Houston, TX, USA
  - **Email:** `[Author primary email]` *(to be filled by human author during submission)*
  - **Corresponding Author:** Joshua Robert
  - **Consortium / Affiliation Note:** No fictitious consortium or unauthorized hospital departments. Strictly affiliated with Texas A&M School of Engineering Medicine.

---

## 2. Character & Word Count Compliance Verification

| Metric | Measured Value | Portal Limit | Status |
| :--- | :--- | :--- | :--- |
| **Abstract Body Characters (excluding spaces)** | **2,253** | $\le$ 2,500 characters | **PASS (Under limit by 247 chars)** |
| **Full Submission Characters (excluding spaces, with Title)** | **2,381** | $\le$ 2,500 characters | **PASS (Under limit by 119 chars)** |
| **Abstract Body Word Count** | 320 words | N/A (Character-based) | Informational |
| **Title Word Count** | 14 words | N/A | Informational |
| **Mandatory Headings Included** | `Aim`, `Material & Methods`, `Results`, `Conclusion` | Required | **PASS (100% compliant)** |

---

## 3. Abstract Text (Formatted for Direct Copy-Paste)

### Title
Auditing Clinical Provenance, Molecular Prognosis, and Endpoint Harmonization in Machine Learning Benchmarks for Post-Prostatectomy Recurrence

### Aim
To evaluate endpoint integrity, molecular missingness, and prognostic discrimination in machine learning benchmarks for post-prostatectomy prostate cancer recurrence, establishing rigorous clinical baseline performance necessary for risk stratification prior to salvage theranostic interventions.

### Material & Methods
We conducted a comprehensive provenance audit of the TCGA-PRAD radical prostatectomy cohort (n=314) linking multi-source clinical XML follow-up records with GDC RNA-sequencing profiles. Primary post-prostatectomy recurrence was defined using an audited multi-source provenance ledger (n=261; 70 events, 191 double-negatives), excluding indeterminate outcomes, non-primary events, and unconfirmed intervals. Predictor discrimination and probability calibration were evaluated for standard clinical staging (age, pT, pN, Gleason score) versus clinical plus a candidate 13-gene transcriptomic panel using stratified 5-fold cross-validation in a complete-case paired cohort (n=187; 44 events). Survival analysis was modeled from initial pathologic diagnosis using Cox proportional hazards and Kaplan-Meier estimation (n=248 with 57 dated events and 191 double-negative censored), strictly excluding surveillance PSA (median 511.0 days post-diagnosis) from baseline predictors to prevent look-ahead bias.

### Results
In the paired molecular cohort (n=187), the clinical staging model achieved an out-of-fold AUROC of 0.7104 (95% CI 0.6208–0.7913), Brier score of 0.1673, and calibration slope of 0.756. Adding the 13-gene transcriptomic panel yielded an AUROC of 0.7028 (95% CI 0.6046–0.7917), Brier score of 0.1672, and calibration slope of 0.511 (paired ΔAUROC -0.0076 [95% CI -0.0735 to +0.0559]; paired ΔBrier -0.0001 [95% CI -0.0198 to +0.0199]). Audit of recurrence provenance resolved 29 biochemical recurrence status/date discrepancies across registry supplements. In strict continuous survival analysis (n=248), 1-, 3-, and 5-year recurrence-free survival rates were 88.1%, 75.2%, and 65.3% (median not reached), with a cross-validated Harrell's C-index of 0.6716 (apparent C-index 0.7011).

### Conclusion
Candidate transcriptomic panels did not provide incremental prognostic discrimination over standard clinicopathologic staging for post-prostatectomy recurrence. Accurate baseline risk stratification for patient selection in emerging theranostic and salvage pathways requires rigorous endpoint harmonization, strict prevention of immortal time bias, and complete-case reporting rather than uncritical reliance on high-dimensional genomic features.

---

## 4. Mandatory Submission Disclosures & Policy Statements

### Financial & Competing Interests Disclosure
> **Declaration:** The author declares that they have no competing financial interests, commercial affiliations, or personal relationships that could have appeared to influence the work reported in this abstract.

### Funding Statement
> **Funding:** This research was conducted as an investigator-initiated academic study with no dedicated commercial or industry grant support.

### Ethical Oversight & Data Provenance Statement
> **Data Access & Compliance:** The results published here are based in part upon data generated by The Cancer Genome Atlas (TCGA) Research Network (TCGA-PRAD). Clinical records and processed molecular RNA sequencing profiles were obtained from open, de-identified public repositories hosted by the National Cancer Institute Genomic Data Commons (GDC) under NIH/NCI data use policies. In accordance with U.S. Federal Policy for the Protection of Human Subjects (Common Rule 45 CFR 46), research utilizing publicly accessible, fully de-identified secondary data is exempt from institutional review board (IRB) review.

### Code & Benchmark Reproducibility Statement
> **Reproducibility:** All parsing pipelines, clinical XML event ledgers, cross-validation scripts, survival models, and evaluation outputs are version-controlled and publicly accessible at [https://github.com/joshdrobert/tcga-prad](https://github.com/joshdrobert/tcga-prad).
