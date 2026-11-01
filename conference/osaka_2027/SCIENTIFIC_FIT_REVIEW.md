# Scientific & Venue-Fit Review: 9th Theranostics World Congress (Osaka 2027)

**Document Purpose:** Independent scientific and methodological audit of project alignment, scope compatibility, risk factors, and submission defensibility for Osaka TWC 2027.  
**Evaluator Role:** Senior Biomedical Methodologist & Conference Submission Auditor  
**Project:** TCGA-PRAD Post-Prostatectomy Clinical & Genomic Recurrence Benchmark  
**Submission Classification:** **READY WITH VENUE-FIT RISK (AUTHOR REVIEW REQUIRED)**

---

## 1. Executive Summary & Verdict

- **Venue:** 9th Theranostics World Congress (TWC 2027), Osaka, Japan
- **Official Portal:** `https://www.theranostics-world-congress.org/abstracts`
- **Submission Deadline:** November 1, 2026, 23:59 CET (**TODAY — 4:59 PM Central**)
- **Target Category:** "Prostate Cancer, including PSMA-Targeted Theranostics" *(or "Artificial Intelligence")*
- **Recommended Action:** Author may submit today, acknowledging venue-fit risk and providing required personal confirmations (COI, funding, institutional IRB determination).
- **Core Scientific Scope:**
  - The study evaluates post-prostatectomy recurrence risk stratification methodology across 314 TCGA-PRAD patients.
  - Standard clinicopathology (AUROC = 0.7104) is compared to a candidate 13-gene transcriptomic panel (AUROC = 0.7028, paired ΔAUROC = -0.0076, 95% CI: [-0.0735, +0.0559]).
  - The 13-gene panel did not demonstrate an AUROC improvement over clinical staging in this matched dataset, and calibration slope degraded from 0.756 to 0.511.
  - Surveillance PSA (median 511.0 days post-diagnosis) was strictly excluded from baseline Cox regression to reduce temporal leakage; however, exact prostatectomy dates remain unavailable in public records, and post-resection pathology is modeled with diagnosis-indexed follow-up.
  - Reconciled 29 registry biochemical recurrence status/date conflicts.

---

## 2. Explicit Venue-Fit Risk Analysis

> [!WARNING]
> **VENUE-FIT RISK RATING: MODERATE-TO-HIGH**
> Do not artificially strengthen or fabricate a link to theranostics. The conference identity and scope must be presented transparently.

### What the Study Does and Does NOT Contain
| Domain | Study Status | Details |
| :--- | :---: | :--- |
| **Recurrence Risk Methodology** | **PRESENT** | Evaluates prognostic models for post-radical prostatectomy recurrence. |
| **Clinical Endpoint Integrity** | **PRESENT** | Audits 314 patients, reconciles 29 BCR conflicts, excludes surveillance PSA to reduce temporal leakage. |
| **PSMA-PET / SPECT Imaging** | **ABSENT** | Zero patient-level imaging scans, SUV measurements, or segmentation models were evaluated. |
| **Radioligand Therapy (RLT)** | **ABSENT** | Zero radiopharmaceutical administration (e.g., [^177Lu]Lu-PSMA-617, [^225Ac]Ac-PSMA), dosimetry, or treatment response data. |
| **Theranostic Patient Selection** | **UNTESTED** | Study did **not** validate patient selection for salvage theranostic treatments. Its relevance is prospective and methodological only. |

### Why This Research Is Relevant to the Field
In clinical oncology, post-prostatectomy recurrence is the critical transition point where clinicians consider salvage imaging (such as PSMA-PET/CT) and targeted radioligand trials. If prognostic machine learning models deployed at this juncture suffer from temporal leakage (such as surveillance PSA look-ahead) or uncalibrated genomic features, clinical risk stratification will be corrupted. Our study demonstrates that standard clinicopathologic staging remains the primary robust baseline, providing an essential cautionary reference before machine learning pipelines are coupled to theranostic triage.

### Potential Reviewer Feedback & Defensible Position
- **Reviewer Critique:** *"This abstract does not present radiopharmaceutical chemistry, cellular dosimetry, or PSMA imaging."*
- **Defensible Response:** *"The abstract transparently states that PSMA imaging and radioligand therapies were not evaluated. It is submitted under the 'Prostate Cancer' / 'Artificial Intelligence' track to establish the clinicopathologic baseline and methodological safeguards necessary before machine learning models can be used to select patients for salvage theranostic evaluation."*

---

## 3. Compliance and Author Confirmations Required

1. **Conflicts of Interest:** Author must verify their personal financial disclosures directly in the portal.
2. **Funding:** Author must confirm institutional and grant funding status.
3. **IRB Status:** Public, de-identified data do not by themselves document an institutional IRB determination. The author must verify their home institution's human research protection requirements.
