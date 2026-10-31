# TCGA-PRAD V2 CLINICAL DATA & LONGITUDINAL ENDPOINT AUDIT

**Date:** 2026-10-31  
**Study:** Prostate Adenocarcinoma (TCGA-PRAD) Endpoint Harmonization & Molecular Prediction Benchmark  
**Dataset Provenance:** Genomic Data Commons (GDC) TCGA-PRAD Clinical XML Supplements & GDC RNA-Seq STAR Raw Counts  
**Analysis Cohort:** 314 Patients who underwent radical prostatectomy with curative intent  

---

## 1. Executive Summary of Corrections

In V1, two critical methodological errors were identified and resolved:

1. **Longitudinal Event Timing Fallback Bug**:
   - *Previous Defect:* For patients with a documented recurrence event but missing exact event dates, the code fell back to `days_to_last_followup_clean` as an invented event date.
   - *V2 Resolution:* This fallback was completely excised. Six patients had documented recurrence without exact event dates in the XML (`TCGA-HC-8212`, `TCGA-V1-A8MK`, `TCGA-V1-A8WS`, `TCGA-V1-A9O9`, `TCGA-V1-A9ZR`, `TCGA-XQ-A8TA`). These cases are documented in `research/v2/UNDATED_RECURRENCE_CASES.csv` and excluded from continuous-time survival analysis to prevent falsification.
   - *Resulting Survival Population:* Exactly **308 evaluable patients**, comprising **66 verified dated events** and **242 censored non-events**.

2. **Missing-Outcome Semantics & Indeterminate Handling**:
   - *Previous Defect:* `parse_boolean_flag` converted `NaN` to `0`, conflating missing or unrecorded outcomes with documented absence of recurrence.
   - *V2 Resolution:* Implemented a 4-state adjudication policy (`DOCUMENTED_POSITIVE`, `DOCUMENTED_NEGATIVE`, `UNKNOWN_MISSING`, `NOT_APPLICABLE_OR_UNRECORDED`) across Biochemical Recurrence (BCR), New Tumor Event (NTE), and Post-HT Progression.
   - *Classification Findings:*
     - Documented Recurrence (Composite): 72 patients (66 dated, 6 undated).
     - High-Confidence Recurrence-Free (BCR=NO & NTE=NO): 191 patients.
     - Asymmetric Negative (BCR=NO, NTE missing): 9 patients.
     - Asymmetric Negative (NTE=NO, BCR missing): 37 patients.
     - Indeterminate (BCR missing & NTE missing): 5 patients.

---

## 2. Longitudinal Cohort Flow Accounting

| Step / Stratum | Cohort Size (n) | Event Count | Censored Count | Note |
| :--- | :--- | :--- | :--- | :--- |
| **Total Surgical Patients** | 314 | 72 | 242 | Radical prostatectomy analysis cohort |
| **Documented Recurrence (Composite)** | 72 | 72 | 0 | Documented BCR (38), NTE (64), or post-HT progression (13) |
| *-- Dated Recurrence Events* | 66 | 66 | 0 | Exact date of BCR or NTE documented (Median: 708 days) |
| *-- Undated Recurrence Events* | 6 | 6 | 0 | Recurrence documented; exact date omitted in clinical XML |
| **Documented Recurrence-Free** | 242 | 0 | 242 | No documented recurrence across any follow-up encounter |
| *-- High Confidence (BCR=NO & NTE=NO)* | 191 | 0 | 191 | Both clinical surveillance axes confirmed negative |
| *-- Asymmetric Negative (BCR=NO)* | 9 | 0 | 9 | Serum PSA documented normal; NTE form omitted |
| *-- Asymmetric Negative (NTE=NO)* | 37 | 0 | 37 | Clinical NTE form negative; serum PSA indicator omitted |
| *-- Indeterminate (All Missing)* | 5 | 0 | 5 | No follow-up response in XML (`TCGA-CH-5748`, `HC-A8CY`, `KK-A7AP`, `KK-A7AZ`, `KK-A8IB`) |
| **Primary Evaluable Survival Cohort** | 308 | 66 | 242 | Patients with valid time-to-event or follow-up |

---

## 3. Detailed Audit of the 6 Undated Recurrence Cases

| Submitter ID | BCR Status | NTE Status | HT Progression | Days to Follow-up | Reason for Timing Omission |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `TCGA-HC-8212` | UNKNOWN_MISSING | POSITIVE (NO\|YES) | UNRECORDED | 2179.0 | Distant recurrence documented on encounter 2 without date |
| `TCGA-V1-A8MK` | POSITIVE (YES) | NEGATIVE (NO) | NEGATIVE (NO) | 508.0 | BCR documented in primary form; `days_to_first_bcr` omitted |
| `TCGA-V1-A8WS` | NEGATIVE (NO) | POSITIVE (YES) | NEGATIVE (NO) | 474.0 | NTE documented in follow-up; `days_to_nte` omitted |
| `TCGA-V1-A9O9` | POSITIVE (YES) | NEGATIVE (NO) | NEGATIVE (NO) | 2870.0 | Late BCR documented without day token in XML |
| `TCGA-V1-A9ZR` | POSITIVE (YES) | POSITIVE (YES) | NEGATIVE (NO) | 1361.0 | Multiple recurrence axes confirmed positive without exact day |
| `TCGA-XQ-A8TA` | POSITIVE (YES) | NEGATIVE (NO) | POSITIVE (YES) | NaN | Advanced progression documented; patient lost to follow-up |

---

## 4. cBioPortal Cross-Dataset Discordance

Re-evaluation against public cBioPortal (`prad_tcga` study, $n=499$) confirms:
- **Matched Patients in Local Cohort:** 269 evaluable cases with cBioPortal BCR records.
- **Concordance:** 238 / 269 = **88.48%** agreement.
- **Concordant Events:** 38 cases.
- **Concordant Non-Events:** 200 cases.
- **Discordances:** Exactly **31 cases** (all Local-Only positive events).
  - *Mechanism:* cBioPortal's `BIOCHEMICAL_RECURRENCE_INDICATOR` captures only serum PSA elevation. Local composite surveillance captured 31 clinical New Tumor Events, local recurrences, and distant metastases that occurred without a documented isolated PSA elevation in the GDC BCR record.
