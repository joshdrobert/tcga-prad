# TCGA-PRAD Baseline Inventory & System State

**Date:** 2026-11-01  
**Repository:** `https://github.com/joshdrobert/tcga-prad`  
**Git Baseline Commit:** `b2c0f74`  
**Audit Branch:** `audit/2027-conference-final-integrity`  
**Target Conference:** 9th Theranostics World Congress (Osaka 2027)  

---

## 1. Operating Environment & Hardware

- **Operating System:** macOS (Darwin 27.0.1, arm64)
- **CPU:** Apple Silicon (8 CPU cores)
- **Python Version:** 3.13.1 (venv: `/Users/josh/Research/March2027/tcga-prad/.venv/bin/python`)
- **Key Installed Libraries:**
  - `numpy`: 2.5.4
  - `scipy`: 1.18.1
  - `pandas`: 3.1.0
  - `scikit-learn`: 1.9.1
  - `lifelines`: 0.30.0
  - `matplotlib`: 3.11.2
  - `pytest`: 9.2.0

---

## 2. Authoritative Data & Artifact Inventory

| File Path | Description | Rows / Size | SHA256 (Prefix) |
| :--- | :--- | :--- | :--- |
| `data_processed/tcga_prad_recurrence_cohort.csv` | Primary definitive & exploratory recurrence cohort | 314 patients | `d37c8a49...` |
| `data_processed/tcga_prad_clinical_expression_cohort.csv` | Paired clinical and RNA-seq profile dataset | 228 patients | `a5f21bc9...` |
| `research/artifacts/tcga_prad_metrics.json` | Master RFC 8259 JSON metrics manifest | 11.2 KB | `829aef10...` |
| `research/v3.2/BCR_STATUS_DATE_CONFLICTS.csv` | 29 reconciled registry status/date discrepancies | 29 cases | `41f893cd...` |
| `research/v3.2/EVENT_ADJUDICATION_CHANGES.csv` | Corrected provenance ledger | 27 cases | `5b719488...` |
| `research/v3.2/SURVIVAL_SENSITIVITY_ANALYSIS.csv` | 4-tier survival analysis comparison table | 4 rows | `f28a5061...` |
| `research/v3.2/PAIRED_MOLECULAR_COMPARISON.csv` | Primary paired clinical vs 13-gene comparison | 4 models | `6f38aa71...` |
| `research/v3.2/PREDICTOR_TIMING_AUDIT.csv` | Denominator and timing landmark audit | 8 predictors | `429ba260...` |
| `conference/osaka_2027/READY_TO_PASTE_ABSTRACT.txt` | Submission text for Osaka TWC 2027 | 332 words | `9628fc04...` |
| `conference/osaka_2027/PORTAL_READY_2026-11-01.txt` | Pre-formatted handoff for draft A-325 | 332 words | `1b89af52...` |
| `conference/osaka_2027/PORTAL_ENTRY_VALUES.md` | Field-by-field entry values for draft A-325 | Reference guide | `3e7a91fc...` |

---

## 3. Cohort Architecture & Reproduction Gates

- **Total Radical Prostatectomy Patients:** $N = 314$
- **Primary Definitive Recurrence Cohort:** $N = 261$ (70 positive recurrence events, 191 double-negatives)
- **Paired Complete-Case Molecular Cohort:** $N = 187$ (44 events, 143 double-negatives)
- **Continuous-Time Survival Cohort (Strict):** $N = 248$ (57 dated events, 191 double-negatives censored; median RFS not reached, 1-yr RFS 88.1%, 3-yr RFS 75.2%, 5-yr RFS 65.3%)
- **Exploratory Sensitivity Survival Cohort:** $N = 294$ (57 dated events, 237 censored)
- **Permissive Sensitivity Survival Cohort:** $N = 255$ (64 dated events, 191 censored)
