# REPRODUCIBILITY PROTOCOL: TCGA-PRAD (V3.2 AUDITED)

**Project:** Endpoint Harmonization & Molecular Prediction in TCGA-PRAD  
**Git Baseline Commit:** `b2c0f74`  
**Audit Branch:** `audit/2027-conference-final-integrity`  
**Execution Environment:** macOS arm64, Python 3.13, virtualenv  

---

## 1. Quick Reproduction
From the repository root:

```bash
# 1. Activate virtual environment
source .venv/bin/activate

# 2. Run fast smoke test
make smoke

# 3. Run full reproduction pipeline, tests, and claim audit
make reproduce
```

The `make reproduce` command executes:
1. `src/build_recurrence_endpoint.py`: Reconstructs the multi-source clinical XML longitudinal cohort ($N=314$), adjudicates 29 BCR status/date conflicts, separates second primary malignancies ($N=2$), and exports primary definitive ($N=261$) and strict survival ($N=248$) cohorts.
2. `src/cbioportal_endpoint_crosscheck.py`: Audits local definitive recurrence against cBioPortal disease-free status (DFS) and biochemical recurrence (BCR).
3. `src/run_rigorous_validation.py`: Evaluates stratified 5-fold cross-validation on the paired complete-case molecular cohort ($N=187$; 44 events, 143 non-events), fits Kaplan-Meier RFS and Cox proportional hazards models ($N=248$, apparent C = 0.7011, CV C = 0.6716), and exports `research/artifacts/tcga_prad_metrics.json`.
4. `pytest tests/ -v`: Runs all clinical validity, pipeline invariance, conference constraint, and scientific consistency tests.
5. `python scripts/audit_conference_claims.py`: Verifies abstract character bounds and bi-directional consistency against JSON metrics.

---

## 2. Cohort Definitions & Tiers
- **Primary Definitive Recurrence Cohort ($N=261$):** 70 documented prostate recurrence events, 191 double-negative censored (BCR=NO and NTE=NO). Excludes 46 asymmetric negatives and 5 indeterminate cases.
- **Paired Definitive Molecular Cohort ($N=187$):** Complete cases with verified RNA sequencing (44 events, 143 non-events). Pure clinical AUROC = 0.7104, Clinical + 13-gene RNA AUROC = 0.7028 ($\Delta\text{AUROC} = -0.0076$, 95% CI: $[-0.0735, +0.0559]$).
- **Strict Continuous Survival Cohort ($N=248$):** 57 dated recurrence events, 191 double-negative censored. 1-yr RFS = 88.1%, 3-yr RFS = 75.2%, 5-yr RFS = 65.3%.
- **Historical Benchmark Cohort ($N=228$):** Legacy unadjudicated complete-case pairing; retained solely as a historical comparison baseline.

---

## 3. Key Verified Artifacts
- Clinical cohort: `data_processed/tcga_prad_recurrence_cohort.csv`
- Model performance: `research/v3.2/PAIRED_MOLECULAR_COMPARISON.csv`
- 4-Tier survival comparison: `research/v3.2/SURVIVAL_SENSITIVITY_ANALYSIS.csv`
- Discrepancy ledgers: `research/v3.2/BCR_STATUS_DATE_CONFLICTS.csv` and `EVENT_ADJUDICATION_CHANGES.csv`
- Predictor timing audit: `research/v3.2/PREDICTOR_TIMING_AUDIT.csv`
- Master metrics JSON: `research/artifacts/tcga_prad_metrics.json`
- Conference submission packet: `conference/osaka_2027/PORTAL_READY_2026-11-01.txt`
