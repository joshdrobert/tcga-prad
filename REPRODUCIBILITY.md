# REPRODUCIBILITY PROTOCOL: TCGA-PRAD

**Project:** Endpoint Harmonization & Molecular Prediction in TCGA-PRAD  
**Git Branch:** `research/2027-osaka-validation`  
**Base Commit:** `fe0d23a`  
**Execution Environment:** macOS arm64, Python 3.13, uv virtualenv  

---

## 1. Quick Reproduction
From the repository root:

```bash
# 1. Create and activate virtual environment
uv venv .venv --python 3.13
source .venv/bin/activate
uv pip install -r _shared/00_environment/requirements.txt pytest

# 2. Run full reproduction pipeline and test suite
make reproduce
```

The `make reproduce` command executes:
1. `src/build_recurrence_endpoint.py`: Extracts longitudinal XML encounters, cleans follow-up days and composite recurrence events.
2. `src/cbioportal_endpoint_crosscheck.py`: Audits local composite against cBioPortal biochemical recurrence (BCR).
3. `src/run_rigorous_validation.py`: Evaluates 5-fold cross-validation on full (n=314) and paired complete cases (n=228), fits Kaplan-Meier RFS and Cox proportional hazards models, and generates figures and `research/artifacts/tcga_prad_metrics.json`.
4. `pytest tests/test_pipeline.py -v`: Runs all data integrity and invariance unit tests.

## 2. Key Verified Artifacts
- Clinical cohort: `data_processed/tcga_prad_recurrence_cohort.csv`
- Model performance: `tables/table_prad_rigorous_validation.csv`
- Complete-case validation: `tables/table_prad_paired_complete_case_validation.csv`
- Survival analysis: `tables/table_prad_survival_analysis.csv`
- cBioPortal discordance: `tables/table_cbioportal_discordance_breakdown.csv`
- Calibration figure: `figures/prad_calibration.png`
- Kaplan-Meier figure: `figures/prad_km_survival.png`
- Structured JSON: `research/artifacts/tcga_prad_metrics.json`
