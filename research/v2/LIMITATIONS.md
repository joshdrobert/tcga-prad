# TCGA-PRAD RESEARCH METHODOLOGY LIMITATIONS & INTEGRITY AUDIT

## 1. Study Limitations

1. **Observational Retrospective Cohort Design**:
   - TCGA-PRAD is an observational multi-center retrospective genomic cohort collected across dozens of tissue source sites without centralized protocolized follow-up.
   - Surveillance schedules (PSA measurement frequency, radiographic imaging triggers) varied by institution.

2. **Asymmetric Clinical Outcome Documentation**:
   - Of 242 recurrence-free patients, 191 have both negative BCR and negative NTE forms.
   - 46 patients have asymmetric documentation (one form negative while the other was omitted), and 5 patients have missing follow-up forms in GDC XML. While cBioPortal classifies these as Disease-Free, they represent unobserved follow-up ambiguity.

3. **Transcriptomic Missingness (27.4%)**:
   - 86 of 314 patients (27.4%) lack paired RNA sequencing in GDC.
   - Imputing these 86 patients degraded calibration slope from 0.777 to 0.547 without improving discrimination.
   - In paired complete cases ($n=228$), adding 13 transcriptomic features yielded an AUROC of 0.686 vs 0.689 for clinical predictors alone.

4. **Undated Recurrence Events (n=6)**:
   - Six patients had confirmed recurrence without exact day tokens in the clinical XML.
   - Rather than fabricating event times via last follow-up fallback, they are excluded from the primary continuous-time Cox model, leaving $n=308$ evaluable cases.

5. **Internal vs External Validation**:
   - The reported Cox PH Concordance Index is **0.703 (apparent)** and **0.688 (5-fold cross-validated out-of-fold)**.
   - It is an internally validated benchmark and must not be described as externally validated.

---

## 2. Invalidated & Corrected Claims

| Historical / Literature Claim | Verified Finding | Status |
| :--- | :--- | :--- |
| "34 patients discordant with cBioPortal" | Naive arithmetic error ($72 - 38 = 34$). True matched discordance is **31 cases** across 269 evaluable cases. | **CORRECTED** |
| "212 patients lack longitudinal follow-up" | Follow-up tokens were flattened (`995\|1373`) and coerced to NaN. Repaired parser extracted follow-up for **311 of 314 patients**. | **INVALIDATED / REPAIRED** |
| "RNA expression panel improves recurrence prediction" | AUROC was 0.689 (clinical) vs 0.686 (molecular) in paired complete cases ($n=228$). No incremental value detected. | **DISPROVEN / EXPLORATORY** |
| "Median imputation improves statistical power" | Imputation induced severe overconfidence, collapsing calibration slope from 0.777 to 0.547. | **INVALIDATED / WARNED** |
| "TCGA-PRAD supports Theranostic World Congress submission" | TCGA-PRAD is surgical genomics with 0% radioligand therapy or PET imaging. Incompatible with Osaka. | **NO-GO MAINTAINED** |
