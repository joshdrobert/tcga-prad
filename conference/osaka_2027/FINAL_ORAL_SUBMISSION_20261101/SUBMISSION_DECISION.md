# SUBMISSION DECISION: TWC 2027 OSAKA (DRAFT A-325)
**Document ID:** DECISION-TWC2027-20261101  
**Timestamp:** 2026-11-01T16:25:00-06:00  
**Target Venue:** 9th Theranostics World Congress (TWC 2027, Osaka, Japan)  
**Deadline:** Today, November 1, 2026, 23:59 CET = 16:59 CST (America/Chicago)  

---

## 1. Selected Candidate
**Selected:** **Candidate A (Recurrence Endpoint Ascertainment Reverses Apparent Genomic Advantage in Prostate Cancer Prognostication)**.

### Why Candidate A is Scientifically Stronger than Potential Alternatives
- **Candidate B (Molecular Target Biology / FOLH1 PSMA Expression):** Rejected because `FOLH1` was not extracted into the current 13-gene expression panel (`BRCA1`, `SPOP`, `FOXA1`, `MYC`, `BRCA2`, `TP53`, `KLK3`, `ATM`, `ERG`, `NKX3-1`, `AR`, `PTEN`, `CHEK2`). Downloading and processing raw TCGA RNA-seq BAM/counts minutes before the deadline is prohibited by safety protocols and cannot be validated in time. Moreover, equating bulk mRNA to in vivo PET uptake would constitute a scientific flaw.
- **Candidate C (Direct AI Contribution to Theranostic Workflow):** Rejected because the `psma-pet-research` workspace contains only an unexecuted study protocol (`PROTOCOL.md`). No volumetric PET/CT DICOM or NIfTI scans exist locally. Any attempt to report imaging metrics would require data fabrication.
- **Candidate A:** The only candidate with 100% accessible, audited, and independently reproducible empirical evidence. Demonstrates a major methodological finding directly relevant to theranostic patient triage.

---

## 2. Independently Verified Scientific Results
1. **Denominators:** $N=314$ radical prostatectomy, $N=261$ definitive adjudication cohort, $N=187$ complete-case molecular cohort with paired RNA.
2. **Endpoint Discordance:** BCR-status-only identified 25 events, whereas longitudinal follow-up documented 19 additional composite recurrence events (43.2% of the 44 harmonized events).
3. **Biological Risk of Discordant Cases ($N=19$):** 17 of 19 (89.5%) had pathologic stage T3/T4; 12 of 19 (63.2%) had Gleason $\ge 8$ (mean Gleason 8.16).
4. **Frozen Common Folds Discrimination:**
   - Naive BCR: Clinical AUROC = 0.6540; Molecular AUROC = 0.7457; $\Delta\text{AUROC} = +0.0917$ (95% CI: +0.0001 to +0.1924).
   - Harmonized Recurrence: Clinical AUROC = 0.7069; Molecular AUROC = 0.6917; $\Delta\text{AUROC} = -0.0153$ (95% CI: -0.0861 to +0.0509).
   - Difference-in-Differences ($\Delta\Delta\text{AUROC}$): $+0.1070$ (95% CI: +0.0111 to +0.2013).
5. **Sensitivity Comparator:** Outcome-dependent stratified folds confirm the identical reversal (+0.1168 naive vs -0.0076 harmonized).

---

## 3. Evidence Supporting Conference Relevance
- In clinical prostate cancer management, risk stratification models and AI algorithms are deployed to triage patients for salvage therapies, restaging PSMA-PET/CT, and theranostic radioligand therapy.
- Incomplete recurrence ascertainment distorts risk predictions by penalizing standard clinical staging, generating false confidence in unvalidated genomic panels.
- Demonstrating this failure mechanism provides an essential research integrity and quality-control safeguard for the TWC audience.

---

## 4. Realistic Subjective Oral-Selection Probability Assessment
- **Can the 70% Aspirational Target be Supported?**
  **No.** Claiming $\ge 70\%$ probability for an observational secondary data analysis without direct radiopharmaceutical or PET imaging at a single-auditorium theranostics congress is scientifically indefensible.
- **Realistic Oral Probability:** **20% to 28%** (Competitive oral candidate for the dedicated AI/Computational session).
- **Realistic Poster Acceptance Probability:** **75% to 85%**.
- **Limiting Factors:** Single-room format, scarcity of computational slots, and absence of direct in vivo imaging endpoints.

---

## 5. Decision on Draft A-325
- **Decision:** **REPLACE Draft A-325 text with the refined Final Oral Submission package (`conference/osaka_2027/FINAL_ORAL_SUBMISSION_20261101/EXACT_PORTAL_COPY.txt`).**
- **Rationale:** The refined abstract incorporates frozen common-fold metrics, paired bootstrap difference-in-differences intervals, explicit characterization of the 19 discordant cases, eliminates em dashes and the word "underscore", clarifies the BCR vs composite recurrence definition, and ensures 100% compliance with TWC character limits.
- **Git Commit Baseline:** `7175a85` (main branch).
- **CI Workflow Status:** Run `41354730306` (PASSED `✓`).
