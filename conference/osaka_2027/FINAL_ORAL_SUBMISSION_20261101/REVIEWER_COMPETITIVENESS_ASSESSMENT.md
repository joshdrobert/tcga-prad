# REVIEWER COMPETITIVENESS ASSESSMENT: TWC 2027 OSAKA
**Document ID:** COMP-ASSESS-TWC2027-20261101  
**Timestamp:** 2026-11-01T16:25:00-06:00  
**Target Venue:** 9th Theranostics World Congress (TWC 2027, Osaka, Japan)  
**Evaluated Manuscript:** Candidate A (Recurrence Endpoint Ascertainment Reverses Apparent Genomic Advantage in Prostate Cancer Prognostication)

---

## 1. Multi-Perspective Reviewer Panel Assessment

### Reviewer 1: Theranostics Clinician & Nuclear Medicine Oncologist
- **Perspective:** High-volume clinician treating metastatic castration-resistant prostate cancer with [177Lu]Lu-PSMA-617 and [225Ac]Ac-PSMA radioligand therapy.
- **Evaluation of Candidate A:**
  - *Relevance to Practice:* High conceptual relevance, but indirect. Clinical oncologists and nuclear medicine physicians rely heavily on risk models and molecular panels to decide which postoperative patients require early salvage radiation, systemic therapy, or salvage theranostic PSMA-PET imaging. If commercial risk tools or AI models are artificially inflated because they miss aggressive early recurrences, patients are misclassified.
  - *Major Limitation:* No radiopharmaceutical, PET tracer uptake (e.g., PSMA SUVmax), or theranostic therapeutic response data are included. The data are retrospective observational genomics from TCGA.
  - *Oral Recommendation:* **Borderline / Strong Poster**. In a single-plenary congress, clinicians will prefer direct 177Lu or 225Ac prospective clinical trial data, but this provides a cautionary presentation for the AI/Precision Medicine session.

### Reviewer 2: Nuclear Medicine Imaging & Radiochemistry Researcher
- **Perspective:** Focuses on PET/CT imaging biomarkers, target quantification, radiotracer development, and dosimetry.
- **Evaluation of Candidate A:**
  - *Relevance:* Methodological caution on registry data. While the paper does not evaluate molecular imaging or PET segmentation, its demonstration that incomplete clinical follow-up creates illusory biomarker gains is directly transferable to PET radiomics and AI dosimetry.
  - *Major Limitation:* Does not advance radiopharmaceutical chemistry, dosimetry algorithms, or in vivo imaging physics.
  - *Oral Recommendation:* **Poster candidate**; acceptable for Oral only if the conference has a dedicated session on AI/Computational Oncology Quality Control.

### Reviewer 3: Senior Biomedical AI Methodologist & Biostatistician
- **Perspective:** Focuses on machine learning rigor, study design, validation protocols, and epidemiological biases.
- **Evaluation of Candidate A:**
  - *Methodological Rigor:* **Exceptional (Top 5%)**. The use of frozen common folds across naive and adjudicated outcomes, paired bootstrap difference-in-differences resampling ($B=1,500$), complete exclusion of immortal time surveillance PSA, and thorough characterization of the 19 discordant high-grade cases ($89.5\%$ pT3/T4, $63.2\%$ Gleason $\ge 8$) represents exemplary biostatistical integrity.
  - *Novelty:* Unmasks a profound methodological trap: incomplete endpoint ascertainment penalizes the clinical model (by falsely treating aggressive recurrences as non-events), which artificially manufactures apparent incremental discrimination for molecular biomarkers ($+0.0917$). Harmonizing longitudinal outcomes eliminates this artifact ($\Delta\text{AUROC} = -0.0153$).
  - *Oral Recommendation:* **Strong Oral Candidate** for an AI methodology or computational precision medicine session.

---

## 2. Realistic Subjective Probability Appraisal

### Can a 70% Oral-Selection Probability Be Justified?
**No. An oral-selection probability of 70% or greater is scientifically indefensible for this abstract at TWC 2027.**

### Empirical Limiting Factors:
1. **Congress Format:** TWC traditionally utilizes a "single core program, one room" structure with only 25–40 oral slots total for the entire conference, dominated by prospective theranostic trials (e.g., PSMA-617, FAPi, actinium-225).
2. **Lack of Direct Radiopharmaceutical Data:** TCGA-PRAD is an observational genomic cohort; it contains no PET imaging (PSMA SUV), no SPECT dosimetry, and no radioligand therapy outcomes.
3. **General AI Slot Scarcity:** Specialized theranostics conferences allocate at most 2–4 oral slots to pure computational/AI methodology.

### Defensible, Calibrated Probability Breakdown:
- **Probability of Overall Acceptance (Oral or Poster):** **75% – 85%** (Methodological rigor, clear clinical message, well-structured abstract).
- **Probability of Selection for ORAL Presentation:** **20% – 28%** (Realistic competitive probability for the AI/Prostate Cancer session).
- **Conditional Probability of Oral Selection IF an AI Session is Programmed:** **35% – 45%**.

---

## 3. Recommended Actions to Maximize Oral Appeal
1. Frame the abstract around **Theranostic Patient Selection Integrity**: Emphasize that clinical decision algorithms triage patients for salvage PSMA-PET/CT and radioligand therapy; endpoint distortion invalidates this triage.
2. Emphasize the statistically robust Difference-in-Differences metric ($\Delta\Delta\text{AUROC} = +0.1070$, $p < 0.05$).
3. Maintain 100% data integrity without overclaiming.
