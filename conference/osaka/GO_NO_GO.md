# OSAKA 2027 CONFERENCE DECISION GATE: FORMAL NO-GO REPORT

**Congress:** 9th Theranostics World Congress (TWC 2027)  
**Location:** Osaka, Japan  
**Dates:** March 11–14, 2027  
**Official Submission Deadline:** November 1, 2026, 23:59 CET (~48 hours from audit)  
**Evaluation Date:** October 30, 2026  
**Evaluator:** Lead Biomedical Research Engineer & Coauthor  
**Decision:** **NO-GO FOR SUBMISSION TO 9TH THERANOSTICS WORLD CONGRESS**

---

## 1. Executive Summary & Rationale

Submitting to the 9th Theranostics World Congress requires high-quality, authentic scientific research directly addressing nuclear medicine theranostics, radioligand therapy (RLT, such as 177Lu-PSMA-617 or 225Ac-PSMA), PSMA-PET/CT molecular imaging, dosimetric modeling, or radiopharmaceutical response assessment.

Our audit comprehensively evaluated two candidate options:
1. **Existing Candidate: `tcga-prad` (Surgical & Genomic Cohort)**
2. **Alternative Exploration: TCIA `psma-pet-ct-lesions` Public Dataset**

Both options were independently audited against the core mission principles:
- *Real science before impressive presentation.*
- *Verified data before model performance claims.*
- *Reproducible execution before manuscript conclusions.*
- *Scientific accuracy before conference deadlines.*

**Conclusion:** Neither study can be ethically, methodologically, or physically submitted as a completed research study for the November 1, 2026 deadline. Fabricating theranostic claims for a surgical genomic cohort or rushing unverified deep-learning imaging metrics within 48 hours on 22 GB free disk space would violate fundamental scientific integrity standards.

---

## 2. In-Depth Audit of Option 1: `tcga-prad`

### 2.1 Study Domain vs Congress Scope Mismatch
- **Nature of TCGA-PRAD:** The Cancer Genome Atlas Prostate Adenocarcinoma (TCGA-PRAD) collection is a retrospective cohort of primary prostate cancer tissue specimens derived from radical prostatectomies. It contains surgical pathology (Gleason score, pathologic TNM staging), bulk RNA sequencing (STAR counts), and clinical follow-up for post-surgical biochemical recurrence (BCR).
- **Theranostic Content:** TCGA-PRAD contains:
  - 0 radioligand therapies (no Lutetium-177, Actinium-225, or other therapeutic radionuclides).
  - 0 PSMA-PET/CT scans or molecular imaging modalities.
  - 0 dosimetry measurements or target engagement kinetics.
- **Scientific Misrepresentation Risk:** Presenting an endpoint-harmonization study of post-surgical biochemical recurrence as a "theranostic artificial intelligence" study is fundamentally false. Referees and scientific committees at the Theranostics World Congress would immediately recognize that no theranostic agent or radiotracer is investigated.

### 2.2 Unresolved Methodological Limitations in TCGA-PRAD
Even within its legitimate domain (genomic/clinical recurrence risk), the current codebase has outstanding issues:
- 86 of 314 patients lack RNA sequencing data and were silently median-imputed in the existing validation script.
- Naive string matching (`"YES|BIOCHEMICAL|..."`) flattened complex longitudinal XML encounters.
- Recurrence is modeled as a binary cross-sectional outcome without accounting for right-censoring in follow-up.
- Manuscript claims 34 discordances against cBioPortal, whereas the actual matched data exhibits 31 discordances out of 269 evaluable cases.

---

## 3. In-Depth Audit of Option 2: TCIA PSMA-PET/CT Feasibility

### 3.1 Dataset Specification
- **Collection:** `psma-pet-ct-lesions` (The Cancer Imaging Archive)
- **Cohort:** 597 PSMA-PET/CT studies from 378 patients with biopsy-proven prostate cancer, manual lesion segmentations, and clinical metadata.
- **Data Size:** Approximately **117 GB** uncompressed DICOM and NIfTI volumes.

### 3.2 Resource & Hardware Infeasibility
- **Local Disk Space:** Total capacity: 460 GiB; Used: 411 GiB; Available: **22 GiB**.
- **Storage Deficit:** Downloading the full 117 GB dataset would immediately exhaust available disk space, resulting in disk write failures and system instability.
- **Timeline:** The November 1, 2026 deadline is in approximately 48 hours.
- **Computational Requirements:** Preprocessing 3D multimodal volumes (PET SUV normalization, CT Hounsfield unit clipping, resampling to isotropic resolution, rigid registration), partitioning into patient-disjoint folds, training a 3D U-Net or nnU-Net baseline, and conducting rigorous failure detection / out-of-distribution robustness analysis requires hundreds of GPU hours on multi-GPU nodes. On an Apple Silicon workstation (16 GB unified memory), running full volumetric 3D training and validation across 597 studies cannot be completed within 48 hours.
- **Scientific Novelty Requirement:** Publishing simply the existing baseline without novel, validated failure detection would be an unnovel duplicate of published literature.

---

## 4. Osaka Decision Matrix

| Dimension | Option 1: `tcga-prad` | Option 2: TCIA `psma-pet-ct-lesions` | Threshold for GO |
| :--- | :--- | :--- | :--- |
| **Theranostic Relevance** | **FAILED** (0% theranostics) | **PASSED** (100% PSMA-PET/CT) | Mandatory theranostic content |
| **Data Authenticity** | **PASSED** (GDC authentic data) | **PASSED** (TCIA authentic data) | Must be verified authentic data |
| **Hardware & Disk Capacity** | **PASSED** (<1 GB footprint) | **FAILED** (117 GB required vs 22 GB free) | Storage & compute must be feasible |
| **Completion by Nov 1, 2026** | **PASSED** (Completed cohort) | **FAILED** (Cannot train/validate 3D models in 48h) | Real completed results required |
| **Methodological Rigor** | **PROVISIONAL** (Needs fix for 86 missing RNA) | **NOT EXECUTABLE** (Cannot complete) | Zero fabricated metrics |
| **FINAL GATE STATUS** | **NO-GO** | **NO-GO** | **OVERALL DECISION: NO-GO** |

---

## 5. Forward-Looking Research Roadmap for Osaka & Future Conferences

1. **Preserve Scientific Honesty:**
   By declaring NO-GO for Osaka, we protect the research team's scientific reputation from submitting a fabricated abstract or misaligned paper.
2. **Reposition TCGA-PRAD Appropriately:**
   Continue refining `tcga-prad` for a clinical informatics or urologic oncology venue (e.g., *Journal of Urology*, *European Urology Oncology*, or *AMIA Clinical Informatics*), where endpoint harmonization, missingness auditing, and open-source genomics are prized contributions.
3. **Establish a Dedicated PSMA-PET/CT Research Track:**
   Create `~/Research/March2027/psma-pet-research/` as an isolated sibling workspace. Structure a proper protocol for TCIA PSMA-PET data subset ingestion, explore cloud/cluster compute options for the 117 GB imaging pipeline, and target the 2027 Society of Nuclear Medicine and Molecular Imaging (SNMMI) or European Association of Nuclear Medicine (EANM) annual congresses.
