# TCGA-PRAD V2 CLINICAL ENDPOINT DEFINITIONS & HARMONIZATION SPECIFICATION

## 1. Context & Motivation

Clinical recurrence post-radical prostatectomy exhibits multiple disease phenotypes that must not be conflated:
1. **Biochemical Recurrence (BCR):** Detectable and rising serum PSA ($\ge 0.2\text{ ng/mL}$) post-surgery.
2. **New Tumor Event (NTE):** Histologically or radiographically confirmed local recurrence in the prostatic bed, regional lymph node metastasis, or distant metastasis (e.g. bone, lung).
3. **Progression on Adjuvant / Neoadjuvant Hormone Therapy (HT):** Clinical progression while receiving androgen deprivation therapy.

In TCGA-PRAD, these dimensions are recorded across disparate XML nodes and follow-up encounter forms.

---

## 2. Formal Endpoint Definitions

### Component 1: Biochemical Recurrence (BCR)
- **Positive:** GDC XML tag `<biochemical_recurrence>` contains `"YES"`, `"TRUE"`, or `"1"`.
- **Negative:** GDC XML tag `<biochemical_recurrence>` contains `"NO"`, `"FALSE"`, or `"0"`.
- **Missing:** Node absent or value empty/whitespace.
- **Event Timing:** Extracted from `<days_to_first_biochemical_recurrence>`. If multiple tokens exist, the minimum positive value is taken as the primary event onset.

### Component 2: New Tumor Event (NTE)
- **Positive:** `<new_tumor_event_after_initial_treatment>` contains `"YES"`, `"TRUE"`, `"1"`, or multi-encounter string containing `"YES"` (e.g. `"NO|YES"`).
- **Negative:** Value strictly `"NO"`, `"FALSE"`, or `"0"`.
- **Missing:** Node absent or unrecorded.
- **Event Timing:** Extracted from `<days_to_new_tumor_event_after_initial_treatment>`. Minimum positive value is selected.

### Component 3: Post-HT Tumor Progression
- **Positive:** `<tumor_progression_post_ht>` contains `"YES"`.
- **Negative:** Value contains `"NO"`.
- **Not Applicable / Unrecorded:** Value is NaN (patient did not receive hormone therapy or progression status not documented).

---

## 3. Composite Recurrence Endpoint

For binary risk modeling:
$$\text{Composite Recurrence} = \begin{cases} 1 & \text{if BCR is Positive OR NTE is Positive OR HT Progression is Positive} \\ 0 & \text{otherwise} \end{cases}$$

### Continuous Time-to-Event Survival Endpoint
For longitudinal recurrence-free survival (RFS):
- **Event Case ($T_i, \Delta_i = 1$):**
  $$T_i = \min(\{t \in \{\text{days\_to\_bcr}, \text{days\_to\_nte}\} : t > 0\})$$
  *Condition:* Valid only if $\min(\cdot)$ is non-empty. If a patient is composite-positive but neither date exists, the patient is flagged as an **Undated Event** and excluded from exact continuous-time proportional hazards modeling.
- **Censored Case ($T_i, \Delta_i = 0$):**
  $$T_i = \max(\{\text{days\_to\_last\_followup}, \text{days\_to\_death}\})$$
  *Condition:* Must be $> 0$.
