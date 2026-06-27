# Development of an Open-Source Genomic-Clinical Risk Score for Prostate Cancer Using TCGA-PRAD

## Abstract

### Background
TCGA-PRAD provides public clinical and molecular data for reproducible prostate cancer risk stratification.

### Objective
To build an interpretable prostate cancer genomic-clinical risk score using TCGA-PRAD.

### Methods
GDC manifests and open clinical XML supplements were downloaded and parsed into a patient-level clinical cohort. Mortality event counts were assessed before model fitting.

### Results
The current parsed cohort includes 314 patients. Only 4 death events were observed in the parsed clinical supplements, so a mortality model was intentionally not fit.

### Conclusions
The TCGA-PRAD project is ready for endpoint enrichment, especially recurrence/progression endpoints and selected molecular features.

## Introduction

## Methods

### Data Source
TCGA-PRAD open clinical supplement XML files and GDC manifests were downloaded through the GDC API.

### Study Population
Patients with parseable open clinical XML supplements were retained.

### Outcomes
Vital status and overall survival time were parsed where available. Mortality modeling was deferred due to insufficient events.

## Results

Generated outputs include a processed analysis cohort, flattened XML table, missingness table, and clinical characteristics table.

## Discussion

Mortality is not a good first endpoint for the currently parsed PRAD data. Recurrence, progression, grade group, PSA, and molecular predictors should be added before model development.

## Limitations

The current project is clinical-supplement only and does not yet include molecular features or curated recurrence endpoints.

## Conclusion

The project has a reproducible cohort build and a clear next step: endpoint and molecular feature enrichment.

