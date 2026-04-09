# Transdiagnostic Predictors of Treatment Outcomes

This repository contains generic scripts used to identify, model, and validate transdiagnostic clinical and genetic predictors of treatment outcomes (remission and treatment resistance) across major psychiatric disorders (SCZ, BD, MDD).

The code is designed to be cohort-agnostic and reusable across datasets with similar structure.

---

## Repository structure

The main folder contains scripts corresponding to the full analytical workflow:

```
scripts_transdiagnostic_predictors_treatment_outcomes/
├── 00-preprocessing.R
├── 01-univariable_analyses.R
├── 02-multivariable_analyses.R
├── 03-univariable_age_sex_adjusted_for_meta-analyses.R
├── 04-meta-analyses_across_cohorts.R
├── 05-external_validation.R
```

Each script represents one step of the pipeline and can be adapted to different cohorts by modifying only input data and variable mappings.

---

## Analytical workflow

The analysis is structured in six sequential steps:

1. Preprocessing  
2. Univariable screening  
3. Multivariable modelling (within cohort)  
4. Age- and sex-adjusted models for meta-analysis  
5. Meta-analysis across cohorts  
6. External validation  

---

## Scripts description

### 00 – Preprocessing

Builds a cohort-level analytic dataset from raw data.

Main steps:
- Merge multiple raw sources  
- Harmonise sociodemographic and clinical variables  
- Derive grouped predictors (e.g. ordinal categories)  
- Merge polygenic scores and ancestry PCs  
- Residualise PRSs on PCs  
- Define treatment outcomes (Remission, TR1, TR2)  

Note:  
This script is intentionally flexible and should be adapted to each cohort’s raw data structure.

---

### 01 – Univariable analyses

Performs univariable logistic regression analyses for each predictor.

Main features:
- One model per predictor  
- Descriptive summaries by outcome status  
- Odds ratios (OR), 95% CI, and p-values  
- Consistent output tables across outcomes and cohorts  

Typical use:
- Initial screening of candidate predictors (e.g. p < 0.10)

---

### 02 – Multivariable analyses

Implements penalised multivariable models within each cohort.

Main features:
- Elastic Net logistic regression  
- Nested cross-validation (outer AUC, inner lambda selection)  
- Alpha tuning  
- Gelman scaling (2 SD)  
- Model selection based on performance and sparsity  

Output:
- Selected predictors  
- Out-of-fold AUC estimates  
- Final fitted model objects  

---

### 03 – Age- and sex-adjusted univariable models (for meta-analysis)

Generates harmonised effect estimates across cohorts.

Model:
```
outcome ~ predictor + AGE + SEX
```

Main features:
- Cohort-specific adjusted models  
- Standardised log(OR) and SE extraction  
- Output formatted for meta-analysis  

Note:  
The same script can be used for:
- Remission  
- TR1 (treatment-based resistance)  
- TR2 (symptom-based resistance)  

---

### 04 – Meta-analyses across cohorts

Combines cohort-level estimates using random-effects models.

Main features:
- Explicit cross-cohort harmonisation of predictors  
- REML random-effects meta-analysis  
- Heterogeneity metrics (I², τ², Q)  
- Leave-one-out analyses  
- Optional forest plots  

Key concept:  
Predictor correspondence across cohorts is defined explicitly rather than assumed.

---

### 05 – External validation

Validates selected predictors in an independent cohort.

Workflow:

Step 1 – Univariable validation  
- Same logic as 01-univariable_analyses  
- Restricted to predictors selected from meta-analyses:  
  - pooled p < 0.10  
  - consistent direction in ≥2/3 cohorts  

Step 2 – Multivariable validation  
- Standard logistic regression (glm)  
- No penalisation or tuning  
- Joint estimation of selected predictors  

---

## Outcomes

The analyses focus on three harmonised outcomes:

- Remission  
- Treatment resistance (TR1): treatment-based phenotype  
- Treatment resistance (TR2): symptom-based phenotype  

---

## General principles

- Scripts are generic templates, not cohort-specific pipelines  
- Cohort differences are handled at the level of:
  - input data  
  - variable naming  
  - outcome definitions  
- Preprocessing is intentionally flexible to accommodate heterogeneity across studies  
- Downstream steps enforce consistency in modelling and reporting  

---

## Requirements

R (≥ 4.0 recommended)

Key packages:
- dplyr  
- tidyverse  
- glmnet  
- pROC  
- metafor  
- broom  

---

## Notes

- Scripts are intended to be run sequentially  
- File paths must be adapted locally  
- Raw data are not included in this repository  
- Outputs are generated as tables and model objects for downstream use  

---

## Contact

For questions or collaboration, please contact the repository owner.
