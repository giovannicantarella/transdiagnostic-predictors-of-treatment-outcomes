
# Transdiagnostic Predictors of Treatment Outcomes

This repository contains generic scripts used to identify and analyse transdiagnostic clinical and genetic factors associated with treatment outcomes (remission and treatment resistance) across major psychiatric disorders (SCZ, BD, MDD).

The analytical workflow combines cohort-specific association analyses, penalised multivariable prediction models, random-effects meta-analyses across diagnostic groups, and independent replication of selected associations.

The code is designed to be applied across different cohorts and adapted to datasets with similar structure.

---

## Repository structure

The main folder contains scripts corresponding to the full analytical workflow:

```text
scripts_transdiagnostic_predictors_treatment_outcomes/
├── 00-preprocessing.R
├── 01-univariable_analyses.R
├── 02-multivariable_analyses.R
├── 03-univariable_age_sex_adjusted_for_meta-analyses.R
├── 04-meta-analyses_across_cohorts.R
├── 05-independent_replication.R
```

Each script represents one step of the pipeline and can be adapted to different cohorts by modifying input data, outcome definitions, and variable mappings.

---

## Analytical workflow

The analysis is structured in six sequential steps:

1. Preprocessing
2. Univariable screening
3. Multivariable modelling and calibration (within cohort)
4. Age- and sex-adjusted models for meta-analysis
5. Meta-analysis across cohorts
6. Independent replication (univariable and multivariable analyses)

---

## Scripts description

### 00 – Preprocessing

Builds a cohort-level analytic dataset from raw data.

Main steps:
- Merge multiple raw sources
- Harmonise sociodemographic and clinical variables
- Derive grouped predictors (e.g. ordinal categories)
- Merge polygenic scores and ancestry PCs
- Residualise PGSs on PCs
- Define treatment outcomes (Remission, TR1, TR2, where applicable)

Note:
This script is intentionally flexible and should be adapted to each cohort's raw data structure.

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

Note:
Univariable screening was performed before nested cross-validation, rather than separately within each outer training fold. This may introduce some optimism into subsequent out-of-fold performance estimates.

---

### 02 – Multivariable analyses

Implements penalised multivariable models within each discovery cohort.

Main features:
- Elastic-net logistic regression
- Nested cross-validation (outer AUC, inner lambda selection)
- Alpha tuning
- Gelman scaling (2 SD), with scaling parameters estimated on training data
- Model selection based on performance and sparsity
- Out-of-fold predictions for discrimination and calibration assessment

Calibration assessment:
- Uses out-of-fold predicted probabilities
- Divides participants into five groups according to predicted probability
- Calculates mean predicted probability and observed event rate within each group
- Generates calibration plots against the line of perfect calibration

Output:
- Selected predictors
- Out-of-fold AUC estimates and 95% CI
- Final fitted model objects
- Calibration plots (PDF)

Note:
Predictors are selected upstream of nested cross-validation. Alpha selection is also based on outer-CV performance. Consequently, the resulting performance and calibration estimates should not be interpreted as fully independent of all model-selection steps.

---

### 03 – Age- and sex-adjusted univariable models (for meta-analysis)

Generates harmonised effect estimates across cohorts.

General model:

```text
outcome ~ predictor + AGE + SEX
```

When age is the predictor of interest, models are adjusted for sex only. Conversely, when sex is the predictor of interest, models are adjusted for age only.

Main features:
- Cohort-specific adjusted models
- Extraction of log(OR) estimates and standard errors
- Output formatted for meta-analysis

The same modelling approach can be applied to:
- Remission
- TR1 (treatment-based resistance)
- TR2 (symptom-based resistance)

Outcome availability and definitions depend on the cohort.

---

### 04 – Meta-analyses across cohorts

Combines cohort-level estimates using random-effects models.

Main features:
- Explicit cross-cohort harmonisation of predictors
- REML random-effects meta-analysis
- Heterogeneity metrics (I², τ², Q)
- Leave-one-out analyses when three cohorts contribute
- Optional forest plots

Key concept:
Predictor correspondence across cohorts is defined explicitly rather than assumed.

Pooled associations are interpreted alongside the direction of cohort-specific effects and between-cohort heterogeneity. A statistically significant pooled association does not necessarily indicate a stable effect magnitude across diagnoses.

---

### 05 – Independent replication

Tests selected associations in an independent MDD cohort (STAR*D).

The script implements univariable and multivariable replication analyses, using predictors identified from the discovery meta-analyses.

Workflow:

Step 1 – Univariable independent replication
- Descriptive summaries of candidate predictors by outcome status
- Univariable logistic regression for each selected predictor
- Odds ratios (OR), 95% confidence intervals, and p-values
- Overall predictor p-values based on likelihood-ratio tests
- Predictors selected from discovery meta-analyses according to:
  - pooled p < 0.10
  - consistent direction in at least two of the three contributing discovery cohorts
- No additional significance-based screening in the replication cohort

Step 2 – Multivariable independent replication
- Standard multivariable logistic regression (`glm`)
- No penalisation or hyperparameter tuning
- Joint estimation of predefined predictor sets
- Separate complete-case analysis for each model
- Adjusted ORs, 95% confidence intervals, and p-values
- In-sample AUC estimates and 95% confidence intervals

For treatment resistance, predictor sets derived from the TR1 and TR2 discovery meta-analyses are evaluated in separate models against the same treatment-resistance outcome in the replication cohort.

The candidate predictors tested in the univariable analyses and the predictor sets entered into the multivariable models are specified separately. Multivariable predictor inclusion is not automatically determined by univariable significance in the replication cohort.

Output:
- Descriptive statistics by outcome status
- Univariable ORs, 95% confidence intervals, and p-values
- Multivariable adjusted ORs, 95% confidence intervals, and p-values
- In-sample AUC estimates and 95% confidence intervals
- CSV tables containing descriptive, univariable, and multivariable results
- Summary tables of multivariable model performance

Note:
Independent replication assesses whether selected associations are observed in a separate cohort and estimates their joint effects through model refitting. The reported AUCs reflect apparent discrimination within the replication sample, rather than predictive performance obtained by directly applying the original discovery models.

---

## Outcomes

The analyses focus on remission and treatment resistance, harmonised across cohorts where possible.

- Remission: defined using cohort-specific clinical symptom thresholds
- TR1: treatment-based resistance phenotype
- TR2: symptom-based resistance phenotype

TR1 and TR2 were operationalised separately in the SCZ and BD discovery cohorts. In the MDD discovery cohort (GSRD) and independent replication cohort (STAR*D), a single treatment-resistance outcome was defined.

Outcome definitions, assessment periods, and treatment contexts differ across cohorts and should be considered when interpreting cross-diagnostic associations.

---

## General principles

- Scripts are generic templates, not cohort-specific pipelines
- Cohort differences are handled at the level of:
  - input data
  - variable naming
  - outcome definitions
  - predictor sets
- Preprocessing is intentionally flexible to accommodate heterogeneity across studies
- Downstream steps support consistency in modelling and reporting
- Associations across diagnostic groups are assessed using random-effects meta-analysis
- Predictive discrimination is evaluated separately from cross-cohort associations
- Independent replication is based on selected associations and model refitting in an independent cohort
- Multivariable predictor sets are defined explicitly rather than automatically selected according to replication-cohort p-values

---

## Requirements

R (≥ 4.0 recommended)

Key packages:
- dplyr
- tidyr
- tibble
- glmnet
- pROC
- metafor
- broom
- ggplot2
- scales

---

## Notes

- Scripts are intended to be run sequentially
- File paths must be adapted locally
- Raw data are not included in this repository
- Candidate predictor sets and cohort-specific preprocessing must be defined before modelling
- The independent replication script includes both univariable and multivariable analyses
- Multivariable replication models are fitted using standard logistic regression without penalisation
- Outputs include descriptive tables, regression estimates, model summaries, and calibration plots where applicable
- The repository provides adaptable analytical templates rather than an immediately executable reproduction of every cohort-specific analysis

---

## Contact

For questions or collaboration, please contact the repository owner.
