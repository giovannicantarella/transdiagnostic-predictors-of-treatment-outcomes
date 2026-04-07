# Transdiagnostic Predictors of Treatment Outcomes

This repository contains the scripts used to identify, model, and validate transdiagnostic clinical and genetic predictors of treatment outcomes (remission and treatment resistance) across major psychiatric disorders.

## Repository structure

The repository is organised into five main directories reflecting the analytical workflow:

scripts_transdiagnostic_predictors_treatment_outcomes/
├── preprocessing/
├── univariable/
├── multivariable/
├── meta-analyses/
└── validation/

## Preprocessing

Scripts for data cleaning, harmonisation, and variable construction across cohorts.

Files:
- preprocessing_SCZ_cohort.R
- preprocessing_BD_cohort.R
- preprocessing_MDD_cohort.R

Main steps:
- Variable harmonisation across cohorts
- Recoding of categorical predictors
- Creation of derived variables
- Handling of missing data
- Preparation of analysis-ready datasets

## Univariable analyses

Scripts performing univariable logistic regression analyses to screen predictors associated with treatment outcomes.

Files:
- Univariable_SCZ_cohort.R
- Univariable_BD_cohort.R
- Univariable_MDD_cohort.R

Main steps:
- Logistic regression models for each predictor
- Estimation of odds ratios (OR), 95% confidence intervals, and p-values
- Selection of predictors based on predefined significance thresholds

## Multivariable models

Scripts implementing multivariable predictive models within each cohort.

Files:
- Multivariable_SCZ_cohort.R
- Multivariable_BD_cohort.R
- Multivariable_MDD_cohort.R

Main features:
- Elastic Net logistic regression
- Nested cross-validation framework
- Hyperparameter tuning (alpha, lambda.1se)
- Out-of-fold AUC estimation
- Model selection based on performance and parsimony

## Meta-analyses

Scripts combining results across cohorts to identify robust transdiagnostic predictors.

Files:
- meta-analyses_remission_across_cohorts.R
- meta-analyses_TR1_across_cohorts.R
- meta-analyses_TR2_across_cohorts.R
- univariable_age_sex_adjusted_meta.R
- univariable_age_sex_adjusted_for_meta-analyses.R

Main steps:
- Age- and sex-adjusted univariable models
- Random-effects meta-analysis
- Estimation of pooled effect sizes (logOR and OR)
- Assessment of heterogeneity (I², τ²)

## Validation

Scripts for external validation of predictors in independent cohorts.

Files:
- preprocessing_MDD_validation_cohort.R
- univariable_validation_MDD_cohort.R
- multivariable_validation_MDD_cohort.R

Main steps:
- Harmonisation of the validation dataset
- Replication of univariable associations
- Application of multivariable models
- Evaluation of predictive performance (AUC, effect estimates)

## Analytical framework

The analytical workflow follows these steps:

1. Data preprocessing and harmonisation
2. Univariable screening of predictors
3. Multivariable modelling within each cohort
4. Meta-analysis across cohorts
5. External validation in independent samples

## Outcomes

The analyses focus on the following outcomes:
- Remission
- Treatment resistance (TR1: treatment-based; TR2: symptom-based)

## Requirements

R (version 4.0 or higher recommended)

Key packages:
- glmnet
- pROC
- dplyr
- tidyverse
- metafor

## Notes

- Scripts are intended to be run sequentially following the directory structure.
- File paths should be adapted to the local environment.
- Raw data are not included in this repository.

## Author

Giovanni Cantarella  
Research Fellow  
Department of Biomedical and Neuromotor Sciences (DIBINEM)  
University of Bologna, Italy
