
# ============================================================
# Project: Transdiagnostic predictors of treatment outcomes
# Script: Independent replication of treatment-outcome associations
#
# Purpose:
# - replicate selected univariable associations in an independent cohort
# - summarise candidate predictors by outcome status
# - refit multivariable logistic regression models
# - report odds ratios, 95% confidence intervals and p-values
# - estimate in-sample discriminative performance (AUC)
#
# Workflow:
# 1) Univariable independent replication
# 2) Multivariable independent replication
# 3) Export and summarise results
#
# Notes:
# - candidate predictors are selected upstream from the
#   discovery meta-analyses
# - univariable associations are tested without additional
#   outcome-based predictor selection in the replication cohort
# - separate TR1-derived and TR2-derived predictor sets can
#   be tested against a common resistance outcome
# - multivariable models are refitted in the independent cohort
# - original discovery models are not transferred unchanged
# - cohort-specific preprocessing must precede this script
# ============================================================


# ============================================================
# Packages
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(pROC)
})

options(contrasts = c("contr.treatment", "contr.poly"))


# ============================================================
# Input dataset
# ============================================================
# Replace with the preprocessed independent replication dataset.
#
# Required structure:
# - one row per participant
# - binary outcomes coded as 0/1
# - harmonised clinical and sociodemographic predictors
# - available polygenic scores
#
# Cohort-specific exclusions, missing-value recoding,
# transformations and harmonisation must be completed upstream.
# ============================================================

replication_df <- replication_cohort_df


# ============================================================
# Settings
# ============================================================

outcome_remission  <- "Remission"
outcome_resistance <- "Resistance"

# Abort when a requested predictor is absent.

stop_on_missing_predictors <- TRUE

# Minimum complete-case sample size for logistic models.

min_complete_cases <- 30

# Directory for output tables.

output_dir <- "Outputs"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# Predictor sets
# ============================================================
# Replace the illustrative variables below with the actual
# predictors selected from the discovery meta-analyses.
#
# Selection criteria used upstream:
# - pooled association p < 0.10
# - concordant effect direction in at least two of the
#   three discovery cohorts
#
# The same candidate sets are initially used for univariable
# replication and multivariable model specification.
#
# If multivariable models use a narrower predefined subset,
# modify multivariable_predictor_sets separately below.
#
# The examples do not reproduce the exact published models.
# ============================================================

pred_vars_remission <- c(
  "AGE",
  "SEX",
  "EMPLOYMENT_STATUS",
  "EDUCATION_LEVEL_ORD",
  "FUNCTIONING_SCORE",
  "ANXIETY_LEVEL_ORD",
  "DEPRESSIVE_FEATURES_ORD",
  "PSYCHOTIC_FEATURES_BIN",
  "SUICIDALITY_BIN",
  "SOCIAL_WITHDRAWAL_ORD"
)

pred_vars_tr1 <- c(
  "AGE",
  "SEX",
  "EMPLOYMENT_STATUS",
  "EDUCATION_LEVEL_ORD",
  "FUNCTIONING_SCORE",
  "BASELINE_SEVERITY",
  "ANXIETY_LEVEL_ORD",
  "PSYCHOTIC_FEATURES_BIN",
  "SUICIDALITY_BIN"
)

pred_vars_tr2 <- c(
  "AGE",
  "SEX",
  "EMPLOYMENT_STATUS",
  "FUNCTIONING_SCORE",
  "BASELINE_SEVERITY",
  "SOCIAL_WITHDRAWAL_ORD",
  "PRS_EDUCATION_resid_z"
)


# ============================================================
# Outcome and predictor mappings
# ============================================================
# Remission-derived predictors are tested against remission.
#
# TR1-derived and TR2-derived predictors are both tested
# against the same resistance outcome.
#
# Separate model labels preserve the discovery origin
# of each predictor set.
# ============================================================

replication_specs <- list(

  Remission = list(
    outcome = outcome_remission,
    predictors = pred_vars_remission
  ),

  TR1_derived = list(
    outcome = outcome_resistance,
    predictors = pred_vars_tr1
  ),

  TR2_derived = list(
    outcome = outcome_resistance,
    predictors = pred_vars_tr2
  )

)

# The multivariable predictor sets can be specified independently
# without changing the univariable replication candidates.

multivariable_predictor_sets <- list(
  Remission = pred_vars_remission,
  TR1_derived = pred_vars_tr1,
  TR2_derived = pred_vars_tr2
)


# ============================================================
# Helper functions
# ============================================================

# ------------------------------------------------------------
# Prepare binary outcomes
# ------------------------------------------------------------
# Outcomes must be numeric 0/1 or logical.
# The event of interest must be coded as 1.
# ------------------------------------------------------------

prepare_outcome <- function(x, outcome_name) {

  if (is.logical(x)) {
    x <- as.integer(x)
  }

  if (!is.numeric(x) && !is.integer(x)) {
    stop(
      outcome_name,
      ": outcome must be coded as numeric 0/1 or logical."
    )
  }

  vals <- unique(x[!is.na(x)])

  if (!all(vals %in% c(0, 1))) {
    stop(
      outcome_name,
      ": outcome contains values other than 0 and 1."
    )
  }

  as.integer(x)
}


# ------------------------------------------------------------
# Prepare binary predictors
# ------------------------------------------------------------
# Recognised coding:
# No = reference category
# Yes = comparison category
# ------------------------------------------------------------

prepare_binary_factor <- function(x) {

  x_chr <- tolower(trimws(as.character(x)))

  yes_vals <- c(
    "yes", "y", "1", "true", "t",
    "present", "positive"
  )

  no_vals <- c(
    "no", "n", "0", "false", "f",
    "absent", "negative"
  )

  valid <- is.na(x_chr) |
    x_chr %in% c(yes_vals, no_vals)

  if (!all(valid)) {
    stop(
      "Unrecognised binary predictor categories: ",
      paste(unique(x_chr[!valid]), collapse = ", ")
    )
  }

  out <- rep(NA_character_, length(x_chr))

  out[x_chr %in% yes_vals] <- "Yes"
  out[x_chr %in% no_vals] <- "No"

  factor(out, levels = c("No", "Yes"))
}


# ------------------------------------------------------------
# Prepare sex
# ------------------------------------------------------------
# Example coding:
# 0 = Female
# 1 = Male
#
# Verify this coding before analysing a particular cohort.
# ------------------------------------------------------------

prepare_sex <- function(x) {

  x_chr <- tolower(trimws(as.character(x)))

  female_vals <- c("female", "f", "woman", "0")
  male_vals <- c("male", "m", "man", "1")

  valid <- is.na(x_chr) |
    x_chr %in% c(female_vals, male_vals)

  if (!all(valid)) {
    stop("Unrecognised sex categories.")
  }

  out <- rep(NA_character_, length(x_chr))

  out[x_chr %in% female_vals] <- "Female"
  out[x_chr %in% male_vals] <- "Male"

  factor(out, levels = c("Female", "Male"))
}


# ------------------------------------------------------------
# Prepare ordered three-category predictors
# ------------------------------------------------------------
# Low < Moderate < High
#
# Polynomial contrasts:
# .L = linear component
# .Q = quadratic component
#
# exp(beta.L) is not necessarily the direct
# High-versus-Low odds ratio.
# ------------------------------------------------------------

prepare_ordered_3cat <- function(x) {

  x_chr <- tolower(trimws(as.character(x)))

  out <- rep(NA_character_, length(x_chr))

  out[x_chr %in% c("low", "1")] <- "Low"

  out[x_chr %in% c(
    "moderate", "medium", "2"
  )] <- "Moderate"

  out[x_chr %in% c("high", "3")] <- "High"

  valid <- is.na(x_chr) | !is.na(out)

  if (!all(valid)) {
    stop(
      "Unrecognised ordinal categories: ",
      paste(unique(x_chr[!valid]), collapse = ", ")
    )
  }

  ordered(
    out,
    levels = c("Low", "Moderate", "High")
  )
}


# ------------------------------------------------------------
# Check requested predictors
# ------------------------------------------------------------

resolve_predictors <- function(
    data,
    predictors,
    label,
    stop_on_missing = TRUE
) {

  predictors <- unique(predictors)

  missing_pred <- setdiff(
    predictors,
    names(data)
  )

  if (length(missing_pred) > 0) {

    msg <- paste0(
      label,
      " | Missing predictors: ",
      paste(missing_pred, collapse = ", ")
    )

    if (isTRUE(stop_on_missing)) {
      stop(msg)
    } else {
      warning(msg)
    }
  }

  predictors_in <- intersect(
    predictors,
    names(data)
  )

  if (length(predictors_in) == 0) {
    stop(label, ": no predictors available.")
  }

  list(
    available = predictors_in,
    missing = missing_pred
  )
}


# ------------------------------------------------------------
# Extract logistic regression results
# ------------------------------------------------------------
# Coefficient-level OR, 95% Wald CI and p-values.
#
# Categorical predictors may generate multiple coefficients.
# For ordered factors, .L and .Q are reported separately.
# ------------------------------------------------------------

extract_model_table <- function(model) {

  sm <- summary(model)$coefficients

  terms_keep <- setdiff(
    rownames(sm),
    "(Intercept)"
  )

  beta <- as.numeric(
    sm[terms_keep, "Estimate"]
  )

  se <- as.numeric(
    sm[terms_keep, "Std. Error"]
  )

  tibble(
    term = terms_keep,
    logOR = beta,
    SE = se,
    OR = exp(beta),
    CI_low = exp(beta - 1.96 * se),
    CI_high = exp(beta + 1.96 * se),
    p.value = as.numeric(
      sm[terms_keep, "Pr(>|z|)"]
    )
  )
}


# ------------------------------------------------------------
# Overall predictor p-value
# ------------------------------------------------------------
# Likelihood-ratio comparison between a model containing
# the predictor and an intercept-only model.
#
# This is especially useful for categorical predictors
# with more than one regression coefficient.
# ------------------------------------------------------------

overall_predictor_test <- function(model) {

  null_model <- update(
    model,
    formula = . ~ 1
  )

  comparison <- anova(
    null_model,
    model,
    test = "Chisq"
  )

  as.numeric(
    comparison[2, "Pr(>Chi)"]
  )
}


# ------------------------------------------------------------
# Calculate in-sample AUC
# ------------------------------------------------------------

calculate_insample_auc <- function(y, predicted) {

  stopifnot(
    length(y) == length(predicted)
  )

  roc_obj <- pROC::roc(
    response = y,
    predictor = predicted,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )

  auc_ci <- as.numeric(
    pROC::ci.auc(roc_obj)
  )

  list(
    AUC = as.numeric(pROC::auc(roc_obj)),
    AUC_low = auc_ci[1],
    AUC_high = auc_ci[3],
    roc = roc_obj
  )
}


# ============================================================
# Prepare analytic dataset
# ============================================================

df <- replication_df

required_outcomes <- unique(
  vapply(
    replication_specs,
    function(x) x$outcome,
    character(1)
  )
)

missing_outcomes <- setdiff(
  required_outcomes,
  names(df)
)

if (length(missing_outcomes) > 0) {
  stop(
    "Missing outcomes: ",
    paste(missing_outcomes, collapse = ", ")
  )
}

for (out_nm in required_outcomes) {

  df[[out_nm]] <- prepare_outcome(
    df[[out_nm]],
    out_nm
  )
}


# ------------------------------------------------------------
# Predictor type preparation
# ------------------------------------------------------------
# Adapt these lists to the cohort-specific data dictionary.
# ------------------------------------------------------------

if ("SEX" %in% names(df)) {
  df$SEX <- prepare_sex(df$SEX)
}

binary_vars <- c(
  "EMPLOYMENT_STATUS",
  "LIVING_STATUS_BIN",
  "FAMILY_HISTORY_BIN",
  "COMORBIDITY_BIN",
  "SUBSTANCE_USE_BIN",
  "PSYCHOTIC_FEATURES_BIN",
  "SUICIDALITY_BIN"
)

for (v in intersect(binary_vars, names(df))) {
  df[[v]] <- prepare_binary_factor(df[[v]])
}

ordered_3cat_vars <- c(
  "EDUCATION_LEVEL_ORD",
  "ANXIETY_LEVEL_ORD",
  "DEPRESSIVE_FEATURES_ORD",
  "SOCIAL_WITHDRAWAL_ORD"
)

for (v in intersect(ordered_3cat_vars, names(df))) {
  df[[v]] <- prepare_ordered_3cat(df[[v]])
}


# ============================================================
# PART 1: Univariable independent replication
# ============================================================
# For each discovery-derived predictor:
#
# - summarise its distribution by outcome status
# - fit a separate univariable logistic regression
# - estimate OR, 95% CI and p-value
# - calculate an overall predictor p-value
#
# Analyses use the available complete-case sample for
# each outcome-predictor pair.
#
# No additional p-value threshold is applied to determine
# entry into the predefined multivariable models.
# ============================================================


# ------------------------------------------------------------
# Descriptive summaries by outcome status
# ------------------------------------------------------------
# Numeric predictors:
# - N, mean, SD, median, Q1 and Q3
#
# Categorical predictors:
# - counts and percentages within outcome group
#
# Missing predictor values are excluded from summaries;
# N_available indicates the non-missing denominator.
# ------------------------------------------------------------

describe_predictor <- function(
    data,
    outcome_var,
    predictor_var,
    model_label
) {

  d <- data %>%
    select(
      all_of(c(outcome_var, predictor_var))
    ) %>%
    filter(!is.na(.data[[outcome_var]]))

  x <- d[[predictor_var]]
  y <- d[[outcome_var]]

  results <- list()

  for (group in c(0, 1)) {

    x_group <- x[y == group]
    x_valid <- x_group[!is.na(x_group)]

    n_available <- length(x_valid)

    if (is.numeric(x) && !is.factor(x)) {

      results[[length(results) + 1]] <- tibble(
        Model = model_label,
        Outcome = outcome_var,
        Predictor = predictor_var,
        Outcome_status = group,
        Type = "Continuous",
        Category = NA_character_,
        N_available = n_available,
        N_category = NA_integer_,
        Percentage = NA_real_,
        Mean = if (n_available > 0) mean(x_valid) else NA_real_,
        SD = if (n_available > 1) sd(x_valid) else NA_real_,
        Median = if (n_available > 0) median(x_valid) else NA_real_,
        Q1 = if (n_available > 0) {
          as.numeric(quantile(x_valid, 0.25))
        } else {
          NA_real_
        },
        Q3 = if (n_available > 0) {
          as.numeric(quantile(x_valid, 0.75))
        } else {
          NA_real_
        }
      )

    } else {

      categories <- if (is.factor(x)) {
        levels(x)
      } else {
        sort(unique(as.character(x_valid)))
      }

      for (category in categories) {

        n_category <- sum(
          as.character(x_valid) == category
        )

        results[[length(results) + 1]] <- tibble(
          Model = model_label,
          Outcome = outcome_var,
          Predictor = predictor_var,
          Outcome_status = group,
          Type = "Categorical",
          Category = as.character(category),
          N_available = n_available,
          N_category = n_category,
          Percentage = if (n_available > 0) {
            100 * n_category / n_available
          } else {
            NA_real_
          },
          Mean = NA_real_,
          SD = NA_real_,
          Median = NA_real_,
          Q1 = NA_real_,
          Q3 = NA_real_
        )
      }
    }
  }

  bind_rows(results)
}


# ------------------------------------------------------------
# Fit one univariable logistic regression
# ------------------------------------------------------------

fit_univariable_model <- function(
    data,
    outcome_var,
    predictor_var,
    model_label,
    min_n = 30
) {

  dsub <- data %>%
    select(
      all_of(c(outcome_var, predictor_var))
    ) %>%
    filter(complete.cases(.)) %>%
    droplevels()

  if (nrow(dsub) < min_n) {
    stop(
      "Too few complete cases for ",
      predictor_var,
      " (N = ",
      nrow(dsub),
      ")."
    )
  }

  y <- dsub[[outcome_var]]
  x <- dsub[[predictor_var]]

  if (length(unique(y)) < 2) {
    stop("Outcome has no variation.")
  }

  if (length(unique(x)) < 2) {
    stop("Predictor has no variation.")
  }

  form <- reformulate(
    termlabels = predictor_var,
    response = outcome_var
  )

  fit <- glm(
    formula = form,
    data = dsub,
    family = binomial()
  )

  if (!isTRUE(fit$converged)) {
    warning(
      model_label,
      " | ",
      predictor_var,
      ": model did not converge."
    )
  }

  result <- extract_model_table(fit)

  result %>%
    mutate(
      Model = model_label,
      Outcome = outcome_var,
      Predictor = predictor_var,
      N_cc = nrow(dsub),
      Events = sum(y == 1),
      Nonevents = sum(y == 0),
      Overall_p = overall_predictor_test(fit),
      .before = 1
    )
}


# ------------------------------------------------------------
# Run univariable replication analyses
# ------------------------------------------------------------

univariable_rows <- list()
descriptive_rows <- list()
univariable_errors <- list()

for (model_label in names(replication_specs)) {

  spec <- replication_specs[[model_label]]

  outcome_var <- spec$outcome

  resolved <- resolve_predictors(
    data = df,
    predictors = spec$predictors,
    label = model_label,
    stop_on_missing = stop_on_missing_predictors
  )

  for (predictor_var in resolved$available) {

    descriptive_rows[[length(descriptive_rows) + 1]] <-
      describe_predictor(
        data = df,
        outcome_var = outcome_var,
        predictor_var = predictor_var,
        model_label = model_label
      )

    result <- tryCatch(

      fit_univariable_model(
        data = df,
        outcome_var = outcome_var,
        predictor_var = predictor_var,
        model_label = model_label,
        min_n = min_complete_cases
      ),

      error = function(e) {

        univariable_errors[[length(univariable_errors) + 1]] <<-
          tibble(
            Model = model_label,
            Predictor = predictor_var,
            Error = conditionMessage(e)
          )

        NULL
      }
    )

    if (!is.null(result)) {

      univariable_rows[[length(univariable_rows) + 1]] <-
        result
    }
  }
}

univariable_results <- bind_rows(univariable_rows)

descriptive_results <- bind_rows(descriptive_rows)

univariable_error_log <- bind_rows(univariable_errors)


# ============================================================
# Print univariable replication results
# ============================================================

cat("\n============================================================\n")
cat("UNIVARIABLE INDEPENDENT REPLICATION\n")
cat("============================================================\n")

print(
  univariable_results,
  n = Inf,
  width = Inf
)

cat("\nDESCRIPTIVE SUMMARIES\n")

print(
  descriptive_results,
  n = Inf,
  width = Inf
)

if (nrow(univariable_error_log) > 0) {

  cat("\nUNIVARIABLE MODEL ERRORS\n")

  print(
    univariable_error_log,
    n = Inf,
    width = Inf
  )
}


# ============================================================
# PART 2: Multivariable independent replication
# ============================================================
# The candidate sets are defined upstream and are not
# automatically filtered according to their univariable
# p-values in the replication cohort.
#
# Three models are fitted:
# - Remission
# - TR1-derived resistance
# - TR2-derived resistance
#
# Standard logistic regression is used throughout.
# ============================================================


# ------------------------------------------------------------
# Fit one multivariable logistic regression model
# ------------------------------------------------------------

fit_multivariable_model <- function(
    data,
    outcome_var,
    predictors,
    model_label,
    stop_on_missing = TRUE,
    min_n = 30
) {

  resolved <- resolve_predictors(
    data = data,
    predictors = predictors,
    label = model_label,
    stop_on_missing = stop_on_missing
  )

  predictors_in <- resolved$available

  dsub <- data %>%
    select(
      all_of(c(outcome_var, predictors_in))
    ) %>%
    filter(complete.cases(.)) %>%
    droplevels()

  if (nrow(dsub) < min_n) {
    stop(
      model_label,
      ": too few complete-case observations (N = ",
      nrow(dsub),
      ")."
    )
  }

  y <- dsub[[outcome_var]]

  if (length(unique(y)) < 2) {
    stop(
      model_label,
      ": outcome has no variation."
    )
  }

  for (v in predictors_in) {

    if (length(unique(dsub[[v]])) < 2) {
      stop(
        model_label,
        ": predictor has no variation: ",
        v
      )
    }
  }

  form <- reformulate(
    termlabels = predictors_in,
    response = outcome_var
  )

  fit <- glm(
    formula = form,
    data = dsub,
    family = binomial()
  )

  if (!isTRUE(fit$converged)) {
    warning(
      model_label,
      ": logistic regression did not converge."
    )
  }

  out_table <- extract_model_table(fit)

  predicted <- as.numeric(
    predict(
      fit,
      newdata = dsub,
      type = "response"
    )
  )

  auc_result <- calculate_insample_auc(
    y = y,
    predicted = predicted
  )

  list(
    label = model_label,
    outcome_name = outcome_var,
    predictors_used = predictors_in,
    n_total = nrow(data),
    n_cc = nrow(dsub),
    events = sum(y == 1),
    nonevents = sum(y == 0),
    AUC = auc_result$AUC,
    AUC_low = auc_result$AUC_low,
    AUC_high = auc_result$AUC_high,
    fit = fit,
    table = out_table,
    roc = auc_result$roc,
    predicted = predicted
  )
}


# ------------------------------------------------------------
# Run multivariable replication models
# ------------------------------------------------------------

model_objects <- list()

for (model_label in names(replication_specs)) {

  spec <- replication_specs[[model_label]]

  predictors <- multivariable_predictor_sets[[model_label]]

  obj <- fit_multivariable_model(
    data = df,
    outcome_var = spec$outcome,
    predictors = predictors,
    model_label = model_label,
    stop_on_missing = stop_on_missing_predictors,
    min_n = min_complete_cases
  )

  model_objects[[model_label]] <- obj

  cat("\n============================================================\n")
  cat("MULTIVARIABLE INDEPENDENT REPLICATION\n")
  cat("Model: ", model_label, "\n", sep = "")
  cat("Outcome: ", obj$outcome_name, "\n", sep = "")
  cat("N complete-case: ", obj$n_cc, "\n", sep = "")

  cat(
    "Events: ", obj$events,
    " | Non-events: ", obj$nonevents,
    "\n",
    sep = ""
  )

  cat(
    "In-sample AUC: ",
    sprintf("%.3f", obj$AUC),
    " [",
    sprintf("%.3f", obj$AUC_low),
    ", ",
    sprintf("%.3f", obj$AUC_high),
    "]\n",
    sep = ""
  )

  print(
    obj$table,
    n = Inf,
    width = Inf
  )
}


# ============================================================
# PART 3: Summary and export
# ============================================================

# ------------------------------------------------------------
# Multivariable model summary
# ------------------------------------------------------------

replication_summary <- bind_rows(
  lapply(model_objects, function(obj) {

    tibble(
      Model = obj$label,
      Outcome = obj$outcome_name,
      N_cc = obj$n_cc,
      Events = obj$events,
      Nonevents = obj$nonevents,
      AUC = obj$AUC,
      AUC_low = obj$AUC_low,
      AUC_high = obj$AUC_high
    )
  })
)


# ------------------------------------------------------------
# Multivariable coefficient tables
# ------------------------------------------------------------

multivariable_results <- bind_rows(
  lapply(model_objects, function(obj) {

    obj$table %>%
      mutate(
        Model = obj$label,
        Outcome = obj$outcome_name,
        N_cc = obj$n_cc,
        .before = 1
      )
  })
)


# ------------------------------------------------------------
# Export CSV files
# ------------------------------------------------------------

write.csv(
  descriptive_results,
  file.path(
    output_dir,
    "replication_descriptive_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  univariable_results,
  file.path(
    output_dir,
    "replication_univariable_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  multivariable_results,
  file.path(
    output_dir,
    "replication_multivariable_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  replication_summary,
  file.path(
    output_dir,
    "replication_model_summary.csv"
  ),
  row.names = FALSE
)

if (nrow(univariable_error_log) > 0) {

  write.csv(
    univariable_error_log,
    file.path(
      output_dir,
      "replication_univariable_errors.csv"
    ),
    row.names = FALSE
  )
}


# ============================================================
# Final summary
# ============================================================

cat("\n============================================================\n")
cat("INDEPENDENT REPLICATION SUMMARY\n")
cat("============================================================\n")

print(
  replication_summary,
  n = Inf,
  width = Inf
)

cat(
  "\nResults exported to: ",
  output_dir,
  "\n",
  sep = ""
)


# ============================================================
# Stored objects
# ============================================================
#
# descriptive_results:
# - predictor distributions by outcome status
#
# univariable_results:
# - univariable OR, CI, p-values
# - overall predictor p-values
#
# multivariable_results:
# - adjusted OR, CI and p-values
#
# replication_summary:
# - complete-case sample sizes and in-sample AUC
#
# model_objects:
# - fitted multivariable glm objects
# - predicted probabilities and ROC objects
#
# Examples:
#
# model_objects[["Remission"]]$table
# model_objects[["TR1_derived"]]$AUC
# model_objects[["TR2_derived"]]$fit
#
# ============================================================

# ============================================================
# End of independent replication workflow
# ============================================================
