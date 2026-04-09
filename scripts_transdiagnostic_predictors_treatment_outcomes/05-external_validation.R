# ============================================================
# Project: Transdiagnostic predictors of treatment outcomes
# Script: External validation analyses
#
# Workflow:
# This script implements the external validation stage in two steps:
#
# 1) Univariable analyses
# 2) Multivariable analyses
#
# Notes:
# - this script assumes that the validation dataset has already
#   undergone preprocessing
# - any cohort-specific cleaning, exclusions, or handling of
#   outliers should be completed before this step
# ============================================================

# ============================================================
# PART 1: Univariable external validation
# ============================================================
#- summarise selected predictors by outcome status
#    - fit one univariable logistic regression per predictor
#    - this step follows the same workflow used in
#      01-univariable_analyses
#    - the only difference is that predictors are restricted to
#      those selected from the discovery meta-analyses, based on:
#        • pooled effects at p < 0.10
#        • consistency in effect direction in at least 2 out of
#          3 discovery cohorts
#
# See 01-univariable_analyses for implementation details.
# Only the predictor selection rule differs.

# ============================================================
# PART 2: Multivariable external validation
# ============================================================

#    - fit standard multivariable logistic regression models
#    - estimate joint effects of predictors carried forward from
#      the validation univariable step
#    - models are fitted using glm() (no elastic net, no tuning)

suppressPackageStartupMessages({
  library(dplyr)
  library(broom)
  library(tibble)
})

options(contrasts = c("contr.treatment", "contr.poly"))

# ------------------------------------------------------------
# Input dataset
# ------------------------------------------------------------
# Replace this object with the external replication dataset.
# The dataset should already contain one row per participant,
# binary outcomes coded as 0/1, and the predictors intended for
# the multivariable replication models.
# ------------------------------------------------------------

validation_df <- validation_cohort_df

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

resolve_first_present <- function(data, candidates, what = "variable") {
  hit <- candidates[candidates %in% names(data)]
  if (length(hit) == 0) {
    stop(
      "None of the candidate ", what, " names found: ",
      paste(candidates, collapse = ", ")
    )
  }
  hit[1]
}

prepare_binary_factor <- function(x) {
  if (is.factor(x)) {
    return(droplevels(x))
  }
  
  if (is.logical(x)) {
    x <- ifelse(is.na(x), NA, ifelse(x, "Yes", "No"))
    return(factor(x, levels = c("No", "Yes")))
  }
  
  if (is.numeric(x) || is.integer(x)) {
    ux <- sort(unique(x[!is.na(x)]))
    if (length(ux) == 2 && all(ux %in% c(0, 1))) {
      return(factor(ifelse(x == 1, "Yes", "No"), levels = c("No", "Yes")))
    }
  }
  
  x_chr <- as.character(x)
  x_low <- tolower(trimws(x_chr))
  
  yes_vals <- c("yes", "y", "1", "true", "t", "present", "positive")
  no_vals  <- c("no", "n", "0", "false", "f", "absent", "negative")
  
  out <- rep(NA_character_, length(x_low))
  out[x_low %in% yes_vals] <- "Yes"
  out[x_low %in% no_vals]  <- "No"
  
  if (all(is.na(out))) {
    ux <- unique(x_chr[!is.na(x_chr)])
    if (length(ux) == 2) {
      out[x_chr == ux[1]] <- ux[1]
      out[x_chr == ux[2]] <- ux[2]
      return(factor(out))
    }
    return(factor(x_chr))
  }
  
  factor(out, levels = c("No", "Yes"))
}

prepare_gender <- function(x) {
  if (is.factor(x)) {
    levs <- levels(x)
    if (all(c("Female", "Male") %in% levs)) {
      return(relevel(droplevels(x), ref = "Female"))
    }
    return(droplevels(x))
  }
  
  if (is.numeric(x) || is.integer(x)) {
    x2 <- ifelse(
      is.na(x), NA,
      ifelse(x %in% c(1), "Male",
             ifelse(x %in% c(0), "Female", NA))
    )
    return(factor(x2, levels = c("Female", "Male")))
  }
  
  x_chr <- trimws(as.character(x))
  x_low <- tolower(x_chr)
  
  x2 <- ifelse(
    x_low %in% c("male", "m", "man"), "Male",
    ifelse(x_low %in% c("female", "f", "woman"), "Female", x_chr)
  )
  
  factor(x2, levels = c("Female", "Male"))
}

prepare_ordered_3cat <- function(x) {
  if (is.numeric(x) || is.integer(x)) {
    ux <- sort(unique(x[!is.na(x)]))
    if (length(ux) == 3 && all(ux %in% c(1, 2, 3))) {
      return(ordered(
        ifelse(x == 1, "Low", ifelse(x == 2, "Moderate", "High")),
        levels = c("Low", "Moderate", "High")
      ))
    }
  }
  
  x_chr <- trimws(as.character(x))
  x_low <- tolower(x_chr)
  
  x_map <- ifelse(
    x_low %in% c("low", "1"), "Low",
    ifelse(x_low %in% c("moderate", "medium", "2"), "Moderate",
           ifelse(x_low %in% c("high", "3"), "High", x_chr))
  )
  
  vals <- unique(x_map[!is.na(x_map)])
  if (all(vals %in% c("Low", "Moderate", "High"))) {
    return(ordered(x_map, levels = c("Low", "Moderate", "High")))
  }
  
  ordered(x_chr)
}

wald_ci <- function(model) {
  b  <- coef(model)
  se <- sqrt(diag(vcov(model)))
  
  out <- cbind(
    CI_low  = b - 1.96 * se,
    CI_high = b + 1.96 * se
  )
  
  out
}

extract_model_table <- function(model) {
  sm <- summary(model)$coefficients
  ci <- wald_ci(model)
  
  terms_keep <- rownames(sm)
  terms_keep <- terms_keep[terms_keep != "(Intercept)"]
  
  tibble(
    term    = terms_keep,
    logOR   = as.numeric(sm[terms_keep, "Estimate"]),
    OR      = exp(logOR),
    CI_low  = exp(ci[terms_keep, "CI_low"]),
    CI_high = exp(ci[terms_keep, "CI_high"]),
    p.value = as.numeric(sm[terms_keep, "Pr(>|z|)"])
  )
}

print_model_block <- function(label, fit_obj, out_table) {
  cat("\n============================================================\n")
  cat(label, "\n")
  cat("============================================================\n")
  cat("Outcome: ", fit_obj$outcome_name, "\n", sep = "")
  cat("N total: ", fit_obj$n_total, "\n", sep = "")
  cat("N complete-case: ", fit_obj$n_cc, "\n", sep = "")
  cat("Events: ", fit_obj$events, " | Non-events: ", fit_obj$nonevents, "\n", sep = "")
  cat("Predictors entered: ", paste(fit_obj$predictors_used, collapse = ", "), "\n\n", sep = "")
  print(out_table, n = Inf, width = Inf)
  cat("\n")
}

# ------------------------------------------------------------
# Outcome resolution
# ------------------------------------------------------------
# Replace these aliases if the external validation dataset uses
# different outcome names.
# ------------------------------------------------------------

outcome_aliases <- list(
  Resistance = c("Resistance", "TR1", "TR2"),
  Remission  = c("remission", "Remission")
)

outcome_resolved <- list(
  Resistance = resolve_first_present(validation_df, outcome_aliases$Resistance, what = "Resistance outcome"),
  Remission  = resolve_first_present(validation_df, outcome_aliases$Remission,  what = "Remission outcome")
)


# ------------------------------------------------------------
# Predictor sets
# ------------------------------------------------------------
# These are example predictor sets for the external replication
# multivariable models.
#
# Replace the placeholder names below with the variables available
# in the validation dataset. Predictors are assumed to have been
# defined upstream, for example based on the validation univariable
# analyses or on a predefined replication plan.
# ------------------------------------------------------------

pred_vars_model_1 <- c(
  "AGE",
  "SEX",
  "EDUCATION_LEVEL_ORD",
  "ILLNESS_COURSE_INDEX",
  "BASELINE_SEVERITY",
  "ANXIETY_LEVEL_ORD",
  "DEPRESSIVE_FEATURES_ORD",
  "PSYCHOTIC_FEATURES_BIN",
  "SUICIDALITY_BIN",
  "FUNCTIONING_SCORE",
  "SOCIAL_WITHDRAWAL_ORD",
  "PRS_A_resid_z",
  "PRS_B_resid_z"
)

pred_vars_model_2 <- c(
  "AGE",
  "SEX",
  "EMPLOYMENT_STATUS",
  "AGE_AT_ONSET",
  "SUBSTANCE_USE_BIN",
  "PSYCHOTIC_FEATURES_BIN",
  "SUICIDALITY_BIN",
  "FUNCTIONING_SCORE",
  "SOCIAL_WITHDRAWAL_ORD"
)

pred_vars_model_3 <- c(
  "AGE",
  "SEX",
  "EDUCATION_LEVEL_ORD",
  "ILLNESS_COURSE_INDEX",
  "BASELINE_SEVERITY",
  "ANXIETY_LEVEL_ORD",
  "DEPRESSIVE_FEATURES_ORD",
  "PSYCHOTIC_FEATURES_BIN",
  "SUICIDALITY_BIN",
  "FUNCTIONING_SCORE",
  "SOCIAL_WITHDRAWAL_ORD",
  "PRS_A_resid_z",
  "PRS_B_resid_z"
)

model_specs <- list(
  model_1 = list(
    outcome_name = outcome_resolved$Outcome_1,
    predictors = pred_vars_model_1
  ),
  model_2 = list(
    outcome_name = outcome_resolved$Outcome_2,
    predictors = pred_vars_model_2
  ),
  model_3 = list(
    outcome_name = outcome_resolved$Outcome_3,
    predictors = pred_vars_model_3
  )
)

# ------------------------------------------------------------
# Prepare data types
# ------------------------------------------------------------

df <- validation_df

# Outcomes as 0/1 numeric if needed
for (out_nm in unique(vapply(model_specs, function(x) x$outcome_name, character(1)))) {
  if (is.factor(df[[out_nm]])) {
    levs <- levels(df[[out_nm]])
    if (length(levs) == 2) {
      df[[out_nm]] <- as.integer(df[[out_nm]] == levs[2])
    }
  }
}

# Sex / gender
if ("SEX" %in% names(df)) {
  df$SEX <- prepare_gender(df$SEX)
}

# Binary predictors
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

# Ordered 3-category predictors
ordered_3cat_vars <- c(
  "EDUCATION_LEVEL_ORD",
  "ANXIETY_LEVEL_ORD",
  "DEPRESSIVE_FEATURES_ORD",
  "SOCIAL_WITHDRAWAL_ORD"
)

for (v in intersect(ordered_3cat_vars, names(df))) {
  df[[v]] <- prepare_ordered_3cat(df[[v]])
}

# ------------------------------------------------------------
# Model fitter
# ------------------------------------------------------------

fit_multivariable_model <- function(data, outcome_var, predictors, model_label) {
  
  predictors_in <- predictors[predictors %in% names(data)]
  missing_pred  <- setdiff(predictors, predictors_in)
  
  if (length(missing_pred) > 0) {
    cat(
      model_label, " | Missing predictors dropped: ",
      paste(missing_pred, collapse = ", "), "\n", sep = ""
    )
  }
  
  if (!(outcome_var %in% names(data))) {
    stop(model_label, ": outcome not found -> ", outcome_var)
  }
  
  if (length(predictors_in) == 0) {
    stop(model_label, ": no predictors available after missing check.")
  }
  
  vars_needed <- c(outcome_var, predictors_in)
  
  dsub <- data %>%
    dplyr::select(all_of(vars_needed)) %>%
    dplyr::filter(complete.cases(.))
  
  if (nrow(dsub) < 30) {
    stop(model_label, ": too few complete-case observations (n = ", nrow(dsub), ").")
  }
  
  y <- dsub[[outcome_var]]
  
  if (!all(na.omit(unique(y)) %in% c(0, 1))) {
    stop(model_label, ": outcome must be coded as 0/1 after preprocessing.")
  }
  
  if (length(unique(y)) < 2) {
    stop(model_label, ": outcome has no variation in complete-case data.")
  }
  
  form <- as.formula(
    paste0(outcome_var, " ~ ", paste(sprintf("`%s`", predictors_in), collapse = " + "))
  )
  
  fit <- glm(form, data = dsub, family = binomial())
  out_table <- extract_model_table(fit)
  
  list(
    label = model_label,
    outcome_name = outcome_var,
    predictors_used = predictors_in,
    missing_predictors = missing_pred,
    n_total = nrow(data),
    n_cc = nrow(dsub),
    events = sum(y == 1),
    nonevents = sum(y == 0),
    fit = fit,
    table = out_table
  )
}

# ------------------------------------------------------------
# Run
# ------------------------------------------------------------

model_objects <- list()

for (nm in names(model_specs)) {
  spec <- model_specs[[nm]]
  
  obj <- fit_multivariable_model(
    data = df,
    outcome_var = spec$outcome_name,
    predictors = spec$predictors,
    model_label = nm
  )
  
  model_objects[[nm]] <- obj
  print_model_block(nm, obj, obj$table)
}

# model_objects contains the fitted models and the extracted tables