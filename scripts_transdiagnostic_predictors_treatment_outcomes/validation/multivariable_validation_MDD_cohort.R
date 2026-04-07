# ============================================================
# STAR*D external validation | Multivariable logistic regression
# - Standard (non-penalized) logistic regression via glm()
# - Includes only predictors that reached p < 0.10
#   in the previous univariable validation analyses
#   (see: univariable_validation_MDD_cohort.R)
# - Three models:
#     1) Resistance (treatment-based phenotype) - TR1
#     2) Resistance (symptom-based phenotype) - TR2
#     3) Remission
# - Ordinal 3-category predictors are treated as ordered factors
# - Output: console tables with
#     term, logOR, OR, CI_low, CI_high, p.value
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(broom)
  library(tibble)
})

options(contrasts = c("contr.treatment", "contr.poly"))

# ------------------------------------------------------------
# PRECHECK
# ------------------------------------------------------------
if (!exists("STAR")) {
  stop("Object 'STAR' not found in the environment. Load the STAR dataset first.")
}

df0 <- STAR

# ------------------------------------------------------------
# HELPERS
# ------------------------------------------------------------
resolve_first_present <- function(data, candidates, what = "variable") {
  hit <- candidates[candidates %in% names(data)]
  if (length(hit) == 0) {
    stop("None of the candidate ", what, " names found: ",
         paste(candidates, collapse = ", "))
  }
  hit[1]
}

prepare_binary_factor <- function(x, positive_first = FALSE) {
  if (is.factor(x)) {
    x <- droplevels(x)
    return(x)
  }
  
  if (is.logical(x)) {
    x <- ifelse(is.na(x), NA, ifelse(x, "Yes", "No"))
    x <- factor(x, levels = c("No", "Yes"))
    return(x)
  }
  
  if (is.numeric(x) || is.integer(x)) {
    ux <- sort(unique(x[!is.na(x)]))
    if (length(ux) == 2 && all(ux %in% c(0, 1))) {
      x <- factor(ifelse(x == 1, "Yes", "No"), levels = c("No", "Yes"))
      return(x)
    }
  }
  
  x_chr <- as.character(x)
  x_low <- tolower(trimws(x_chr))
  
  yes_vals <- c("yes", "y", "1", "true", "t", "male", "m", "present", "positive")
  no_vals  <- c("no", "n", "0", "false", "f", "female", "woman", "absent", "negative")
  
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
      return(relevel(x, ref = "Female"))
    }
    return(droplevels(x))
  }
  
  if (is.numeric(x) || is.integer(x)) {
    x2 <- ifelse(is.na(x), NA,
                 ifelse(x %in% c(1), "Male",
                        ifelse(x %in% c(0), "Female", NA)))
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
# OUTCOME RESOLUTION
# ------------------------------------------------------------
outcome_aliases <- list(
  Resistance = c("Resistance", "TRD_vs_resp"),
  Remission  = c("remission", "QIDS_w12_rem")
)

outcome_resolved <- list(
  Resistance = resolve_first_present(df0, outcome_aliases$Resistance, what = "Resistance outcome"),
  Remission  = resolve_first_present(df0, outcome_aliases$Remission,  what = "Remission outcome")
)

# ------------------------------------------------------------
# PREDICTOR SETS
# These are the multivariable predictor sets used for the
# STAR*D replication models, corresponding to predictors that
# reached p < 0.10 in the previous univariable validation step.
#
# If one STAR column has a slightly different name in your object,
# edit ONLY these vectors.
# ------------------------------------------------------------
pred_vars_Resistance_TR1 <- c(
  "Age",
  "Gender",
  "Education_3cat",
  "N_previous_episodes_log",
  "HRSD_total",
  "Anxiety_3cat",
  "Depressive_symptoms_3cat",
  "Psychotic_symptoms_bin",
  "Suicidal_ideation_bin",
  "WSAS_baseline",
  "Social_Withdrawal_3cat",
  "DEPR07_PRS_resid_z",
  "EDUC03_PRS_resid_z"
)

pred_vars_Resistance_TR2 <- c(
  "Age",
  "Gender",
  "Employed",
  "Age_at_onset",
  "Substance_Use_Problem",
  "Psychotic_symptoms_bin",
  "Suicidal_ideation_bin",
  "WSAS_baseline",
  "Social_Withdrawal_3cat"
)

pred_vars_Remission <- c(
  "Age",
  "Gender",
  "Education_3cat",
  "N_previous_episodes_log",
  "HRSD_total",
  "Anxiety_3cat",
  "Depressive_symptoms_3cat",
  "Psychotic_symptoms_bin",
  "Suicidal_ideation_bin",
  "WSAS_baseline",
  "Social_Withdrawal_3cat",
  "DEPR07_PRS_resid_z",
  "EDUC03_PRS_resid_z"
)

model_specs <- list(
  Resistance_TR1 = list(
    outcome_name = outcome_resolved$Resistance,
    predictors = pred_vars_Resistance_TR1
  ),
  Resistance_TR2 = list(
    outcome_name = outcome_resolved$Resistance,
    predictors = pred_vars_Resistance_TR2
  ),
  Remission = list(
    outcome_name = outcome_resolved$Remission,
    predictors = pred_vars_Remission
  )
)

# ------------------------------------------------------------
# PREPARE DATA TYPES
# ------------------------------------------------------------
df <- df0

# outcomes as 0/1 numeric if needed
for (out_nm in unique(vapply(model_specs, function(x) x$outcome_name, character(1)))) {
  if (is.factor(df[[out_nm]])) {
    levs <- levels(df[[out_nm]])
    if (length(levs) == 2) {
      df[[out_nm]] <- as.integer(df[[out_nm]] == levs[2])
    }
  }
}

# gender
if ("Gender" %in% names(df)) {
  df$Gender <- prepare_gender(df$Gender)
}

# binary predictors
binary_vars <- c(
  "Employed",
  "Married_bin",
  "Living_conditions",
  "Family_History_mood",
  "Psychiatric_Comorbidity",
  "Substance_Use_Problem",
  "Psychotic_symptoms_bin",
  "Suicidal_ideation_bin"
)

for (v in intersect(binary_vars, names(df))) {
  df[[v]] <- prepare_binary_factor(df[[v]])
}

# ordered 3-category predictors
ordered_3cat_vars <- c(
  "Education_3cat",
  "Anxiety_3cat",
  "Depressive_symptoms_3cat",
  "Social_Withdrawal_3cat"
)

for (v in intersect(ordered_3cat_vars, names(df))) {
  df[[v]] <- prepare_ordered_3cat(df[[v]])
}

# ------------------------------------------------------------
# MODEL FITTER
# ------------------------------------------------------------
fit_multivariable_model <- function(data, outcome_var, predictors, model_label) {
  
  predictors_in <- predictors[predictors %in% names(data)]
  missing_pred  <- setdiff(predictors, predictors_in)
  
  if (length(missing_pred) > 0) {
    cat(model_label, " | Missing predictors dropped: ",
        paste(missing_pred, collapse = ", "), "\n", sep = "")
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
# RUN
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

