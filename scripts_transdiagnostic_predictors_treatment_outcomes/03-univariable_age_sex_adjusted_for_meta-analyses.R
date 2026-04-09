# ============================================================
# Project: Transdiagnostic predictors of treatment outcomes
# Script: Univariable age- and sex-adjusted models for meta-analyses
#
# Purpose:
# - fit cohort-specific univariable logistic regression models adjusted
#   for age and sex
# - return harmonised effect estimates that can be used as inputs for
#   downstream meta-analyses across cohorts
#
# General model:
#   outcome ~ predictor + AGE + SEX
#
# Additional note:
# - AGE and SEX themselves are estimated from minimally adjusted models:
#     outcome ~ AGE + SEX
#   so that AGE is adjusted for SEX and SEX is adjusted for AGE
#
# This same script can be adapted for the three phenotypes considered:
# - remission
# - TR1 (treatment-based phenotype)
# - TR2 (symptom-based phenotype)
#
# Notes:
# - the script assumes that each cohort has already undergone the main
#   preprocessing workflow
# - any cohort-specific cleaning, exclusions, or handling of outliers
#   should be completed before this step
# - the goal here is to produce comparable cohort-level effect estimates,
#   not to enforce identical raw coding across studies
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(broom)
  library(purrr)
})

options(stringsAsFactors = FALSE)

# ============================================================
# Input cohorts
# ============================================================
# Replace these objects with the three cohort-level analytic datasets.
# Each dataset should already contain:
# - one row per participant
# - harmonised predictors
# - binary outcomes coded as 0/1
# - age and sex variables or variables that can be mapped to them
# ============================================================

cohort_list <- list(
  COHORT_1 = cohort1_df,
  COHORT_2 = cohort2_df,
  COHORT_3 = cohort3_df
)

# ============================================================
# Outcome harmonisation
# ============================================================
# This object maps the phenotype of interest to the outcome variable name
# used in each cohort.
#
# Update only this section when moving from remission to TR1 or TR2.
# The rest of the script can remain unchanged.
# ============================================================

phenotype_map <- list(
  remission = c(
    COHORT_1 = "REMISSION_OUTCOME_NAME_1",
    COHORT_2 = "REMISSION_OUTCOME_NAME_2",
    COHORT_3 = "REMISSION_OUTCOME_NAME_3"
  ),
  TR1 = c(
    COHORT_1 = "TR1_OUTCOME_NAME_1",
    COHORT_2 = "TR1_OUTCOME_NAME_2",
    COHORT_3 = "TR1_OUTCOME_NAME_3"
  ),
  TR2 = c(
    COHORT_1 = "TR2_OUTCOME_NAME_1",
    COHORT_2 = "TR2_OUTCOME_NAME_2",
    COHORT_3 = "TR2_OUTCOME_NAME_3"
  )
)

# Select which phenotype to run
phenotype_to_run <- "remission"

# ============================================================
# Predictor lists by cohort
# ============================================================
# Predictors do not need to be identical across cohorts, but they should
# follow the same broad harmonisation logic whenever possible.
#
# These vectors should include the candidate predictors intended for the
# cross-cohort meta-analytic step.
# ============================================================

predictor_map <- list(
  COHORT_1 = c(
    "AGE",
    "SEX",
    "EMPLOYMENT",
    "EDUCATION_3CAT_ORD",
    "SEVERITY_3CAT_ORD",
    "FUNCTIONING",
    "PRS_1_resid_z"
  ),
  COHORT_2 = c(
    "AGE",
    "SEX",
    "EMPLOYMENT",
    "EDUCATION_3CAT_ORD",
    "FUNCTIONING",
    "ANXIETY_3CAT_ORD",
    "PRS_1_resid_z"
  ),
  COHORT_3 = c(
    "AGE",
    "SEX",
    "EMPLOYMENT",
    "FUNCTIONING",
    "SEVERITY_3CAT_ORD",
    "PRS_1_resid_z"
  )
)

# ============================================================
# Optional variable classification
# ============================================================
# These lists are useful when a variable should be handled in a specific
# way regardless of its storage mode in the dataset.
# ============================================================

force_binary <- c(
  "SEX",
  "EMPLOYMENT"
)

force_ordinal <- c(
  "EDUCATION_3CAT_ORD",
  "SEVERITY_3CAT_ORD",
  "ANXIETY_3CAT_ORD"
)

# ============================================================
# Shared helpers
# ============================================================

gelman_2sd <- function(x) {
  x <- as.numeric(x)
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  
  if (!is.finite(s) || s == 0) {
    return(rep(NA_real_, length(x)))
  }
  
  (x - m) / (2 * s)
}

to_binary_01 <- function(x) {
  if (is.logical(x)) {
    return(ifelse(is.na(x), NA_real_, ifelse(x, 1, 0)))
  }
  
  if (is.numeric(x)) {
    return(as.numeric(x))
  }
  
  x0 <- x
  if (is.factor(x0)) x0 <- as.character(x0)
  xl <- tolower(x0)
  
  yes_vals <- c("yes", "y", "true", "t", "1", "present", "positive")
  no_vals  <- c("no", "n", "false", "f", "0", "absent", "negative")
  
  out <- rep(NA_real_, length(x0))
  out[xl %in% yes_vals] <- 1
  out[xl %in% no_vals]  <- 0
  
  if (all(is.na(out))) {
    u <- sort(unique(x0[!is.na(x0)]))
    if (length(u) == 2) {
      out[x0 == u[1]] <- 0
      out[x0 == u[2]] <- 1
    }
  }
  
  out
}

infer_type <- function(data, varname) {
  if (varname %in% force_binary) return("binary")
  if (varname %in% force_ordinal) return("ordinal")
  if (grepl("_resid_z$", varname)) return("continuous")
  
  x <- data[[varname]]
  
  if (is.ordered(x)) return("ordinal")
  if (is.factor(x)) return("categorical")
  if (is.numeric(x) || is.integer(x)) return("continuous")
  
  "categorical"
}

level_label_from_term_ord <- function(pred, term, ord_spec) {
  if (!pred %in% names(ord_spec)) return(NA_character_)
  
  levs <- ord_spec[[pred]]
  ref  <- levs[1]
  
  k <- suppressWarnings(as.integer(sub(paste0("^", pred), "", term)))
  if (!is.na(k) && k >= 2 && k <= length(levs)) {
    return(paste0(levs[k], " vs ", ref))
  }
  
  lvl <- sub(paste0("^", pred), "", term)
  if (nzchar(lvl)) return(paste0(lvl, " vs ", ref))
  
  NA_character_
}

level_label_from_term_cat <- function(data, pred, term) {
  if (!pred %in% names(data)) return(NA_character_)
  
  f <- data[[pred]]
  if (!is.factor(f)) return(NA_character_)
  
  levs <- levels(f)
  if (length(levs) < 2) return(NA_character_)
  
  ref <- levs[1]
  lvl <- sub(paste0("^", pred), "", term)
  
  if (!nzchar(lvl)) return(NA_character_)
  paste0(lvl, " vs ", ref)
}

fit_one_adjusted <- function(data, outcome, predictor) {
  f <- as.formula(paste0(outcome, " ~ ", predictor, " + AGE + SEX"))
  vars_needed <- all.vars(f)
  
  if (!all(vars_needed %in% names(data))) return(NULL)
  
  d <- data[, vars_needed, drop = FALSE]
  idx <- complete.cases(d)
  
  if (sum(idx) < 30) return(NULL)
  
  mod <- try(glm(f, data = data[idx, , drop = FALSE], family = binomial()), silent = TRUE)
  if (inherits(mod, "try-error")) return(NULL)
  
  tt <- broom::tidy(mod) %>%
    dplyr::filter(term != "(Intercept)", term != "AGE", term != "SEX") %>%
    dplyr::mutate(
      outcome   = outcome,
      predictor = predictor,
      n         = nobs(mod),
      logOR     = estimate,
      SE        = std.error,
      OR        = exp(logOR),
      CI_lower  = exp(logOR - 1.96 * SE),
      CI_upper  = exp(logOR + 1.96 * SE),
      p_value   = p.value
    ) %>%
    dplyr::select(outcome, predictor, term, n, logOR, SE, OR, CI_lower, CI_upper, p_value)
  
  if (nrow(tt) == 0) return(NULL)
  tt
}

fit_age_sex_minadj <- function(data, outcome, sex_contrast_label = "Male vs Female") {
  f <- as.formula(paste0(outcome, " ~ AGE + SEX"))
  vars_needed <- all.vars(f)
  
  if (!all(vars_needed %in% names(data))) return(NULL)
  
  d <- data[, vars_needed, drop = FALSE]
  idx <- complete.cases(d)
  
  if (sum(idx) < 30) return(NULL)
  
  mod <- try(glm(f, data = data[idx, , drop = FALSE], family = binomial()), silent = TRUE)
  if (inherits(mod, "try-error")) return(NULL)
  
  tt <- broom::tidy(mod) %>%
    dplyr::filter(term %in% c("AGE", "SEX")) %>%
    dplyr::mutate(
      outcome   = outcome,
      predictor = term,
      n         = nobs(mod),
      logOR     = estimate,
      SE        = std.error,
      OR        = exp(logOR),
      CI_lower  = exp(logOR - 1.96 * SE),
      CI_upper  = exp(logOR + 1.96 * SE),
      p_value   = p.value,
      contrast_nice = dplyr::if_else(term == "SEX", sex_contrast_label, NA_character_)
    ) %>%
    dplyr::select(outcome, predictor, term, contrast_nice, n, logOR, SE, OR, CI_lower, CI_upper, p_value)
  
  if (nrow(tt) == 0) return(NULL)
  tt
}

prepare_cohort_for_models <- function(data, predictors, outcome_name) {
  df <- data
  
  # Standard age/sex variables expected by the downstream models.
  # Adapt here only if a cohort uses different names.
  if (!"AGE" %in% names(df)) stop("AGE variable not found in cohort dataset.")
  if (!"SEX" %in% names(df)) stop("SEX variable not found in cohort dataset.")
  
  # Coerce SEX to 0/1 if needed
  df$SEX <- to_binary_01(df$SEX)
  
  predictors_in <- predictors[predictors %in% names(df)]
  
  # Identify variable classes
  ord_pred <- predictors_in[predictors_in %in% force_ordinal]
  bin_pred <- predictors_in[predictors_in %in% force_binary]
  cont_pred <- predictors_in[!predictors_in %in% c(ord_pred, bin_pred)]
  
  # Continuous predictors are Gelman-rescaled
  cont_pred_num <- cont_pred[sapply(cont_pred, function(v) is.numeric(df[[v]]) || is.integer(df[[v]]))]
  for (v in cont_pred_num) {
    df[[paste0(v, "_g")]] <- gelman_2sd(df[[v]])
  }
  cont_pred_g <- paste0(cont_pred_num, "_g")
  
  # Binary predictors are coerced to 0/1 when useful
  for (v in bin_pred) {
    if (v %in% names(df)) {
      df[[v]] <- to_binary_01(df[[v]])
    }
  }
  
  # Ordinal predictors are explicitly ordered when appropriate
  ord_spec <- list()
  for (v in ord_pred) {
    if (!v %in% names(df)) next
    
    x <- df[[v]]
    if (is.factor(x) || is.character(x)) {
      levs <- unique(as.character(x[!is.na(x)]))
      if (all(c("Low", "Moderate", "High") %in% levs)) {
        ord_spec[[v]] <- c("Low", "Moderate", "High")
        df[[v]] <- ordered(df[[v]], levels = ord_spec[[v]])
        df[[v]] <- droplevels(df[[v]])
        if (nlevels(df[[v]]) >= 2) {
          contrasts(df[[v]]) <- contr.treatment(nlevels(df[[v]]), base = 1)
        }
      }
    }
  }
  
  cat_pred <- setdiff(
    predictors_in,
    c(cont_pred_num, ord_pred, bin_pred, "AGE", "SEX")
  )
  
  for (v in cat_pred) {
    if (v %in% names(df) && (is.factor(df[[v]]) || is.character(df[[v]]))) {
      df[[v]] <- droplevels(factor(df[[v]]))
      if (nlevels(df[[v]]) >= 2) {
        contrasts(df[[v]]) <- contr.treatment(nlevels(df[[v]]), base = 1)
      }
    }
  }
  
  list(
    data = df,
    all_predictors = unique(c(cont_pred_g, bin_pred, ord_pred, cat_pred)),
    ord_spec = ord_spec,
    cont_pred_g = cont_pred_g,
    bin_pred = bin_pred,
    ord_pred = ord_pred,
    cat_pred = cat_pred,
    outcome_name = outcome_name
  )
}

run_one_cohort <- function(data, cohort_name, outcome_name, predictors) {
  prep <- prepare_cohort_for_models(data, predictors, outcome_name)
  df <- prep$data
  
  adjusted_results <- bind_rows(
    purrr::map_dfr(prep$all_predictors, function(pred) {
      fit_one_adjusted(df, outcome_name, pred)
    }),
    fit_age_sex_minadj(df, outcome_name, sex_contrast_label = "Male vs Female")
  ) %>%
    dplyr::mutate(
      study = cohort_name,
      predictor_type = dplyr::case_when(
        predictor %in% c("AGE", "SEX") ~ "demographic (minimally adjusted: AGE + SEX)",
        predictor %in% prep$cont_pred_g ~ "continuous (Gelman 2SD)",
        predictor %in% prep$bin_pred ~ "binary",
        predictor %in% prep$ord_pred ~ "ordinal (treatment contrasts)",
        predictor %in% prep$cat_pred ~ "categorical (treatment contrasts)",
        TRUE ~ "other"
      ),
      contrast_nice = dplyr::case_when(
        predictor %in% prep$ord_pred ~ purrr::map2_chr(
          predictor, term, ~ level_label_from_term_ord(.x, .y, prep$ord_spec)
        ),
        predictor %in% prep$cat_pred ~ purrr::map2_chr(
          predictor, term, ~ level_label_from_term_cat(df, .x, .y)
        ),
        predictor == "SEX" ~ dplyr::coalesce(contrast_nice, "Male vs Female"),
        TRUE ~ NA_character_
      ),
      p_value_print = format.pval(p_value, digits = 8, eps = 1e-12)
    ) %>%
    dplyr::arrange(outcome, predictor, term)
  
  meta_ready <- adjusted_results %>%
    dplyr::mutate(
      effect_id = ifelse(
        is.na(contrast_nice),
        predictor,
        paste0(predictor, " | ", contrast_nice)
      ),
      yi  = logOR,
      sei = SE
    ) %>%
    dplyr::select(study, outcome, effect_id, predictor, term, contrast_nice, n, yi, sei, predictor_type)
  
  list(
    adjusted_univariate = adjusted_results,
    meta_ready = meta_ready
  )
}

# ============================================================
# Run selected phenotype across cohorts
# ============================================================

if (!phenotype_to_run %in% names(phenotype_map)) {
  stop("phenotype_to_run not found in phenotype_map.")
}

results_by_cohort <- list()

for (cohort_name in names(cohort_list)) {
  cohort_df <- cohort_list[[cohort_name]]
  
  if (!cohort_name %in% names(phenotype_map[[phenotype_to_run]])) {
    message("Skipping ", cohort_name, ": no outcome mapping provided for phenotype ", phenotype_to_run)
    next
  }
  
  outcome_name <- phenotype_map[[phenotype_to_run]][[cohort_name]]
  predictors <- predictor_map[[cohort_name]]
  
  if (is.null(predictors)) {
    message("Skipping ", cohort_name, ": no predictor list provided.")
    next
  }
  
  if (!outcome_name %in% names(cohort_df)) {
    message("Skipping ", cohort_name, ": outcome ", outcome_name, " not found in dataset.")
    next
  }
  
  results_by_cohort[[cohort_name]] <- run_one_cohort(
    data = cohort_df,
    cohort_name = cohort_name,
    outcome_name = outcome_name,
    predictors = predictors
  )
}

# ============================================================
# Final objects created
# ============================================================
# results_by_cohort[[<COHORT_NAME>]]$adjusted_univariate
# results_by_cohort[[<COHORT_NAME>]]$meta_ready
#
# Example:
# results_by_cohort[["COHORT_1"]]$meta_ready
# ============================================================