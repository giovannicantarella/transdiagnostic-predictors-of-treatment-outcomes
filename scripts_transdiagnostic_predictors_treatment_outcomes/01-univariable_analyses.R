############################################################
# Project: Transdiagnostic predictors of treatment outcomes
# Script: Generic univariable analyses
#
# Purpose:
# - summarise candidate predictors by binary outcome status
# - run one univariable logistic regression per predictor
# - return effect sizes as odds ratios with 95% confidence intervals
# - provide a consistent table structure across cohorts and outcomes
#
# Notes:
# - this script is written as a generic template for cohort-level
#   analytic datasets that already underwent preprocessing
# - any cohort-specific data cleaning (e.g. handling of outliers,
#   harmonisation of labels, exclusions based on quality criteria)
#   should be completed before this step
# - the script assumes binary outcomes coded as 0/1
############################################################

library(dplyr)

options(stringsAsFactors = FALSE)

# ============================================================
# Input dataset
# ============================================================
# Replace this object with the analytic dataset for the cohort of interest.
# The dataset should already contain:
# - one row per participant
# - harmonised predictors
# - binary outcomes coded as 0/1
# ============================================================

analysis_df <- cohort_df

# ============================================================
# Outcomes and predictors
# ============================================================
# Replace these placeholders with the outcomes and predictors available
# in the current cohort. The same script can be reused across cohorts by
# updating only these vectors.
#
# Suggested use:
# - outcomes: remission and/or resistance phenotypes
# - predictors: sociodemographic, clinical, and genetic variables
# ============================================================

outcomes <- c("OUTCOME_1", "OUTCOME_2")

predictors <- c(
  "AGE",
  "SEX",
  "EMPLOYMENT",
  "EDUCATION_3CAT_ORD",
  "FUNCTIONING",
  "SEVERITY_3CAT_ORD",
  "SOCIAL_WITHDRAWAL_3CAT_ORD",
  "PRS_1_resid_z",
  "PRS_2_resid_z"
)

# ============================================================
# Optional predictor type overrides
# ============================================================
# In many datasets, variable classes are already informative
# (numeric, factor, ordered factor). Still, it is often useful to force
# some variables to be treated in a specific way regardless of storage.
#
# Examples:
# - binary variables stored as 0/1 numeric can be forced to categorical
# - grouped severity measures can be forced to ordinal
# ============================================================

force_binary <- c(
  "SEX",
  "EMPLOYMENT"
)

force_ordinal <- c(
  "EDUCATION_3CAT_ORD",
  "SEVERITY_3CAT_ORD",
  "SOCIAL_WITHDRAWAL_3CAT_ORD"
)

# ============================================================
# Helper functions
# ============================================================

round3 <- function(x) {
  ifelse(is.na(x), NA, round(x, 3))
}

# Summary string for continuous predictors
mean_sd_str <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  sprintf("%.3f ± %.3f", mean(x), sd(x))
}

# Summary string for categorical / ordinal predictors
cat_dist_str <- function(x, y, yval) {
  xi <- x[y == yval]
  tab <- table(xi, useNA = "no")
  denom <- sum(tab)
  
  if (denom == 0) return(NA_character_)
  
  parts <- sprintf(
    "%s: %d (%.1f%%)",
    names(tab),
    as.integer(tab),
    100 * as.numeric(tab) / denom
  )
  
  paste(parts, collapse = " | ")
}

# Clean text to avoid issues in exports
clean_text <- function(s) {
  s <- as.character(s)
  s <- gsub("\t", " ", s)
  s <- gsub("\r", " ", s)
  s <- gsub("\n", " ", s)
  s <- gsub(" +", " ", s)
  trimws(s)
}

# Make character columns safer for spreadsheet export
clean_for_export <- function(df, sep = ";") {
  out <- df
  chr_cols <- names(out)[vapply(out, is.character, logical(1))]
  
  for (cc in chr_cols) {
    x <- clean_text(out[[cc]])
    x <- gsub(sep, " | ", x, fixed = TRUE)
    out[[cc]] <- x
  }
  
  out
}

# Infer predictor type using object class plus optional overrides
get_type <- function(df, v) {
  if (v %in% force_binary) return("categorical")
  if (v %in% force_ordinal) return("ordinal")
  if (grepl("_resid_z$", v)) return("continuous")
  
  if (is.ordered(df[[v]])) return("ordinal")
  if (is.factor(df[[v]])) return("categorical")
  if (is.numeric(df[[v]]) || is.integer(df[[v]])) return("continuous")
  
  "categorical"
}

# Set modelling reference level
# For ordinal variables, the default behaviour here is to use the first
# level as reference. If a different reference is preferred, change this
# function rather than editing the modelling code later.
set_reference <- function(x, type) {
  
  if (type == "continuous") return(x)
  
  if (is.ordered(x)) {
    levs <- levels(x)
    x <- factor(as.character(x), levels = levs, ordered = FALSE)
  } else if (!is.factor(x)) {
    x <- factor(x)
  }
  
  if (nlevels(x) < 2) return(x)
  
  stats::relevel(x, ref = levels(x)[1])
}

# Wald confidence intervals on the log-odds scale
wald_ci <- function(model) {
  b  <- coef(model)
  se <- sqrt(diag(vcov(model)))
  
  low  <- b - 1.96 * se
  high <- b + 1.96 * se
  
  out <- cbind(low, high)
  rownames(out) <- names(b)
  out
}

# Extract category label from model term
extract_level <- function(varname, term) {
  t <- term
  t <- gsub("^`", "", t)
  t <- gsub("`$", "", t)
  t <- sub(paste0("^", varname), "", t)
  t <- sub(paste0("^`", varname, "`"), "", t)
  clean_text(t)
}

# ============================================================
# Main function
# ============================================================
# For each predictor:
# 1. keep complete cases on outcome and predictor
# 2. describe predictor separately for outcome = 0 and outcome = 1
# 3. fit a univariable logistic model
# 4. extract coefficient-level ORs and 95% CIs
# 5. store the global likelihood ratio statistic for the predictor
#
# Output structure:
# - one row for each continuous predictor
# - one row per non-reference level for categorical / ordinal predictors
# ============================================================

build_univariable_table <- function(df, outcome, predictors) {
  
  if (!outcome %in% names(df)) {
    stop("Outcome not found: ", outcome)
  }
  
  y_all <- df[[outcome]]
  
  if (!all(na.omit(unique(y_all)) %in% c(0, 1))) {
    stop("Outcome ", outcome, " is not coded as 0/1.")
  }
  
  n0 <- sum(y_all == 0, na.rm = TRUE)
  n1 <- sum(y_all == 1, na.rm = TRUE)
  
  out_col0 <- paste0(outcome, "=0 (N=", n0, ")")
  out_col1 <- paste0(outcome, "=1 (N=", n1, ")")
  
  results_list <- list()
  
  for (v in predictors) {
    
    if (!v %in% names(df)) {
      message("Skipping missing predictor: ", v)
      next
    }
    
    type <- get_type(df, v)
    
    tmp <- df %>%
      dplyr::select(dplyr::all_of(c(outcome, v))) %>%
      dplyr::filter(!is.na(.data[[outcome]]), !is.na(.data[[v]]))
    
    n_total <- nrow(tmp)
    if (n_total < 5) next
    
    y <- tmp[[outcome]]
    x <- tmp[[v]]
    
    if (length(unique(x)) < 2) next
    
    # Descriptive summaries
    if (type == "continuous") {
      dist0 <- mean_sd_str(x[y == 0])
      dist1 <- mean_sd_str(x[y == 1])
    } else {
      dist0 <- cat_dist_str(x, y, 0)
      dist1 <- cat_dist_str(x, y, 1)
    }
    
    # Prepare variable for modelling
    x_mod <- set_reference(x, type)
    tmp[[v]] <- x_mod
    
    # Fit models for LRT
    full_formula  <- as.formula(paste0(outcome, " ~ `", v, "`"))
    null_formula  <- as.formula(paste0(outcome, " ~ 1"))
    
    model_full <- glm(full_formula, data = tmp, family = binomial())
    model_null <- glm(null_formula, data = tmp, family = binomial())
    
    lrt <- anova(model_null, model_full, test = "LRT")
    lrt_chi2 <- as.numeric(lrt$Deviance[2])
    
    sm <- summary(model_full)$coefficients
    ci <- wald_ci(model_full)
    
    ref_level <- if (type == "continuous") "" else clean_text(levels(tmp[[v]])[1])
    
    # Continuous predictor: one coefficient only
    if (type == "continuous") {
      
      beta    <- as.numeric(sm[2, "Estimate"])
      pval    <- as.numeric(sm[2, "Pr(>|z|)"])
      or_val  <- exp(beta)
      ci_low  <- exp(ci[2, 1])
      ci_high <- exp(ci[2, 2])
      
      row <- data.frame(
        Variable = clean_text(v),
        Type = type,
        `Reference Level` = ref_level,
        N_total = n_total,
        Statistic = beta,
        LRT_Chi2 = lrt_chi2,
        OR = or_val,
        CI_lower = ci_low,
        CI_upper = ci_high,
        p_value = pval,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      
      row[[out_col0]] <- clean_text(dist0)
      row[[out_col1]] <- clean_text(dist1)
      
      row <- row[, c(
        "Variable", out_col0, out_col1, "Type", "Reference Level", "N_total",
        "Statistic", "LRT_Chi2", "OR", "CI_lower", "CI_upper", "p_value"
      )]
      
      results_list[[length(results_list) + 1]] <- row
      
    } else {
      
      # Categorical / ordinal predictor: one row for each non-reference level
      rn <- rownames(sm)
      rn <- rn[rn != "(Intercept)"]
      
      for (term in rn) {
        
        beta    <- as.numeric(sm[term, "Estimate"])
        pval    <- as.numeric(sm[term, "Pr(>|z|)"])
        or_val  <- exp(beta)
        ci_low  <- exp(ci[term, 1])
        ci_high <- exp(ci[term, 2])
        
        level <- extract_level(v, term)
        var_label <- clean_text(paste0(v, ": ", level))
        
        row <- data.frame(
          Variable = var_label,
          Type = type,
          `Reference Level` = ref_level,
          N_total = n_total,
          Statistic = beta,
          LRT_Chi2 = lrt_chi2,
          OR = or_val,
          CI_lower = ci_low,
          CI_upper = ci_high,
          p_value = pval,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )
        
        row[[out_col0]] <- clean_text(dist0)
        row[[out_col1]] <- clean_text(dist1)
        
        row <- row[, c(
          "Variable", out_col0, out_col1, "Type", "Reference Level", "N_total",
          "Statistic", "LRT_Chi2", "OR", "CI_lower", "CI_upper", "p_value"
        )]
        
        results_list[[length(results_list) + 1]] <- row
      }
    }
  }
  
  out <- bind_rows(results_list)
  
  num_cols <- c("N_total", "Statistic", "LRT_Chi2", "OR", "CI_lower", "CI_upper", "p_value")
  
  for (cc in intersect(num_cols, names(out))) {
    out[[cc]] <- round3(as.numeric(out[[cc]]))
  }
  
  clean_for_export(out, sep = ";")
}

# ============================================================
# Run univariable analyses for all outcomes
# ============================================================
# The result is a named list of tables, one per outcome.
# ============================================================

univariable_results <- lapply(outcomes, function(outcome_name) {
  build_univariable_table(
    df = analysis_df,
    outcome = outcome_name,
    predictors = predictors
  )
})

names(univariable_results) <- outcomes

# ============================================================
# Optional export
# ============================================================
# One semicolon-separated csv per outcome.
# This is convenient for spreadsheet software and for repository sharing.
# ============================================================

for (outcome_name in names(univariable_results)) {
  write.table(
    univariable_results[[outcome_name]],
    file = paste0("univariable_", outcome_name, ".csv"),
    sep = ";",
    quote = TRUE,
    row.names = FALSE,
    na = ""
  )
}

# ============================================================
# Optional console preview
# ============================================================

for (outcome_name in names(univariable_results)) {
  cat("\n\n==============================\n")
  cat("Outcome:", outcome_name, "\n")
  cat("==============================\n")
  print(utils::head(univariable_results[[outcome_name]], 10))
}