############################################################
# Project: Transdiagnostic predictors of treatment outcomes
# Script: Univariable analyses within each cohort
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

analysis_df <- cohort_df

# ============================================================
# Outcomes and predictors
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
# Minimum data requirements
# ============================================================
# These thresholds are pragmatic, user-defined safeguards rather than
# fixed statistical rules. They should be reviewed and adapted depending
# on the cohort size, outcome prevalence, predictor distribution, and
# modelling aims.
# ============================================================

min_complete_cases <- 5
min_unique_values  <- 2

# ============================================================
# Helper functions
# ============================================================

round3 <- function(x) {
  ifelse(is.na(x), NA, round(x, 3))
}

mean_sd_str <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  sprintf("%.3f ± %.3f", mean(x), sd(x))
}

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

clean_text <- function(s) {
  s <- as.character(s)
  s <- gsub("\t", " ", s)
  s <- gsub("\r", " ", s)
  s <- gsub("\n", " ", s)
  s <- gsub(" +", " ", s)
  trimws(s)
}

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

get_type <- function(df, v) {
  if (v %in% force_binary) return("categorical")
  if (v %in% force_ordinal) return("ordinal")
  if (grepl("_resid_z$", v)) return("continuous")
  
  if (is.ordered(df[[v]])) return("ordinal")
  if (is.factor(df[[v]])) return("categorical")
  if (is.numeric(df[[v]]) || is.integer(df[[v]])) return("continuous")
  
  "categorical"
}

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

wald_ci <- function(model) {
  b  <- coef(model)
  se <- sqrt(diag(vcov(model)))
  
  low  <- b - 1.96 * se
  high <- b + 1.96 * se
  
  out <- cbind(low, high)
  rownames(out) <- names(b)
  out
}

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
      select(all_of(c(outcome, v))) %>%
      filter(!is.na(.data[[outcome]]), !is.na(.data[[v]]))
    
    n_total <- nrow(tmp)
    
    # Skip predictors with too few complete observations
    if (n_total < min_complete_cases) next
    
    y <- tmp[[outcome]]
    x <- tmp[[v]]
    
    # Skip predictors with insufficient variation
    if (length(unique(x)) < min_unique_values) next
    
    # Descriptive summaries
    if (type == "continuous") {
      dist0 <- mean_sd_str(x[y == 0])
      dist1 <- mean_sd_str(x[y == 1])
    } else {
      dist0 <- cat_dist_str(x, y, 0)
      dist1 <- cat_dist_str(x, y, 1)
    }
    
    x_mod <- set_reference(x, type)
    tmp[[v]] <- x_mod
    
    full_formula <- as.formula(paste0(outcome, " ~ `", v, "`"))
    null_formula <- as.formula(paste0(outcome, " ~ 1"))
    
    model_full <- glm(full_formula, data = tmp, family = binomial())
    model_null <- glm(null_formula, data = tmp, family = binomial())
    
    lrt <- anova(model_null, model_full, test = "LRT")
    lrt_chi2 <- as.numeric(lrt$Deviance[2])
    
    sm <- summary(model_full)$coefficients
    ci <- wald_ci(model_full)
    
    ref_level <- if (type == "continuous") "" else clean_text(levels(tmp[[v]])[1])
    
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
# Run analyses
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
# Export
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
# Preview
# ============================================================

for (outcome_name in names(univariable_results)) {
  cat("\n\n==============================\n")
  cat("Outcome:", outcome_name, "\n")
  cat("==============================\n")
  print(utils::head(univariable_results[[outcome_name]], 10))
}
