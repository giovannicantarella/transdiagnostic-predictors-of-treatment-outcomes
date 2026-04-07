# =========================================================
# UNIVARIATE TABLES (STAR*D external validation sample)
# - Outcomes: Resistance, remission
# - Predictors included here were restricted to variables showing
#   meta-analytic associations at p < 0.10 and a consistent direction
#   of effect in at least 2 out of 3 discovery cohorts
# - Predictors include clinical, sociodemographic, and PRS residualized z-scores
# - Output formatted as Excel-safe ";" CSV strings
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

# -------------------------
# Outcomes and predictors
# -------------------------
outcomes <- c("Resistance", "remission")

predictors_STAR <- c(
  # Sociodemographic variables
  "Age", "Gender", "Ethnicity",
  "Employed", "Married_bin", "Education_3cat", "Living_conditions",
  
  # Clinical history
  "Age_at_onset", "N_previous_episodes_log",
  "Family_History_mood",
  
  # Comorbidity / risk
  "Psychiatric_Comorbidity",
  "Substance_Use_Problem",
  
  # Baseline severity / symptom domains
  "HRSD_total",
  "Anxiety_3cat",
  "Depressive_symptoms_3cat",
  "Psychotic_symptoms_bin",
  "Suicidal_ideation_bin",
  
  # Functioning / social withdrawal
  "WSAS_baseline",
  "Social_Withdrawal_3cat",
  
  # Polygenic risk scores (residualized + z-standardized)
  "ADHD05_PRS_resid_z",
  "BODY04_PRS_resid_z",
  "BP_PRS_resid_z",
  "DEPR07_PRS_resid_z",
  "EDUC03_PRS_resid_z",
  "SCZ_PRS_resid_z"
)

# Variables that must be treated as categorical even if stored numerically
force_binary <- c(
  "Gender", "Ethnicity",
  "Employed", "Married_bin", "Living_conditions",
  "Family_History_mood",
  "Psychiatric_Comorbidity", "Substance_Use_Problem",
  "Psychotic_symptoms_bin", "Suicidal_ideation_bin"
)

# Optional manual reference levels
reference_levels <- list(
  Family_History_mood = "No",
  Psychiatric_Comorbidity = "No",
  Substance_Use_Problem = "No",
  Psychotic_symptoms_bin = "No",
  Suicidal_ideation_bin = "No"
)

# -------------------------
# Helpers
# -------------------------
round3 <- function(x) ifelse(is.na(x), NA, round(x, 3))

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

clean_for_excel <- function(df, sep = ";") {
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
  if (is.numeric(df[[v]]) || is.integer(df[[v]])) return("continuous")
  return("categorical")
}

set_reference <- function(x, varname, type) {
  if (!type %in% c("categorical", "ordinal")) return(x)
  
  if (is.ordered(x)) {
    x <- factor(as.character(x), levels = levels(x), ordered = FALSE)
  } else if (!is.factor(x)) {
    x <- factor(x)
  }
  
  lv <- levels(x)
  
  if (varname %in% names(reference_levels)) {
    ref <- reference_levels[[varname]]
    if (!is.na(ref) && ref %in% lv) {
      return(stats::relevel(x, ref = ref))
    }
  }
  
  if (type == "ordinal") {
    return(stats::relevel(x, ref = tail(lv, 1)))
  }
  
  stats::relevel(x, ref = lv[1])
}

wald_ci <- function(model) {
  b <- coef(model)
  se <- sqrt(diag(vcov(model)))
  low <- b - 1.96 * se
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

# -------------------------
# Build one univariate table
# -------------------------
build_uni_table <- function(df, outcome, predictors) {
  
  if (!outcome %in% names(df)) stop("Outcome not found: ", outcome)
  
  y_all <- df[[outcome]]
  if (!all(na.omit(unique(y_all)) %in% c(0, 1))) {
    stop("Outcome ", outcome, " is not coded as 0/1.")
  }
  
  N0 <- sum(y_all == 0, na.rm = TRUE)
  N1 <- sum(y_all == 1, na.rm = TRUE)
  
  out_col0 <- paste0(outcome, "=0 (N=", N0, ")")
  out_col1 <- paste0(outcome, "=1 (N=", N1, ")")
  
  n_levels_ok <- function(x) {
    x <- x[!is.na(x)]
    length(unique(as.character(x))) >= 2
  }
  
  res_list <- list()
  
  for (v in predictors) {
    
    if (!v %in% names(df)) stop("Predictor not found in STAR: ", v)
    
    type <- get_type(df, v)
    
    tmp <- df %>%
      dplyr::select(all_of(outcome), all_of(v)) %>%
      dplyr::filter(!is.na(.data[[outcome]]), !is.na(.data[[v]]))
    
    n_total <- nrow(tmp)
    if (n_total < 5) next
    
    y <- tmp[[outcome]]
    x <- tmp[[v]]
    
    if (type %in% c("categorical", "ordinal")) {
      if (!n_levels_ok(x)) next
    } else if (type == "continuous") {
      xx <- x[!is.na(x)]
      if (length(xx) < 2) next
      if (isTRUE(all.equal(stats::sd(xx), 0))) next
    }
    
    if (type == "continuous") {
      dist0 <- mean_sd_str(x[y == 0])
      dist1 <- mean_sd_str(x[y == 1])
    } else {
      dist0 <- cat_dist_str(x, y, 0)
      dist1 <- cat_dist_str(x, y, 1)
    }
    
    x_mod <- set_reference(x, v, type)
    
    if (type %in% c("categorical", "ordinal")) {
      if (length(unique(as.character(x_mod[!is.na(x_mod)]))) < 2) next
    }
    
    tmp[[v]] <- x_mod
    
    f1 <- as.formula(paste0(outcome, " ~ `", v, "`"))
    m1 <- glm(f1, data = tmp, family = binomial())
    m0 <- glm(as.formula(paste0(outcome, " ~ 1")), data = tmp, family = binomial())
    
    lrt <- anova(m0, m1, test = "LRT")
    LRT_Chi2 <- as.numeric(lrt$Deviance[2])
    
    sm <- summary(m1)$coefficients
    ci <- wald_ci(m1)
    
    ref_level <- if (type == "continuous") "" else clean_text(levels(tmp[[v]])[1])
    
    if (type == "continuous") {
      
      beta <- as.numeric(sm[2, "Estimate"])
      pval <- as.numeric(sm[2, "Pr(>|z|)"])
      OR <- exp(beta)
      CI_low <- exp(ci[2, 1])
      CI_high <- exp(ci[2, 2])
      
      row <- data.frame(
        Variable = clean_text(v),
        Type = type,
        `Reference Level` = ref_level,
        N_total = n_total,
        Statistic = beta,
        LRT_Chi2 = LRT_Chi2,
        OR = OR,
        CI_lower = CI_low,
        CI_upper = CI_high,
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
      
      res_list[[length(res_list) + 1]] <- row
      
    } else {
      
      rn <- rownames(sm)
      rn <- rn[rn != "(Intercept)"]
      
      for (term in rn) {
        
        beta <- as.numeric(sm[term, "Estimate"])
        pval <- as.numeric(sm[term, "Pr(>|z|)"])
        OR <- exp(beta)
        CI_low <- exp(ci[term, 1])
        CI_high <- exp(ci[term, 2])
        
        level <- extract_level(v, term)
        var_label <- clean_text(paste0(v, ": ", level))
        
        row <- data.frame(
          Variable = var_label,
          Type = type,
          `Reference Level` = ref_level,
          N_total = n_total,
          Statistic = beta,
          LRT_Chi2 = LRT_Chi2,
          OR = OR,
          CI_lower = CI_low,
          CI_upper = CI_high,
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
        
        res_list[[length(res_list) + 1]] <- row
      }
    }
  }
  
  out <- bind_rows(res_list)
  
  num_cols <- c("N_total", "Statistic", "LRT_Chi2", "OR", "CI_lower", "CI_upper", "p_value")
  for (cc in intersect(num_cols, names(out))) {
    out[[cc]] <- round3(as.numeric(out[[cc]]))
  }
  
  clean_for_excel(out, sep = ";")
}

# -------------------------
# Export as ";" CSV string
# -------------------------
to_csv_string <- function(df, sep = ";") {
  df2 <- clean_for_excel(df, sep = sep)
  con <- textConnection("out", "w", local = TRUE)
  write.table(df2, file = con, sep = sep, quote = TRUE, row.names = FALSE, na = "")
  close(con)
  paste(out, collapse = "\n")
}

# -------------------------
# Build tables
# -------------------------
tab_Resistance <- build_uni_table(STAR, "Resistance", predictors_STAR)
tab_remission  <- build_uni_table(STAR, "remission", predictors_STAR)

cat("\n\n===== TABLE: Resistance (CSV ; ) =====\n")
cat(to_csv_string(tab_Resistance, sep = ";"))

cat("\n\n===== TABLE: remission (CSV ; ) =====\n")
cat(to_csv_string(tab_remission, sep = ";"))