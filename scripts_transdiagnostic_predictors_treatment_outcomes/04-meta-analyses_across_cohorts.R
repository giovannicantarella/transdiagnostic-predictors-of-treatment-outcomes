# ============================================================
# Project: Transdiagnostic predictors of treatment outcomes
# Script: Meta-analyses across cohorts
#
# Purpose:
# - combine age- and sex-adjusted cohort-specific effect estimates
# - harmonise corresponding predictors across cohorts
# - run random-effects meta-analyses
# - summarise heterogeneity and leave-one-out results
# - optionally export forest plots
#
# This same script can be adapted for the three phenotypes considered:
# - remission
# - TR1 (treatment-based phenotype)
# - TR2 (symptom-based phenotype)
#
# Inputs:
# - meta-ready effect estimates produced by:
#     univariable_age_sex_adjusted_for_meta-analyses.R
#
# General principles:
# - matching of cohort-specific effects is explicit rather than implicit
# - direction harmonisation is handled by optional sign flipping
# - the same workflow can be reused across phenotypes by changing only:
#   1) phenotype_to_run
#   2) outcome mapping
#   3) correspondence map
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(metafor)
  library(tibble)
})

# ============================================================
# Required input object
# ============================================================
# This script expects the object results_by_cohort to already exist,
# created by the upstream age- and sex-adjusted script.
# ============================================================

if (!exists("results_by_cohort")) {
  stop("Missing required object 'results_by_cohort'. Run the upstream univariable script first.")
}

# ============================================================
# Phenotype selection
# ============================================================
# Change only this value to move between remission, TR1 and TR2.
# ============================================================

phenotype_to_run <- "remission"

# ============================================================
# Helpers
# ============================================================

`%||%` <- function(a, b) if (!is.null(a)) a else b

ek <- function(predictor, contrast = NA_character_) {
  if (is.na(contrast) || contrast == "" || contrast == "NA") {
    predictor
  } else {
    paste0(predictor, " | ", contrast)
  }
}

harmonize_meta_ready <- function(df, study_name) {
  df %>%
    mutate(
      study = study_name,
      outcome = as.character(outcome),
      predictor = as.character(predictor),
      contrast_nice = ifelse(
        is.na(contrast_nice) | contrast_nice == "NA",
        NA_character_,
        as.character(contrast_nice)
      ),
      effect_id = ifelse(
        is.na(contrast_nice),
        predictor,
        paste0(predictor, " | ", contrast_nice)
      ),
      yi = as.numeric(yi),
      sei = as.numeric(sei),
      n = suppressWarnings(as.integer(n))
    ) %>%
    select(study, outcome, predictor, contrast_nice, effect_id, n, yi, sei, predictor_type)
}

pick_effect <- function(data, study_name, outcome_name, key) {
  if (is.null(key) || is.na(key) || key == "") return(NULL)
  
  data %>%
    filter(
      study == study_name,
      outcome == outcome_name,
      effect_id == key
    )
}

apply_flip <- function(df, flip_flag = FALSE) {
  if (is.null(df) || nrow(df) == 0) return(df)
  if (isTRUE(flip_flag)) {
    df <- df %>% mutate(yi = -yi)
  }
  df
}

build_effect_from_map <- function(spec, effect_label, pooled_data, outcome_map) {
  effect_rows <- purrr::imap_dfr(outcome_map, function(outcome_name, study_name) {
    d <- apply_flip(
      pick_effect(pooled_data, study_name, outcome_name, spec[[study_name]]),
      isTRUE(spec[[paste0("flip_", study_name)]])
    )
    
    if (is.null(d) || nrow(d) == 0) return(tibble())
    
    d %>%
      mutate(
        flipped = isTRUE(spec[[paste0("flip_", study_name)]])
      )
  })
  
  effect_rows %>%
    mutate(
      meta_effect = effect_label,
      note = spec$note %||% NA_character_
    ) %>%
    select(
      meta_effect, study, outcome, predictor, contrast_nice, effect_id,
      n, yi, sei, predictor_type, flipped, note
    )
}

run_meta_reml <- function(dat) {
  dat <- dat %>% filter(!is.na(yi), !is.na(sei))
  
  if (nrow(dat) == 0) return(NULL)
  
  per_study <- dat %>%
    mutate(
      OR = exp(yi),
      CI_lower = exp(yi - 1.96 * sei),
      CI_upper = exp(yi + 1.96 * sei),
      z = yi / sei,
      p = 2 * pnorm(-abs(z))
    )
  
  if (nrow(dat) == 1) {
    single <- per_study %>%
      transmute(
        estimate = yi,
        se = sei,
        ci.lb = yi - 1.96 * sei,
        ci.ub = yi + 1.96 * sei,
        OR = OR,
        OR_ci.lb = CI_lower,
        OR_ci.ub = CI_upper,
        pval = p
      )
    
    return(list(
      k = 1,
      per_study = per_study,
      model = NULL,
      pooled = single,
      heterogeneity = tibble(
        tau2 = NA_real_,
        tau = NA_real_,
        I2 = NA_real_,
        H2 = NA_real_,
        Q = NA_real_,
        Q_df = NA_real_,
        Q_pval = NA_real_
      ),
      leave1out = NULL
    ))
  }
  
  model <- metafor::rma(yi = yi, sei = sei, data = dat, method = "REML")
  sm <- summary(model)
  
  pooled <- tibble(
    estimate = as.numeric(sm$b),
    se = sm$se,
    z = sm$zval,
    pval = sm$pval,
    ci.lb = sm$ci.lb,
    ci.ub = sm$ci.ub,
    OR = exp(as.numeric(sm$b)),
    OR_ci.lb = exp(sm$ci.lb),
    OR_ci.ub = exp(sm$ci.ub)
  )
  
  heterogeneity <- tibble(
    tau2 = sm$tau2,
    tau = sqrt(sm$tau2),
    I2 = sm$I2,
    H2 = sm$H2,
    Q = sm$QE,
    Q_df = sm$k - 1,
    Q_pval = sm$QEp
  )
  
  loo <- if (nrow(dat) >= 3) {
    as.data.frame(leave1out(model)) %>%
      tibble::rownames_to_column("left_out_study")
  } else {
    NULL
  }
  
  list(
    k = nrow(dat),
    per_study = per_study,
    model = model,
    pooled = pooled,
    heterogeneity = heterogeneity,
    leave1out = loo
  )
}

print_meta_result <- function(res, label) {
  cat("\n", strrep("=", 88), "\n", sep = "")
  cat("META-EFFECT:", label, "\n")
  cat(strrep("=", 88), "\n", sep = "")
  
  if (is.null(res)) {
    cat("No analysable data.\n")
    return(invisible(NULL))
  }
  
  cat("Per-study estimates:\n")
  print(
    res$per_study %>%
      select(study, outcome, predictor, contrast_nice, n, yi, sei, OR, CI_lower, CI_upper, p, flipped),
    n = Inf,
    width = Inf
  )
  
  if (res$k == 1) {
    cat("\nOnly one study available: pooled meta-analysis not run.\n")
    print(res$pooled, n = Inf, width = Inf)
    return(invisible(NULL))
  }
  
  cat("\nPooled random-effects model (REML):\n")
  print(res$pooled, n = Inf, width = Inf)
  
  cat("\nHeterogeneity:\n")
  print(res$heterogeneity, n = Inf, width = Inf)
  
  if (!is.null(res$leave1out)) {
    cat("\nLeave-one-out analysis:\n")
    print(res$leave1out, n = Inf, width = Inf)
  }
  
  invisible(NULL)
}

# ============================================================
# Optional forest plotting helpers
# ============================================================

.calc_or_xlim <- function(yi, sei, pad = 0.15, max_ratio = 80) {
  lo <- exp(min(yi - 1.96 * sei, na.rm = TRUE))
  hi <- exp(max(yi + 1.96 * sei, na.rm = TRUE))
  
  if (is.finite(lo) && is.finite(hi) && lo > 0) {
    if (hi / lo > max_ratio) {
      mid <- sqrt(lo * hi)
      lo  <- mid / sqrt(max_ratio)
      hi  <- mid * sqrt(max_ratio)
    }
  }
  
  span <- hi - lo
  lo2  <- max(lo - pad * span, 1e-4)
  hi2  <- hi + pad * span
  
  c(lo2, hi2)
}

run_meta_logOR_forest <- function(dat,
                                  label = "",
                                  method = "REML",
                                  show_weights = TRUE,
                                  show_or_ci_col = TRUE,
                                  cex = 0.85,
                                  digits = 2,
                                  max_ratio = 80) {
  dat <- dat %>% filter(!is.na(yi), !is.na(sei))
  k <- nrow(dat)
  if (k < 2) return(invisible(NULL))
  
  meta <- metafor::rma(yi = yi, sei = sei, data = dat, method = method)
  
  w_lab <- NULL
  if (isTRUE(show_weights)) {
    vi <- dat$sei^2
    tau2 <- if (is.finite(meta$tau2)) meta$tau2 else 0
    w_raw <- 1 / (vi + tau2)
    w_pct <- round(w_raw / sum(w_raw) * 100, 1)
    w_lab <- sprintf("%.1f%%", w_pct)
  }
  
  or_ci <- NULL
  if (isTRUE(show_or_ci_col)) {
    or_ci <- sprintf(
      paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"),
      exp(dat$yi),
      exp(dat$yi - 1.96 * dat$sei),
      exp(dat$yi + 1.96 * dat$sei)
    )
  }
  
  or_xlim  <- .calc_or_xlim(dat$yi, dat$sei, max_ratio = max_ratio)
  xlim_log <- log(or_xlim)
  
  extra_right <- diff(xlim_log) * 0.80
  xlim_plot   <- c(xlim_log[1], xlim_log[2] + extra_right)
  
  rows    <- k:1
  top_row <- k + 1.2
  
  col1_x <- xlim_log[2] + diff(xlim_log) * 0.10
  col2_x <- xlim_log[2] + diff(xlim_log) * 0.55
  
  metafor::forest(
    meta,
    slab = dat$study,
    xlim = xlim_plot,
    alim = xlim_log,
    atransf = exp,
    at = log(c(0.25, 0.5, 1, 2, 4, 8)),
    rows = rows,
    cex = cex,
    annotate = FALSE,
    header = FALSE,
    xlab = "Odds Ratio",
    col = "red",
    border = "red",
    mlab = ""
  )
  
  title(main = as.character(label)[1], cex.main = 0.90, line = 0.3)
  abline(v = log(1), lty = 2)
  
  if (isTRUE(show_or_ci_col)) {
    text(col1_x, rows, or_ci, pos = 4, cex = cex)
    text(col1_x, top_row, "OR [95% CI]", pos = 4, cex = cex, font = 2)
  }
  
  if (isTRUE(show_weights)) {
    text(col2_x, rows, w_lab, pos = 4, cex = cex)
    text(col2_x, top_row, "Weight", pos = 4, cex = cex, font = 2)
  }
  
  invisible(meta)
}

# ============================================================
# Gather meta-ready input across cohorts
# ============================================================

meta_input <- purrr::imap_dfr(results_by_cohort, function(x, cohort_name) {
  harmonize_meta_ready(x$meta_ready, cohort_name)
})

# ============================================================
# Outcome mapping for the selected phenotype
# ============================================================
# Update only here when switching phenotype.
# ============================================================

outcome_map_by_phenotype <- list(
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

if (!phenotype_to_run %in% names(outcome_map_by_phenotype)) {
  stop("phenotype_to_run not found in outcome_map_by_phenotype.")
}

outcome_map <- outcome_map_by_phenotype[[phenotype_to_run]]

# ============================================================
# Manual correspondence map
# ============================================================
# This is the section to edit when deciding which cohort-specific effects
# represent the same conceptual predictor across studies.
#
# The key idea is explicit harmonisation:
# - you choose the effect_id from each cohort
# - you optionally flip its sign if a reverse direction is needed
#
# Examples below are placeholders and should be replaced with the actual
# effect IDs available in the meta-ready tables.
# ============================================================

correspondence_map <- list(
  
  "Age (older)" = list(
    COHORT_1 = ek("AGE"),
    COHORT_2 = ek("AGE"),
    COHORT_3 = ek("AGE"),
    flip_COHORT_1 = FALSE,
    flip_COHORT_2 = FALSE,
    flip_COHORT_3 = FALSE,
    note = "Age effect"
  ),
  
  "Gender (male)" = list(
    COHORT_1 = ek("SEX", "Male vs Female"),
    COHORT_2 = ek("SEX", "Male vs Female"),
    COHORT_3 = ek("SEX", "Male vs Female"),
    flip_COHORT_1 = FALSE,
    flip_COHORT_2 = FALSE,
    flip_COHORT_3 = FALSE,
    note = "Male versus female"
  ),
  
  "Employment (yes)" = list(
    COHORT_1 = ek("EMPLOYMENT"),
    COHORT_2 = ek("EMPLOYMENT"),
    COHORT_3 = ek("EMPLOYMENT"),
    flip_COHORT_1 = FALSE,
    flip_COHORT_2 = FALSE,
    flip_COHORT_3 = FALSE,
    note = "Employed versus not employed"
  ),
  
  "Functioning (better)" = list(
    COHORT_1 = ek("FUNCTIONING_g"),
    COHORT_2 = ek("FUNCTIONING_g"),
    COHORT_3 = ek("FUNCTIONING_g"),
    flip_COHORT_1 = FALSE,
    flip_COHORT_2 = FALSE,
    flip_COHORT_3 = FALSE,
    note = "Higher functioning"
  ),
  
  "Education (moderate vs low)" = list(
    COHORT_1 = ek("EDUCATION_3CAT_ORD", "Moderate vs Low"),
    COHORT_2 = ek("EDUCATION_3CAT_ORD", "Moderate vs Low"),
    COHORT_3 = ek("EDUCATION_3CAT_ORD", "Moderate vs Low"),
    flip_COHORT_1 = FALSE,
    flip_COHORT_2 = FALSE,
    flip_COHORT_3 = FALSE,
    note = "Moderate versus low education"
  ),
  
  "Education (high vs low)" = list(
    COHORT_1 = ek("EDUCATION_3CAT_ORD", "High vs Low"),
    COHORT_2 = ek("EDUCATION_3CAT_ORD", "High vs Low"),
    COHORT_3 = ek("EDUCATION_3CAT_ORD", "High vs Low"),
    flip_COHORT_1 = FALSE,
    flip_COHORT_2 = FALSE,
    flip_COHORT_3 = FALSE,
    note = "High versus low education"
  ),
  
  "Symptom severity (moderate vs low)" = list(
    COHORT_1 = ek("SEVERITY_3CAT_ORD", "Moderate vs Low"),
    COHORT_2 = ek("SEVERITY_3CAT_ORD", "Moderate vs Low"),
    COHORT_3 = ek("SEVERITY_3CAT_ORD", "Moderate vs Low"),
    flip_COHORT_1 = FALSE,
    flip_COHORT_2 = FALSE,
    flip_COHORT_3 = FALSE,
    note = "Moderate versus low severity"
  ),
  
  "Symptom severity (high vs low)" = list(
    COHORT_1 = ek("SEVERITY_3CAT_ORD", "High vs Low"),
    COHORT_2 = ek("SEVERITY_3CAT_ORD", "High vs Low"),
    COHORT_3 = ek("SEVERITY_3CAT_ORD", "High vs Low"),
    flip_COHORT_1 = FALSE,
    flip_COHORT_2 = FALSE,
    flip_COHORT_3 = FALSE,
    note = "High versus low severity"
  ),
  
  "PRS 1" = list(
    COHORT_1 = ek("PRS_1_resid_z_g"),
    COHORT_2 = ek("PRS_1_resid_z_g"),
    COHORT_3 = ek("PRS_1_resid_z_g"),
    flip_COHORT_1 = FALSE,
    flip_COHORT_2 = FALSE,
    flip_COHORT_3 = FALSE,
    note = "Residualised polygenic score"
  )
)

# ============================================================
# Build meta-analysis dataset
# ============================================================

effects_long <- purrr::imap_dfr(
  correspondence_map,
  ~ build_effect_from_map(.x, .y, meta_input, outcome_map = outcome_map)
)

cat("\n===== CHECK: matched rows by meta-effect and study =====\n")
print(effects_long %>% count(meta_effect, study), n = Inf)

# ============================================================
# Run meta-analyses
# ============================================================

effect_split <- effects_long %>%
  group_by(meta_effect) %>%
  group_split()

results_meta <- list()

for (i in seq_along(effect_split)) {
  dat_i <- effect_split[[i]]
  label_i <- unique(dat_i$meta_effect)
  
  res_i <- run_meta_reml(dat_i)
  results_meta[[label_i]] <- res_i
  
  print_meta_result(res_i, label_i)
}

# ============================================================
# Compact summary table
# ============================================================

meta_summary <- purrr::imap_dfr(results_meta, function(res, nm) {
  if (is.null(res)) {
    return(tibble(
      meta_effect = nm,
      k = 0,
      estimate = NA_real_,
      se = NA_real_,
      pval = NA_real_,
      ci.lb = NA_real_,
      ci.ub = NA_real_,
      OR = NA_real_,
      OR_ci.lb = NA_real_,
      OR_ci.ub = NA_real_,
      tau2 = NA_real_,
      I2 = NA_real_,
      Q = NA_real_,
      Q_df = NA_real_,
      Q_pval = NA_real_
    ))
  }
  
  tibble(
    meta_effect = nm,
    k = res$k,
    estimate = res$pooled$estimate[1],
    se = res$pooled$se[1],
    pval = res$pooled$pval[1],
    ci.lb = res$pooled$ci.lb[1],
    ci.ub = res$pooled$ci.ub[1],
    OR = res$pooled$OR[1],
    OR_ci.lb = res$pooled$OR_ci.lb[1],
    OR_ci.ub = res$pooled$OR_ci.ub[1],
    tau2 = res$heterogeneity$tau2[1],
    I2 = res$heterogeneity$I2[1],
    Q = res$heterogeneity$Q[1],
    Q_df = res$heterogeneity$Q_df[1],
    Q_pval = res$heterogeneity$Q_pval[1]
  )
})

cat("\n===== META-ANALYSIS SUMMARY =====\n")
print(meta_summary, n = Inf, width = Inf)

# ============================================================
# Optional forest plots
# ============================================================
# Example:
# pdf("Meta_Forest_Plots.pdf", width = 13, height = 10, onefile = TRUE)
# effect_list <- effects_long %>% group_by(meta_effect) %>% group_split()
# for (d in effect_list) {
#   if (nrow(d %>% filter(!is.na(yi), !is.na(sei))) >= 2) {
#     run_meta_logOR_forest(
#       dat = d,
#       label = unique(d$meta_effect),
#       method = "REML",
#       cex = 0.90
#     )
#   }
# }
# dev.off()
# ============================================================