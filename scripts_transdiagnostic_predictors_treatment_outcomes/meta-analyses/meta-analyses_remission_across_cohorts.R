# ============================================================
# Remission meta-analysis across CATIE, STEP-BD, and GSRD
# - This script depends on:
#     02_univariable_age_sex_adjusted_meta.R
#   which creates:
#     CATIE_meta_ready
#     STEP_meta_ready
#     GSRD_meta_ready
#
# - Outcome harmonization for remission:
#     CATIE  -> m1_m3_rem
#     STEP   -> remission
#     GSRD   -> remission
#
# - Inputs are age- and sex-adjusted univariable effect estimates
#   (yi = logOR, sei = SE) from the previous script
#
# - Analysis includes:
#     * manual cross-cohort effect matching
#     * optional direction flipping
#     * random-effects meta-analysis (REML)
#     * heterogeneity statistics
#     * leave-one-out analysis when k >= 3
#     * forest plot export
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(metafor)
  library(tibble)
})

# ------------------------------------------------------------
# 0) Load upstream script if needed
# ------------------------------------------------------------
# source("02_univariable_age_sex_adjusted_meta.R")

required_objects <- c("CATIE_meta_ready", "STEP_meta_ready", "GSRD_meta_ready")
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1), inherits = TRUE)]

if (length(missing_objects) > 0) {
  stop(
    "Missing required objects: ",
    paste(missing_objects, collapse = ", "),
    ". Run or source '02_univariable_age_sex_adjusted_meta.R' first."
  )
}

# ------------------------------------------------------------
# 1) Helpers
# ------------------------------------------------------------
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
      contrast_nice = ifelse(is.na(contrast_nice) | contrast_nice == "NA", NA_character_, as.character(contrast_nice)),
      effect_id = ifelse(is.na(contrast_nice), predictor, paste0(predictor, " | ", contrast_nice)),
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

build_effect_from_map <- function(spec, effect_label, pooled_data) {
  d_catie <- apply_flip(
    pick_effect(pooled_data, "CATIE", "m1_m3_rem", spec$CATIE),
    isTRUE(spec$flip_CATIE)
  )
  
  d_step <- apply_flip(
    pick_effect(pooled_data, "STEP", "remission", spec$STEP),
    isTRUE(spec$flip_STEP)
  )
  
  d_gsrd <- apply_flip(
    pick_effect(pooled_data, "GSRD", "remission", spec$GSRD),
    isTRUE(spec$flip_GSRD)
  )
  
  bind_rows(d_catie, d_step, d_gsrd) %>%
    mutate(
      meta_effect = effect_label,
      flipped = case_when(
        study == "CATIE" ~ isTRUE(spec$flip_CATIE),
        study == "STEP"  ~ isTRUE(spec$flip_STEP),
        study == "GSRD"  ~ isTRUE(spec$flip_GSRD),
        TRUE ~ FALSE
      ),
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

# ------------------------------------------------------------
# 2) Forest helpers
# ------------------------------------------------------------
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

plot_forest_grid_pages <- function(effects_long,
                                   ncol = 2,
                                   n_per_page = 4,
                                   cex = 0.90,
                                   method = "REML",
                                   ask = FALSE,
                                   page_title = "Forest plots",
                                   show_effect_list_footer = TRUE) {
  effect_list <- effects_long %>%
    group_by(meta_effect) %>%
    group_split()
  
  effect_list2 <- Filter(function(d) {
    d <- dplyr::filter(d, !is.na(yi), !is.na(sei))
    nrow(d) >= 2
  }, effect_list)
  
  message("Plottabili (k>=2): ", length(effect_list2), " su ", length(effect_list), " effects totali")
  
  if (length(effect_list2) == 0) {
    message("Nessun effect con k>=2: niente forest da plottare.")
    return(invisible(NULL))
  }
  
  pages <- split(effect_list2, ceiling(seq_along(effect_list2) / n_per_page))
  n_pages <- length(pages)
  
  oldpar <- par(no.readonly = TRUE)
  on.exit({
    par(oldpar)
    par(ask = FALSE)
  }, add = TRUE)
  
  par(ask = isTRUE(ask))
  
  for (p in seq_along(pages)) {
    chunk <- pages[[p]]
    nr <- ceiling(length(chunk) / ncol)
    
    par(
      mfrow = c(nr, ncol),
      mar   = c(3.8, 5.2, 2.6, 1.2),
      oma   = c(2.4, 0, 2.4, 0)
    )
    
    eff_names <- character(0)
    
    for (d in chunk) {
      eff <- as.character(unique(d$meta_effect))[1]
      eff_names <- c(eff_names, eff)
      
      d2 <- dplyr::filter(d, !is.na(yi), !is.na(sei))
      if (nrow(d2) < 2) {
        plot.new()
        text(0.5, 0.55, eff, cex = 0.95, font = 2)
        text(0.5, 0.40, "SKIP: k < 2", cex = 0.90)
      } else {
        run_meta_logOR_forest(d2, label = eff, method = method, cex = cex)
      }
    }
    
    mtext(sprintf("%s — page %d/%d", page_title, p, n_pages),
          outer = TRUE, side = 3, cex = 1.10, font = 2, line = 0.8)
    
    if (isTRUE(show_effect_list_footer)) {
      footer_txt <- paste(eff_names, collapse = " | ")
      footer_txt <- paste(strwrap(footer_txt, width = 140), collapse = "\n")
      mtext(footer_txt, outer = TRUE, side = 1, cex = 0.75, line = 0.6)
    }
  }
  
  invisible(NULL)
}

export_forest_pages <- function(effects_long,
                                file_pdf = "Meta_Remission_Forest_Plots.pdf",
                                file_tiff = NULL,
                                ncol = 2,
                                n_per_page = 4,
                                cex = 0.90,
                                method = "REML",
                                width_in = 13,
                                height_in = 10,
                                page_title = "Forest plots",
                                open_pdf_after = TRUE) {
  pdf_path <- normalizePath(file.path(getwd(), file_pdf), winslash = "/", mustWork = FALSE)
  message("Working directory: ", normalizePath(getwd(), winslash = "/"))
  message("PDF target path:  ", pdf_path)
  
  grDevices::pdf(pdf_path, width = width_in, height = height_in, onefile = TRUE)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  
  plot_forest_grid_pages(
    effects_long = effects_long,
    ncol = ncol,
    n_per_page = n_per_page,
    cex = cex,
    method = method,
    ask = FALSE,
    page_title = page_title,
    show_effect_list_footer = TRUE
  )
  
  grDevices::dev.off()
  
  if (!file.exists(pdf_path)) stop("PDF not found where expected: ", pdf_path)
  sz <- file.info(pdf_path)$size
  message("PDF saved: ", pdf_path, " (", round(sz / 1024, 1), " KB)")
  
  if (isTRUE(open_pdf_after)) {
    try(utils::browseURL(pdf_path), silent = TRUE)
  }
  
  if (!is.null(file_tiff)) {
    tiff_path <- normalizePath(file.path(getwd(), file_tiff), winslash = "/", mustWork = FALSE)
    message("TIFF target path: ", tiff_path)
    
    grDevices::tiff(
      tiff_path,
      width = width_in,
      height = height_in,
      units = "in",
      res = 600,
      compression = "lzw"
    )
    on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
    
    plot_forest_grid_pages(
      effects_long = effects_long,
      ncol = ncol,
      n_per_page = n_per_page,
      cex = cex,
      method = method,
      ask = FALSE,
      page_title = page_title,
      show_effect_list_footer = TRUE
    )
    
    grDevices::dev.off()
    
    if (!file.exists(tiff_path)) stop("TIFF not found where expected: ", tiff_path)
    sz2 <- file.info(tiff_path)$size
    message("TIFF saved: ", tiff_path, " (", round(sz2 / 1024, 1), " KB)")
  }
  
  invisible(TRUE)
}

# ------------------------------------------------------------
# 3) Harmonize upstream meta-ready objects
# ------------------------------------------------------------
CATIE_meta <- harmonize_meta_ready(CATIE_meta_ready, "CATIE")
STEP_meta  <- harmonize_meta_ready(STEP_meta_ready,  "STEP")
GSRD_meta  <- harmonize_meta_ready(GSRD_meta_ready,  "GSRD")

meta_input <- bind_rows(CATIE_meta, STEP_meta, GSRD_meta)

# ------------------------------------------------------------
# 4) Manual correspondence map for remission meta-analysis
#    Edit only this section if you want to change mappings
# ------------------------------------------------------------
corresp_remission <- list(
  
  "Age (older)" = list(
    CATIE = ek("AGE"),
    STEP  = ek("AGE"),
    GSRD  = ek("AGE"),
    flip_CATIE = FALSE, flip_STEP = FALSE, flip_GSRD = FALSE,
    note = "Older age"
  ),
  
  "Gender (male)" = list(
    CATIE = ek("SEX", "Male vs Female"),
    STEP  = ek("SEX", "Male vs Female"),
    GSRD  = ek("SEX", "Male vs Female"),
    flip_CATIE = FALSE, flip_STEP = FALSE, flip_GSRD = FALSE,
    note = "Male versus female"
  ),
  
  "Age at onset (higher)" = list(
    CATIE = ek("AGE_ONSET_EMO_clean_g"),
    STEP  = ek("Age_first_depressive_episode_clean_g"),
    GSRD  = ek("Age.at.Onset_clean_g"),
    flip_CATIE = FALSE, flip_STEP = FALSE, flip_GSRD = FALSE,
    note = "Definitions are similar but not identical"
  ),
  
  "Anxiety (moderate vs low)" = list(
    CATIE = ek("Anxiety_3cat_ord", "Moderate vs Low"),
    STEP  = ek("CMF_Anx_3cat_ord", "Moderate vs Low"),
    GSRD  = ek("InnerTension_3cat_ord", "Moderate vs Low"),
    flip_CATIE = FALSE, flip_STEP = FALSE, flip_GSRD = FALSE,
    note = "Inner tension used as anxiety proxy in GSRD"
  ),
  
  "Anxiety (high vs low)" = list(
    CATIE = ek("Anxiety_3cat_ord", "High vs Low"),
    STEP  = ek("CMF_Anx_3cat_ord", "High vs Low"),
    GSRD  = ek("InnerTension_3cat_ord", "High vs Low"),
    flip_CATIE = FALSE, flip_STEP = FALSE, flip_GSRD = FALSE,
    note = "Inner tension used as anxiety proxy in GSRD"
  ),
  
  "Depressive symptoms (higher)" = list(
    CATIE = ek("Calgary_Depressed_Mood_bin"),
    STEP  = ek("MADRS_Total_Baseline_g"),
    GSRD  = ek("MADRS.Retrospective_g"),
    flip_CATIE = FALSE, flip_STEP = FALSE, flip_GSRD = FALSE,
    note = "Different depressive symptom measures"
  ),
  
  "Symptom severity (moderate vs low)" = list(
    CATIE = ek("CGI_Severity_3cat_ord", "Moderate vs Low"),
    STEP  = ek("CGI_baseline_3cat_ord", "Moderate vs Low"),
    GSRD  = ek("MADRS.Retrospective_g"),
    flip_CATIE = FALSE, flip_STEP = FALSE, flip_GSRD = FALSE,
    note = "GSRD uses MADRS retrospective as proxy"
  ),
  
  "Symptom severity (high vs low)" = list(
    CATIE = ek("CGI_Severity_3cat_ord", "High vs Low"),
    STEP  = ek("CGI_baseline_3cat_ord", "High vs Low"),
    GSRD  = ek("MADRS.Retrospective_g"),
    flip_CATIE = FALSE, flip_STEP = FALSE, flip_GSRD = FALSE,
    note = "GSRD uses MADRS retrospective as proxy"
  ),
  
  "Comorbidity (yes)" = list(
    CATIE = ek("Comorbidity"),
    STEP  = ek("Comorbidity"),
    GSRD  = ek("Psychiatric_Comorbidity"),
    flip_CATIE = FALSE, flip_STEP = FALSE, flip_GSRD = FALSE,
    note = "Psychiatric comorbidity"
  ),
  
  "Employment (yes)" = list(
    CATIE = ek("Employment"),
    STEP  = ek("Employment"),
    GSRD  = ek("Employed"),
    flip_CATIE = FALSE, flip_STEP = FALSE, flip_GSRD = FALSE,
    note = "Employed versus not employed"
  ),
  
  "Functioning (better)" = list(
    CATIE = ek("Functioning_g"),
    STEP  = ek("GAF_THIS_MONTH_g"),
    GSRD  = ek("Functioning_g"),
    flip_CATIE = FALSE,
    flip_STEP = FALSE,
    flip_GSRD = TRUE,
    note = "GSRD functioning direction reversed for harmonization"
  ),
  
  "Married (yes)" = list(
    CATIE = ek("Married"),
    STEP  = ek("Married"),
    GSRD  = ek("Marital.Status"),
    flip_CATIE = FALSE, flip_STEP = FALSE, flip_GSRD = FALSE,
    note = "Married versus not married"
  ),
  
  "Living alone (yes)" = list(
    CATIE = ek("SIGOTHER"),
    STEP  = ek("LIVEALONE"),
    GSRD  = ek("Living.with.others"),
    flip_CATIE = TRUE,
    flip_STEP = FALSE,
    flip_GSRD = TRUE,
    note = "CATIE and GSRD require directional harmonization"
  ),
  
  "Ethnicity" = list(
    CATIE = ek("Ethnicity"),
    STEP  = ek("Ethnicity"),
    GSRD  = ek("Ethnicity"),
    flip_CATIE = FALSE, flip_STEP = FALSE, flip_GSRD = FALSE,
    note = "Ethnicity coding may differ across datasets"
  ),
  
  "Education (moderate vs low)" = list(
    CATIE = ek("EDUC_YRS_g"),
    STEP  = ek("Education7_3cat_ord", "Moderate vs Low"),
    GSRD  = ek("Education_3cat_ord", "Moderate vs Low"),
    flip_CATIE = FALSE, flip_STEP = FALSE, flip_GSRD = FALSE,
    note = "CATIE continuous versus STEP/GSRD ordinal"
  ),
  
  "Education (high vs low)" = list(
    CATIE = ek("EDUC_YRS_g"),
    STEP  = ek("Education7_3cat_ord", "High vs Low"),
    GSRD  = ek("Education_3cat_ord", "High vs Low"),
    flip_CATIE = FALSE, flip_STEP = FALSE, flip_GSRD = FALSE,
    note = "CATIE continuous versus STEP/GSRD ordinal"
  ),
  
  "Income (low vs high)" = list(
    CATIE = ek("INCOME_cat", "Low vs High"),
    STEP  = ek("INCOME_3cat_ord", "High vs Low"),
    GSRD  = NA_character_,
    flip_CATIE = FALSE,
    flip_STEP = TRUE,
    flip_GSRD = FALSE,
    note = "STEP flipped to harmonize toward lower income"
  ),
  
  "BMI (medium vs low / higher BMI)" = list(
    CATIE = ek("BMI_3cat_ord", "Medium vs Low"),
    STEP  = NA_character_,
    GSRD  = ek("bmi_g"),
    flip_CATIE = FALSE, flip_STEP = FALSE, flip_GSRD = FALSE,
    note = "Proxy harmonization across available cohorts"
  ),
  
  "BMI (high vs low / higher BMI)" = list(
    CATIE = ek("BMI_3cat_ord", "High vs Low"),
    STEP  = NA_character_,
    GSRD  = ek("bmi_g"),
    flip_CATIE = FALSE, flip_STEP = FALSE, flip_GSRD = FALSE,
    note = "Proxy harmonization across available cohorts"
  ),
  
  "Social withdrawal / poorer social functioning (moderate vs low)" = list(
    CATIE = ek("Social_Withdrawal_3cat_ord", "Moderate vs Low"),
    STEP  = ek("Rel_friends_3cat_ord", "Moderate vs Low"),
    GSRD  = ek("SheehanSocial_3cat_ord", "Moderate vs Low"),
    flip_CATIE = FALSE,
    flip_STEP = TRUE,
    flip_GSRD = TRUE,
    note = "STEP and GSRD flipped to reflect poorer social functioning"
  ),
  
  "Social withdrawal / poorer social functioning (high vs low)" = list(
    CATIE = ek("Social_Withdrawal_3cat_ord", "High vs Low"),
    STEP  = ek("Rel_friends_3cat_ord", "High vs Low"),
    GSRD  = ek("SheehanSocial_3cat_ord", "High vs Low"),
    flip_CATIE = FALSE,
    flip_STEP = TRUE,
    flip_GSRD = TRUE,
    note = "STEP and GSRD flipped to reflect poorer social functioning"
  ),
  
  "Family history (yes)" = list(
    CATIE = ek("Family_History"),
    STEP  = ek("family_history_BD"),
    GSRD  = ek("Family.history.BD"),
    flip_CATIE = FALSE, flip_STEP = FALSE, flip_GSRD = FALSE,
    note = "Not fully identical definitions"
  ),
  
  "Substance abuse (yes)" = list(
    CATIE = ek("Substance_Abuse"),
    STEP  = ek("Substance_Abuse_baseline"),
    GSRD  = ek("Any.Abuse"),
    flip_CATIE = FALSE, flip_STEP = FALSE, flip_GSRD = FALSE,
    note = "Definitions may not be fully identical"
  ),
  
  "Number of previous episodes (moderate vs low)" = list(
    CATIE = ek("N_prev_episodes_3cat_ord", "Moderate vs Low"),
    STEP  = ek("MANPHASE_4cat_ord", "Moderate vs Low"),
    GSRD  = ek("N_dep_episodes_3cat_ord", "Moderate vs Low"),
    flip_CATIE = FALSE, flip_STEP = FALSE, flip_GSRD = FALSE,
    note = "STEP uses lifetime manic phases proxy"
  ),
  
  "Number of previous episodes (high vs low)" = list(
    CATIE = ek("N_prev_episodes_3cat_ord", "High vs Low"),
    STEP  = ek("MANPHASE_4cat_ord", "High vs Low"),
    GSRD  = ek("N_dep_episodes_3cat_ord", "High vs Low"),
    flip_CATIE = FALSE, flip_STEP = FALSE, flip_GSRD = FALSE,
    note = "STEP uses lifetime manic phases proxy"
  ),
  
  "Psychotic symptoms (yes / higher)" = list(
    CATIE = ek("PANSS_Psychotic_g"),
    STEP  = ek("Psychotic_symptoms_bin"),
    GSRD  = ek("Psychotic_MADRS"),
    flip_CATIE = FALSE, flip_STEP = FALSE, flip_GSRD = FALSE,
    note = "Measures differ across cohorts"
  ),
  
  "Suicidality (yes)" = list(
    CATIE = ek("Suicidality_bin"),
    STEP  = ek("Suicide_Attempt"),
    GSRD  = ek("SuicidalRisk_bin"),
    flip_CATIE = FALSE, flip_STEP = FALSE, flip_GSRD = FALSE,
    note = "Constructs differ somewhat across cohorts"
  ),
  
  "PRS ADHD05" = list(
    CATIE = ek("ADHD05_PRS_resid_z_g"),
    STEP  = ek("ADHD05_PRS_resid_z_g"),
    GSRD  = ek("ADHD05_resid_g"),
    flip_CATIE = FALSE, flip_STEP = FALSE, flip_GSRD = FALSE,
    note = "Residualized polygenic score / residual score"
  ),
  
  "PRS BODY04" = list(
    CATIE = ek("BODY04_PRS_resid_z_g"),
    STEP  = ek("BODY04_PRS_resid_z_g"),
    GSRD  = ek("BODY04_resid_g"),
    flip_CATIE = FALSE, flip_STEP = FALSE, flip_GSRD = FALSE,
    note = "Residualized polygenic score / residual score"
  ),
  
  "PRS BP" = list(
    CATIE = ek("BP_PRS_resid_z_g"),
    STEP  = ek("BP_PRS_resid_z_g"),
    GSRD  = ek("BP_resid_g"),
    flip_CATIE = FALSE, flip_STEP = FALSE, flip_GSRD = FALSE,
    note = "Residualized polygenic score / residual score"
  ),
  
  "PRS DEPR07" = list(
    CATIE = ek("DEPR07_PRS_resid_z_g"),
    STEP  = ek("DEPR07_PRS_resid_z_g"),
    GSRD  = ek("DEPR07_resid_g"),
    flip_CATIE = FALSE, flip_STEP = FALSE, flip_GSRD = FALSE,
    note = "Residualized polygenic score / residual score"
  ),
  
  "PRS EDUC03" = list(
    CATIE = ek("EDUC03_PRS_resid_z_g"),
    STEP  = ek("EDUC03_PRS_resid_z_g"),
    GSRD  = ek("EDUC03_resid_g"),
    flip_CATIE = FALSE, flip_STEP = FALSE, flip_GSRD = FALSE,
    note = "Residualized polygenic score / residual score"
  ),
  
  "PRS SCZ" = list(
    CATIE = ek("SCZ_PRS_resid_z_g"),
    STEP  = ek("SCZ_PRS_resid_z_g"),
    GSRD  = ek("SCZ_resid_g"),
    flip_CATIE = FALSE, flip_STEP = FALSE, flip_GSRD = FALSE,
    note = "Residualized polygenic score / residual score"
  )
)

# ------------------------------------------------------------
# 5) Build meta-analysis dataset
# ------------------------------------------------------------
effects_remission <- purrr::imap_dfr(
  corresp_remission,
  ~ build_effect_from_map(.x, .y, meta_input)
)

cat("\n===== CHECK: matched rows by meta-effect and study =====\n")
print(effects_remission %>% count(meta_effect, study), n = Inf)

# ------------------------------------------------------------
# 6) Run meta-analyses
# ------------------------------------------------------------
effect_split <- effects_remission %>%
  group_by(meta_effect) %>%
  group_split()

results_remission <- list()

for (i in seq_along(effect_split)) {
  dat_i <- effect_split[[i]]
  label_i <- unique(dat_i$meta_effect)
  
  res_i <- run_meta_reml(dat_i)
  results_remission[[label_i]] <- res_i
  
  print_meta_result(res_i, label_i)
}

# ------------------------------------------------------------
# 7) Compact summary object
# ------------------------------------------------------------
meta_summary_remission <- purrr::imap_dfr(results_remission, function(res, nm) {
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

cat("\n===== META-ANALYSIS SUMMARY: REMISSION =====\n")
print(meta_summary_remission, n = Inf, width = Inf)

# ------------------------------------------------------------
# 8) Export forest plots
# ------------------------------------------------------------
export_forest_pages(
  effects_long = effects_remission,
  file_pdf = "Meta_Remission_Forest_Plots.pdf",
  ncol = 2,
  n_per_page = 4,
  cex = 0.90,
  method = "REML",
  width_in = 13,
  height_in = 10,
  page_title = "Outcome: Remission — Forest plots"
)