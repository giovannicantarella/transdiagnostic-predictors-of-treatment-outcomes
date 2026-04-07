# ============================================================
# CATIE / STEP-BD / GSRD | Age- and sex-adjusted univariate models
#
# - Each candidate predictor is tested in an adjusted univariate model:
#     Outcome ~ predictor + AGE + SEX
#
# - Continuous predictors are Gelman-rescaled within each dataset:
#     X_g = (X - mean(X)) / (2 * SD(X))
#
# - Binary predictors are entered as coded
# - Ordinal predictors are entered as ordered factors with treatment contrasts
#
# - AGE and SEX are also estimated via minimally adjusted models:
#     Outcome ~ AGE + SEX
#   so that AGE is adjusted for SEX and SEX is adjusted for AGE
#
# - The script returns:
#     * dataset-specific adjusted univariate tables
#     * meta-analysis-ready tables (yi = logOR, sei = SE)
#
# Outcome naming used here:
# - CATIE:   TR1, TR2, remission
# - STEP:    TR1, TR2, remission
# - GSRD:    Resistance, remission
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(broom)
  library(stringr)
  library(purrr)
})

# ============================================================
# Shared helpers
# ============================================================

gelman_2sd <- function(x) {
  x <- as.numeric(x)
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - m) / (2 * s)
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

fit_one <- function(data, outcome, predictor) {
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

# ============================================================
# CATIE
# Outcomes: TR1, TR2, remission
# ============================================================

df <- CATIE

df$AGE <- df$AGE
df$SEX <- ifelse(
  df$GENDER == "Male", 1,
  ifelse(df$GENDER == "Female", 0, NA)
)

if ("BMI_3cat_ord" %in% names(df)) {
  df$BMI_3cat_ord <- ordered(
    df$BMI_3cat_ord,
    levels = c("Normal", "Overweight", "Obese"),
    labels = c("Low", "Medium", "High")
  )
}

cont_vars <- c(
  "EDUC_YRS",
  "AGE_ONSET_EMO_clean",
  "Functioning",
  "PANSS_Psychotic",
  "ADHD05_PRS_resid_z",
  "DEPR07_PRS_resid_z",
  "EDUC03_PRS_resid_z",
  "SCZ_PRS_resid_z",
  "BP_PRS_resid_z",
  "BODY04_PRS_resid_z"
)

bin_vars <- c(
  "Ethnicity",
  "SIGOTHER",
  "Married",
  "Employment",
  "Family_History",
  "Substance_Abuse",
  "Comorbidity",
  "Calgary_Depressed_Mood_bin",
  "Suicidality_bin"
)

ord_spec <- list(
  INCOME_cat                 = c("High", "Low"),
  BMI_3cat_ord               = c("Low", "Medium", "High"),
  Social_Withdrawal_3cat_ord = c("Low", "Moderate", "High"),
  N_prev_episodes_3cat_ord   = c("Low", "Moderate", "High"),
  CGI_Severity_3cat_ord      = c("Low", "Moderate", "High"),
  Anxiety_3cat_ord           = c("Low", "Moderate", "High")
)

for (v in names(ord_spec)) {
  if (!v %in% names(df)) next
  df[[v]] <- ordered(df[[v]], levels = ord_spec[[v]])
  df[[v]] <- droplevels(df[[v]])
  if (nlevels(df[[v]]) >= 2) {
    contrasts(df[[v]]) <- contr.treatment(nlevels(df[[v]]), base = 1)
  }
}

cont_vars <- cont_vars[cont_vars %in% names(df)]
for (v in cont_vars) {
  df[[paste0(v, "_g")]] <- gelman_2sd(df[[v]])
}

outcomes <- c("TR1", "TR2", "remission")
outcomes <- outcomes[outcomes %in% names(df)]

CATIE_age_sex_minadj <- purrr::map_dfr(
  outcomes,
  ~ fit_age_sex_minadj(df, .x, sex_contrast_label = "Male vs Female")
)

cont_pred <- paste0(cont_vars, "_g")
bin_pred  <- bin_vars[bin_vars %in% names(df)]
ord_pred  <- names(ord_spec)
ord_pred  <- ord_pred[ord_pred %in% names(df)]

all_predictors <- c(cont_pred, bin_pred, ord_pred)

CATIE_adjusted_univariate <- bind_rows(
  purrr::map_dfr(outcomes, function(out) {
    purrr::map_dfr(all_predictors, function(pred) fit_one(df, out, pred))
  }),
  CATIE_age_sex_minadj
) %>%
  dplyr::mutate(
    study = "CATIE",
    predictor_type = dplyr::case_when(
      predictor %in% c("AGE", "SEX") ~ "demographic (minimally adjusted: AGE + SEX)",
      predictor %in% cont_pred       ~ "continuous (Gelman 2SD)",
      predictor %in% bin_pred        ~ "binary",
      predictor %in% ord_pred        ~ "ordinal (treatment contrasts)",
      TRUE                           ~ "other"
    ),
    contrast_nice = dplyr::case_when(
      predictor %in% ord_pred ~ purrr::map2_chr(
        predictor, term, ~ level_label_from_term_ord(.x, .y, ord_spec)
      ),
      predictor == "SEX" ~ dplyr::coalesce(contrast_nice, "Male vs Female"),
      TRUE ~ NA_character_
    ),
    p_value_print = format.pval(p_value, digits = 8, eps = 1e-12)
  ) %>%
  dplyr::arrange(outcome, predictor, term)

CATIE_meta_ready <- CATIE_adjusted_univariate %>%
  dplyr::mutate(
    effect_id = ifelse(
      is.na(contrast_nice),
      predictor,
      paste0(predictor, " | ", contrast_nice)
    ),
    yi  = logOR,
    sei = SE
  ) %>%
  dplyr::select(study, outcome, effect_id, predictor, contrast_nice, n, yi, sei, predictor_type)

# ============================================================
# STEP-BD
# Outcomes: TR1, TR2, remission
# ============================================================

df <- STEP

if (!"AGE_AT_RANDOMIZATION" %in% names(df)) stop("Missing AGE_AT_RANDOMIZATION in STEP")
df$AGE <- df$AGE_AT_RANDOMIZATION

if (!"GENDER" %in% names(df)) stop("Missing GENDER in STEP")
df$SEX <- ifelse(
  df$GENDER %in% c("Male", "M", 1), 1,
  ifelse(df$GENDER %in% c("Female", "F", 0), 0, NA)
)

outcomes <- c("TR1", "TR2", "remission")
outcomes <- outcomes[outcomes %in% names(df)]
if (length(outcomes) == 0) stop("None of the outcomes (TR1, TR2, remission) found in STEP")

STEP_age_sex_minadj <- purrr::map_dfr(
  outcomes,
  ~ fit_age_sex_minadj(df, .x, sex_contrast_label = "Male vs Female")
)

df$Female <- ifelse(
  df$GENDER %in% c("Female", "F", 0), 1,
  ifelse(df$GENDER %in% c("Male", "M", 1), 0, NA)
)

bin_vars <- c(
  "Female",
  "Ethnicity",
  "Employment",
  "Married",
  "LIVEALONE",
  "family_history_BD",
  "Substance_Abuse_baseline",
  "Comorbidity",
  "Suicide_Attempt",
  "Psychotic_symptoms_bin"
)
bin_vars <- bin_vars[bin_vars %in% names(df)]

cont_vars <- c(
  "Age_first_depressive_episode_clean",
  "GAF_THIS_MONTH",
  "MADRS_Total_Baseline",
  "ADHD05_PRS_resid_z",
  "DEPR07_PRS_resid_z",
  "EDUC03_PRS_resid_z",
  "SCZ_PRS_resid_z",
  "BP_PRS_resid_z",
  "BODY04_PRS_resid_z"
)
cont_vars <- cont_vars[cont_vars %in% names(df)]

for (v in cont_vars) {
  df[[paste0(v, "_g")]] <- gelman_2sd(df[[v]])
}
cont_pred <- paste0(cont_vars, "_g")

ord_spec <- list(
  Education7_3cat_ord   = c("Low", "Moderate", "High"),
  INCOME_3cat_ord       = c("Low", "Moderate", "High"),
  MANPHASE_4cat_ord     = c("Low", "Moderate", "High", "Very high"),
  CGI_baseline_3cat_ord = c("Low", "Moderate", "High"),
  Rel_friends_3cat_ord  = c("Low", "Moderate", "High"),
  CMF_Anx_3cat_ord      = c("Low", "Moderate", "High")
)

ord_pred <- names(ord_spec)
ord_pred <- ord_pred[ord_pred %in% names(df)]

for (v in ord_pred) {
  df[[v]] <- ordered(df[[v]], levels = ord_spec[[v]])
  df[[v]] <- droplevels(df[[v]])
  if (nlevels(df[[v]]) >= 2) {
    contrasts(df[[v]]) <- contr.treatment(nlevels(df[[v]]), base = 1)
  }
}

all_predictors <- c(cont_pred, bin_vars, ord_pred)

STEP_adjusted_univariate <- bind_rows(
  purrr::map_dfr(outcomes, function(out) {
    purrr::map_dfr(all_predictors, function(pred) fit_one(df, out, pred))
  }),
  STEP_age_sex_minadj
) %>%
  dplyr::mutate(
    study = "STEP",
    predictor_type = dplyr::case_when(
      predictor %in% c("AGE", "SEX") ~ "demographic (minimally adjusted: AGE + SEX)",
      predictor %in% cont_pred       ~ "continuous (Gelman 2SD)",
      predictor %in% bin_vars        ~ "binary",
      predictor %in% ord_pred        ~ "ordinal (treatment contrasts)",
      TRUE                           ~ "other"
    ),
    contrast_nice = dplyr::case_when(
      predictor %in% ord_pred ~ purrr::map2_chr(
        predictor, term, ~ level_label_from_term_ord(.x, .y, ord_spec)
      ),
      predictor == "SEX" ~ dplyr::coalesce(contrast_nice, "Male vs Female"),
      TRUE ~ NA_character_
    ),
    p_value_print = format.pval(p_value, digits = 6, eps = 1e-12)
  ) %>%
  dplyr::arrange(outcome, predictor, term)

STEP_meta_ready <- STEP_adjusted_univariate %>%
  dplyr::mutate(
    effect_id = ifelse(
      is.na(contrast_nice),
      predictor,
      paste0(predictor, " | ", contrast_nice)
    ),
    yi  = logOR,
    sei = SE
  ) %>%
  dplyr::select(study, outcome, effect_id, predictor, contrast_nice, n, yi, sei, predictor_type)

# ============================================================
# GSRD
# Outcomes: Resistance, remission
# ============================================================

df <- GSRD

predictors_in <- c(
  "Sex",
  "Marital.Status",
  "Living.with.others",
  "Ethnicity",
  "Family.history.MDD",
  "Family.history.BD",
  "Psychotic_MADRS",
  "Any.Abuse",
  "Psychiatric_Comorbidity",
  "Age",
  "Age.at.Onset_clean",
  "MADRS.Retrospective",
  "InnerTension_3cat_ord",
  "Education_3cat_ord",
  "Employed",
  "bmi",
  "SheehanWork_3cat_ord",
  "SheehanSocial_3cat_ord",
  "SheehanFamily_3cat_ord",
  "Functioning",
  "N_dep_episodes_3cat_ord",
  "SuicidalRisk_bin",
  "ADHD05_resid",
  "DEPR07_resid",
  "EDUC03_resid",
  "SCZ_resid",
  "BP_resid",
  "BODY04_resid"
)

outcomes <- c("Resistance", "remission")
outcomes <- outcomes[outcomes %in% names(df)]
if (length(outcomes) == 0) stop("GSRD: none of the outcomes found (Resistance, remission).")

predictors_in <- predictors_in[predictors_in %in% names(df)]
if (length(predictors_in) == 0) stop("GSRD: none of the listed predictors found.")

if (!"Age" %in% names(df)) stop("GSRD: missing Age")
df$AGE <- df$Age

if (!"Sex" %in% names(df)) stop("GSRD: missing Sex")
df$SEX <- NA_real_
df$SEX[df$Sex %in% c("Male", "M", "male", "m", 1, "1")] <- 1
df$SEX[df$Sex %in% c("Female", "F", "female", "f", 0, "0")] <- 0
if (all(is.na(df$SEX)) && is.numeric(df$Sex)) df$SEX <- df$Sex

GSRD_age_sex_minadj <- purrr::map_dfr(
  outcomes,
  ~ fit_age_sex_minadj(df, .x, sex_contrast_label = "Male vs Female")
)

predictors_model <- setdiff(predictors_in, c("Age", "Sex"))

ord_spec <- list(
  InnerTension_3cat_ord   = c("Low", "Moderate", "High"),
  Education_3cat_ord      = c("Low", "Moderate", "High"),
  SheehanWork_3cat_ord    = c("Low", "Moderate", "High"),
  SheehanSocial_3cat_ord  = c("Low", "Moderate", "High"),
  SheehanFamily_3cat_ord  = c("Low", "Moderate", "High"),
  N_dep_episodes_3cat_ord = c("Low", "Moderate", "High")
)
ord_pred <- intersect(names(ord_spec), predictors_model)

bin_pred <- intersect(
  c(
    "Living.with.others",
    "Marital.Status",
    "Ethnicity",
    "Family.history.MDD",
    "Family.history.BD",
    "Psychotic_MADRS",
    "Any.Abuse",
    "Psychiatric_Comorbidity",
    "Employed",
    "SuicidalRisk_bin"
  ),
  predictors_model
)

cont_pred <- setdiff(predictors_model, c(ord_pred, bin_pred))
cont_pred_num <- cont_pred[sapply(cont_pred, function(v) is.numeric(df[[v]]))]
cat_pred <- setdiff(predictors_model, c(ord_pred, bin_pred, cont_pred_num))

for (v in ord_pred) {
  df[[v]] <- ordered(df[[v]], levels = ord_spec[[v]])
  df[[v]] <- droplevels(df[[v]])
  if (nlevels(df[[v]]) >= 2) {
    contrasts(df[[v]]) <- contr.treatment(nlevels(df[[v]]), base = 1)
  }
}

for (v in bin_pred) {
  df[[v]] <- to_binary_01(df[[v]])
}

for (v in cont_pred_num) {
  df[[paste0(v, "_g")]] <- gelman_2sd(df[[v]])
}
cont_pred_g <- paste0(cont_pred_num, "_g")

for (v in cat_pred) {
  df[[v]] <- droplevels(factor(df[[v]]))
  if (nlevels(df[[v]]) >= 2) {
    contrasts(df[[v]]) <- contr.treatment(nlevels(df[[v]]), base = 1)
  }
}

all_predictors <- c(cont_pred_g, bin_pred, ord_pred, cat_pred)

GSRD_adjusted_univariate <- bind_rows(
  purrr::map_dfr(outcomes, function(out) {
    purrr::map_dfr(all_predictors, function(pred) fit_one(df, out, pred))
  }),
  GSRD_age_sex_minadj
) %>%
  dplyr::mutate(
    study = "GSRD",
    predictor_type = dplyr::case_when(
      predictor %in% c("AGE", "SEX") ~ "demographic (minimally adjusted: AGE + SEX)",
      predictor %in% cont_pred_g     ~ "continuous (Gelman 2SD)",
      predictor %in% bin_pred        ~ "binary",
      predictor %in% ord_pred        ~ "ordinal (treatment contrasts)",
      predictor %in% cat_pred        ~ "categorical (treatment contrasts)",
      TRUE                           ~ "other"
    ),
    contrast_nice = dplyr::case_when(
      predictor %in% ord_pred ~ purrr::map2_chr(
        predictor, term, ~ level_label_from_term_ord(.x, .y, ord_spec)
      ),
      predictor %in% cat_pred ~ purrr::map2_chr(
        predictor, term, ~ level_label_from_term_cat(df, .x, .y)
      ),
      predictor == "SEX" ~ dplyr::coalesce(contrast_nice, "Male vs Female"),
      TRUE ~ NA_character_
    ),
    p_value_print = format.pval(p_value, digits = 8, eps = 1e-12)
  ) %>%
  dplyr::arrange(outcome, predictor, term)

GSRD_meta_ready <- GSRD_adjusted_univariate %>%
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

# ============================================================
# Final objects created
# ============================================================
# CATIE_adjusted_univariate
# STEP_adjusted_univariate
# GSRD_adjusted_univariate
#
#     CATIE_meta_ready
#     STEP_meta_ready
#     GSRD_meta_ready
#
