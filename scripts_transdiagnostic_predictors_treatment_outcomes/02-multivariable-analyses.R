
# ============================================================
# Project: Transdiagnostic predictors of treatment outcomes
# Script: Multivariable modelling workflow within each cohort
#
# Purpose:
# - fit penalised multivariable logistic regression models
# - tune alpha using nested cross-validation
# - estimate out-of-fold discrimination
# - select a sparse model using lambda.1se
# - assess calibration using out-of-fold predictions
#
# Workflow:
# - predictors are assumed to have been preselected upstream
#   (for example from prior evidence, harmonisation choices,
#   or univariable screening rules defined elsewhere)
# - OUTER CV provides out-of-fold predictions used to estimate AUC
# - INNER CV selects lambda.1se within each outer training fold
# - the final model is refit on the full complete-case dataset
# - calibration plots are generated using out-of-fold predictions
#
# Notes:
# - this script is generic and should be adapted only at the level of:
#   1) input dataset
#   2) outcomes
#   3) predictor sets
# - any cohort-specific preprocessing, exclusions, outlier handling,
#   or recoding should already have been completed before this step
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(glmnet)
  library(pROC)
  library(ggplot2)
  library(scales)
})

set.seed(1234)
options(contrasts = c("contr.treatment", "contr.poly"))

# ============================================================
# Input dataset
# ============================================================
# Replace this with the analytic dataset for the cohort of interest.
# The dataset should contain one row per participant and already include
# the outcomes and predictors to be modelled.
# ============================================================

analysis_df <- cohort_df

# ============================================================
# Settings
# ============================================================

alpha_grid   <- seq(0.10, 0.30, by = 0.05)
nfolds_outer <- 5
nfolds_inner <- 5
auc_tol      <- 0.01

# Optional variables to exclude globally
drop_vars <- character(0)

# Whether ordered 3-level predictors should be encoded
# with polynomial contrasts
use_poly_for_ordinals <- TRUE
ordinal_level_order_3 <- c("Low", "Moderate", "High")

# Calibration settings
n_calibration_bins <- 5
calibration_dir <- "Outputs"

# ============================================================
# Predictor sets by outcome
# ============================================================
# Replace these placeholders with the predictor sets you want to evaluate.
# This structure allows one different predictor set per outcome, but the
# same workflow is applied throughout.
#
# Example:
# outcome_map <- list(
#   remission = c("AGE", "SEX", "SEVERITY_3CAT_ORD", "FUNCTIONING"),
#   resistance = c("AGE", "SEX", "SOCIAL_WITHDRAWAL_3CAT_ORD", "PRS_1_resid_z")
# )
# ============================================================

pred_vars_outcome1 <- c(
  "AGE",
  "SEX",
  "EMPLOYMENT",
  "SEVERITY_3CAT_ORD",
  "FUNCTIONING",
  "PRS_1_resid_z"
)

pred_vars_outcome2 <- c(
  "AGE",
  "SEX",
  "SOCIAL_WITHDRAWAL_3CAT_ORD",
  "ANXIETY_3CAT_ORD",
  "PRS_2_resid_z"
)

outcome_map <- list(
  OUTCOME_1 = pred_vars_outcome1,
  OUTCOME_2 = pred_vars_outcome2
)

outcome_map <- lapply(
  outcome_map,
  function(v) setdiff(unique(v), drop_vars)
)

cat("Elastic net alpha-tuning (nested CV + Gelman rescaling, lambda.1se)\n")
cat("Alpha grid:", paste(alpha_grid, collapse = ", "), "\n")
cat("Outer folds:", nfolds_outer,
    " | Inner folds:", nfolds_inner,
    " | auc_tol:", auc_tol, "\n")
cat("Ordinal polynomial contrasts:", use_poly_for_ordinals, "\n\n")

# ============================================================
# Helper functions
# ============================================================

# Convert common binary encodings to 0/1 integer
to01 <- function(x) {
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x)) return(as.integer(x))

  x <- as.character(x)

  as.integer(
    x %in% c(
      "1", "TRUE", "T", "true",
      "Yes", "YES", "case", "Case"
    )
  )
}

# Create approximately stratified folds for binary outcomes
make_strat_folds <- function(y, K) {
  y <- as.integer(y)

  if (length(unique(y)) < 2) {
    stop("make_strat_folds: y must contain both classes.")
  }

  idx1 <- which(y == 1)
  idx0 <- which(y == 0)

  foldid <- integer(length(y))

  foldid[idx1] <- sample(
    rep(seq_len(K), length.out = length(idx1))
  )

  foldid[idx0] <- sample(
    rep(seq_len(K), length.out = length(idx0))
  )

  foldid
}

# Choose the final alpha using:
# 1) prefer models with at least one selected predictor
# 2) highest outer AUC
# 3) within auc_tol, choose smallest number of selected predictors
# 4) tie-break by higher AUC, then lower alpha
choose_alpha <- function(df_summary, auc_tol = 0.01) {

  df_summary <- df_summary %>%
    mutate(
      n_selected = as.integer(n_selected),
      AUC = as.numeric(AUC)
    )

  df_pos <- df_summary %>%
    filter(n_selected > 0)

  if (nrow(df_pos) == 0) {

    cat("[WARN] All models have n_selected == 0. Falling back to highest AUC.\n")

    auc_best <- max(df_summary$AUC, na.rm = TRUE)

    keep <- df_summary %>%
      filter(AUC >= (auc_best - auc_tol))

    return(
      keep %>%
        arrange(desc(AUC), alpha) %>%
        slice(1)
    )
  }

  auc_best <- max(df_pos$AUC, na.rm = TRUE)

  keep <- df_pos %>%
    filter(AUC >= (auc_best - auc_tol))

  keep %>%
    arrange(n_selected, desc(AUC), alpha) %>%
    slice(1)
}

# Apply polynomial contrasts to ordered 3-level predictors
# This is useful when ordinal grouped predictors are modelled
# as ordered factors
apply_poly_contrasts <- function(
    df,
    predictors_in,
    level_order_3 = c("Low", "Moderate", "High")
) {

  ord_vars <- intersect(
    grep("_3cat_ord$", predictors_in, value = TRUE),
    names(df)
  )

  if (length(ord_vars) == 0) return(df)

  for (v in ord_vars) {

    if (is.numeric(df[[v]])) next

    x <- as.character(df[[v]])

    if (!is.null(level_order_3) &&
        length(level_order_3) == 3) {

      vals <- unique(na.omit(x))

      if (length(vals) > 0 &&
          all(vals %in% level_order_3)) {

        df[[v]] <- ordered(
          x,
          levels = level_order_3
        )

      } else {

        df[[v]] <- ordered(df[[v]])
      }

    } else {

      df[[v]] <- ordered(df[[v]])
    }

    k <- length(levels(df[[v]]))

    if (k >= 2) {
      contrasts(df[[v]]) <- contr.poly(k)
    }
  }

  df
}

# Identify columns that should be treated as continuous
# for Gelman scaling.
# Binary dummy columns are not rescaled.
get_continuous_cols <- function(X) {

  X <- as.matrix(X)

  is_cont <- logical(ncol(X))

  for (j in seq_len(ncol(X))) {

    x <- X[, j]
    ux <- unique(x[is.finite(x)])

    if (length(ux) <= 1) {
      is_cont[j] <- FALSE
      next
    }

    if (all(ux %in% c(0, 1))) {
      is_cont[j] <- FALSE
      next
    }

    is_cont[j] <- TRUE
  }

  which(is_cont)
}

# Fit Gelman scaling parameters on training data only
gelman_fit <- function(X, cont_idx) {

  X <- as.matrix(X)

  mu <- rep(0, ncol(X))
  sc <- rep(1, ncol(X))

  if (length(cont_idx) > 0) {

    mu[cont_idx] <- colMeans(
      X[, cont_idx, drop = FALSE],
      na.rm = TRUE
    )

    sdv <- apply(
      X[, cont_idx, drop = FALSE],
      2,
      sd,
      na.rm = TRUE
    )

    sc[cont_idx] <- 2 * sdv

    sc[cont_idx][
      !is.finite(sc[cont_idx]) |
      sc[cont_idx] == 0
    ] <- 1
  }

  list(
    center = mu,
    scale = sc
  )
}

# Apply previously estimated Gelman scaling
gelman_apply <- function(X, sc) {

  X <- as.matrix(X)

  sweep(
    sweep(X, 2, sc$center, "-"),
    2,
    sc$scale,
    "/"
  )
}

# ============================================================
# Nested cross-validation for one alpha
# ============================================================
# This function:
# - uses outer folds to produce out-of-fold predictions
# - performs inner CV within each outer training fold to choose lambda.1se
# - returns the outer-CV AUC and prediction vector
# ============================================================

nested_auc_alpha <- function(
    X_raw,
    y,
    alpha,
    foldid_outer,
    nfolds_inner,
    cont_idx
) {

  n <- nrow(X_raw)

  stopifnot(
    length(y) == n,
    length(foldid_outer) == n
  )

  prob_oof <- rep(NA_real_, n)

  for (k in sort(unique(foldid_outer))) {

    idx_te <- which(foldid_outer == k)
    idx_tr <- which(foldid_outer != k)

    Xtr_raw <- X_raw[idx_tr, , drop = FALSE]
    Xte_raw <- X_raw[idx_te, , drop = FALSE]

    ytr <- y[idx_tr]

    # Fold-safe scaling: estimated only on outer training data
    sc <- gelman_fit(Xtr_raw, cont_idx)

    Xtr <- gelman_apply(Xtr_raw, sc)
    Xte <- gelman_apply(Xte_raw, sc)

    storage.mode(Xtr) <- "double"
    storage.mode(Xte) <- "double"

    foldid_inner <- make_strat_folds(
      ytr,
      K = nfolds_inner
    )

    cv_in <- cv.glmnet(
      x = Xtr,
      y = ytr,
      family = "binomial",
      alpha = alpha,
      foldid = foldid_inner,
      standardize = FALSE
    )

    fit_in <- glmnet(
      x = Xtr,
      y = ytr,
      family = "binomial",
      alpha = alpha,
      lambda = cv_in$lambda.1se,
      standardize = FALSE
    )

    prob_oof[idx_te] <- as.numeric(
      predict(
        fit_in,
        newx = Xte,
        type = "response"
      )
    )
  }

  roc_obj <- pROC::roc(
    y,
    prob_oof,
    quiet = TRUE
  )

  auc_val <- as.numeric(
    pROC::auc(roc_obj)
  )

  auc_ci <- as.numeric(
    pROC::ci.auc(roc_obj)
  )

  list(
    AUC = auc_val,
    AUC_low = auc_ci[1],
    AUC_high = auc_ci[3],
    prob_oof = prob_oof
  )
}

# ============================================================
# Fit one elastic net model for one alpha
# ============================================================
# This function:
# - filters to the requested outcome and predictors
# - performs complete-case analysis for the selected variables
# - builds the model matrix
# - estimates nested-CV AUC for the chosen alpha
# - refits the final model on the full complete-case sample
# ============================================================

fit_enet_one_alpha <- function(
    data,
    outcome_var,
    predictors,
    alpha,
    foldid_outer,
    nfolds_inner,
    use_poly_for_ordinals,
    ordinal_level_order_3,
    drop_vars
) {

  predictors_in <- unique(predictors)

  predictors_in <- setdiff(
    predictors_in,
    drop_vars
  )

  missing_pred <- setdiff(
    predictors_in,
    names(data)
  )

  predictors_in <- setdiff(
    predictors_in,
    missing_pred
  )

  if (!(outcome_var %in% names(data))) {
    stop("Outcome not found in dataset: ", outcome_var)
  }

  if (length(predictors_in) == 0) {
    stop("No predictors available after filtering.")
  }

  df_raw <- data %>%
    select(all_of(c(outcome_var, predictors_in)))

  n_total <- nrow(df_raw)

  df <- df_raw %>%
    drop_na()

  if (nrow(df) == 0) {
    stop("No complete-case rows after drop_na().")
  }

  if (isTRUE(use_poly_for_ordinals)) {

    df <- apply_poly_contrasts(
      df,
      predictors_in,
      level_order_3 = ordinal_level_order_3
    )
  }

  y <- to01(df[[outcome_var]])

  if (length(unique(y)) < 2) {
    stop("Outcome has no variation after complete-case filtering.")
  }

  n_cc <- nrow(df)

  n_events <- sum(y == 1)
  n_nonevents <- sum(y == 0)

  X_raw <- model.matrix(
    reformulate(predictors_in),
    data = df
  )[, -1, drop = FALSE]

  if (length(foldid_outer) != n_cc) {
    stop(
      "foldid_outer length does not match the number of complete-case observations."
    )
  }

  cont_idx <- get_continuous_cols(X_raw)

  nested <- nested_auc_alpha(
    X_raw = X_raw,
    y = y,
    alpha = alpha,
    foldid_outer = foldid_outer,
    nfolds_inner = nfolds_inner,
    cont_idx = cont_idx
  )

  # Final full-data fit using the same scaling logic
  sc_full <- gelman_fit(
    X_raw,
    cont_idx
  )

  X <- gelman_apply(
    X_raw,
    sc_full
  )

  storage.mode(X) <- "double"

  foldid_inner_full <- make_strat_folds(
    y,
    K = nfolds_inner
  )

  cvfit <- cv.glmnet(
    x = X,
    y = y,
    family = "binomial",
    alpha = alpha,
    foldid = foldid_inner_full,
    standardize = FALSE
  )

  fit <- glmnet(
    x = X,
    y = y,
    family = "binomial",
    alpha = alpha,
    lambda = cvfit$lambda.1se,
    standardize = FALSE
  )

  b <- as.matrix(coef(fit))

  selected <- setdiff(
    rownames(b)[as.numeric(b[, 1]) != 0],
    "(Intercept)"
  )

  list(
    outcome = outcome_var,
    alpha = alpha,
    n_total = n_total,
    n_cc = n_cc,
    events = n_events,
    nonevents = n_nonevents,
    predictors_used = predictors_in,
    missing_pred = missing_pred,
    lambda_min = cvfit$lambda.min,
    lambda_1se = cvfit$lambda.1se,
    n_selected = length(selected),
    selected_terms = selected,
    AUC = nested$AUC,
    AUC_low = nested$AUC_low,
    AUC_high = nested$AUC_high,
    prob_oof = nested$prob_oof,
    y = y,
    fit = fit,
    cvfit = cvfit,
    X = X
  )
}

# ============================================================
# Print chosen model
# ============================================================
# Prints a compact summary of the selected model, including the
# non-zero coefficients sorted by absolute magnitude.
# ============================================================

print_chosen_model <- function(pick_obj, chosen_row) {

  cat(
    "\nCHOSEN | ", pick_obj$outcome,
    " | alpha=", sprintf("%.2f", chosen_row$alpha),
    " | AUC(outer)=", sprintf("%.3f", chosen_row$AUC),
    " [", sprintf("%.3f", chosen_row$AUC_low),
    ", ", sprintf("%.3f", chosen_row$AUC_high), "]",
    " | n_sel=", chosen_row$n_selected,
    " | lambda.1se(full fit)=",
    sprintf("%.6f", pick_obj$lambda_1se),
    "\n",
    sep = ""
  )

  b <- as.matrix(
    coef(pick_obj$fit)
  )

  b_df <- data.frame(
    term = rownames(b),
    beta = as.numeric(b[, 1]),
    stringsAsFactors = FALSE
  )

  b_df <- b_df[
    b_df$beta != 0,
    ,
    drop = FALSE
  ]

  b_df <- b_df[
    order(-abs(b_df$beta)),
    ,
    drop = FALSE
  ]

  cat("Non-zero coefficients (sorted by absolute beta):\n")

  if (nrow(b_df) == 0) {

    cat("(NONE)\n")

  } else {

    print(
      b_df,
      row.names = FALSE
    )
  }

  cat("Summary:\n")

  cat(" - N_total: ", pick_obj$n_total, "\n", sep = "")
  cat(" - N_cc: ", pick_obj$n_cc, "\n", sep = "")

  cat(
    " - Events: ", pick_obj$events,
    " | Nonevents: ", pick_obj$nonevents,
    "\n",
    sep = ""
  )

  cat(
    " - Predictors used: ",
    length(pick_obj$predictors_used),
    "\n",
    sep = ""
  )

  cat(
    " - Design matrix columns: ",
    ncol(pick_obj$X),
    "\n",
    sep = ""
  )

  cat(
    " - lambda.min: ",
    sprintf("%.6f", pick_obj$lambda_min),
    "\n",
    sep = ""
  )

  cat(
    " - lambda.1se: ",
    sprintf("%.6f", pick_obj$lambda_1se),
    "\n",
    sep = ""
  )

  cat(
    " - n_nonzero (excluding intercept): ",
    pick_obj$n_selected,
    "\n",
    sep = ""
  )

  if (length(pick_obj$selected_terms) > 0) {

    cat(" - Selected terms:\n")

    cat(
      "   ",
      paste(pick_obj$selected_terms, collapse = ", "),
      "\n",
      sep = ""
    )

  } else {

    cat(" - Selected terms: NONE\n")
  }
}

# ============================================================
# Run models
# ============================================================
# For each outcome:
# - check available predictors
# - define complete-case dataset for that predictor set
# - build one common outer fold assignment
# - fit all alpha values
# - choose the best alpha
# - store both the full summary and the chosen object
# ============================================================

all_rows <- list()
chosen_rows <- list()
enet_objects <- list()

missing_log_env <- new.env(
  parent = emptyenv()
)

missing_log_env$seen <- list()

for (outcome_var in names(outcome_map)) {

  if (!(outcome_var %in% names(analysis_df))) {

    cat(
      outcome_var,
      " | SKIP: outcome not found in analysis_df.\n",
      sep = ""
    )

    next
  }

  preds <- outcome_map[[outcome_var]]

  missing_pred0 <- setdiff(
    unique(preds),
    names(analysis_df)
  )

  if (length(missing_pred0) > 0 &&
      !isTRUE(missing_log_env$seen[[outcome_var]])) {

    cat(
      outcome_var,
      " | Missing predictors dropped: ",
      paste(sort(unique(missing_pred0)), collapse = ", "),
      "\n",
      sep = ""
    )

    missing_log_env$seen[[outcome_var]] <- TRUE
  }

  preds_in <- unique(preds)

  preds_in <- setdiff(
    preds_in,
    drop_vars
  )

  preds_in <- intersect(
    preds_in,
    names(analysis_df)
  )

  if (length(preds_in) == 0) {

    cat(
      outcome_var,
      " | SKIP: no predictors available after filtering.\n",
      sep = ""
    )

    next
  }

  df_tmp <- analysis_df %>%
    select(all_of(c(outcome_var, preds_in))) %>%
    drop_na()

  if (nrow(df_tmp) == 0) {

    cat(
      outcome_var,
      " | SKIP: no complete-case rows.\n",
      sep = ""
    )

    next
  }

  if (isTRUE(use_poly_for_ordinals)) {

    df_tmp <- apply_poly_contrasts(
      df_tmp,
      preds_in,
      level_order_3 = ordinal_level_order_3
    )
  }

  y_tmp <- to01(
    df_tmp[[outcome_var]]
  )

  if (length(unique(y_tmp)) < 2) {

    cat(
      outcome_var,
      " | SKIP: outcome has no variation after complete-case filtering.\n",
      sep = ""
    )

    next
  }

  foldid_outer <- make_strat_folds(
    y_tmp,
    K = nfolds_outer
  )

  res_list <- lapply(alpha_grid, function(a) {

    tryCatch(

      fit_enet_one_alpha(
        data = analysis_df,
        outcome_var = outcome_var,
        predictors = preds,
        alpha = a,
        foldid_outer = foldid_outer,
        nfolds_inner = nfolds_inner,
        use_poly_for_ordinals = use_poly_for_ordinals,
        ordinal_level_order_3 = ordinal_level_order_3,
        drop_vars = drop_vars
      ),

      error = function(e) {

        cat(
          outcome_var,
          " | alpha=", a,
          " | ERROR: ", conditionMessage(e),
          "\n",
          sep = ""
        )

        NULL
      }
    )
  })

  res_list <- Filter(
    Negate(is.null),
    res_list
  )

  if (length(res_list) == 0) next

  df_sum <- bind_rows(
    lapply(res_list, function(x) {

      data.frame(
        outcome = x$outcome,
        alpha = x$alpha,
        outer_folds = nfolds_outer,
        inner_folds = nfolds_inner,
        n_predictors = length(x$predictors_used),
        N_total = x$n_total,
        N_cc = x$n_cc,
        events = x$events,
        nonevents = x$nonevents,
        lambda_1se = x$lambda_1se,
        n_selected = x$n_selected,
        AUC = x$AUC,
        AUC_low = x$AUC_low,
        AUC_high = x$AUC_high,
        stringsAsFactors = FALSE
      )
    })
  )

  for (i in seq_len(nrow(df_sum))) {

    cat(
      df_sum$outcome[i],
      " | alpha=", df_sum$alpha[i],
      " | N_cc=", df_sum$N_cc[i],
      " | events=", df_sum$events[i],
      " | AUC(outer)=", sprintf("%.3f", df_sum$AUC[i]),
      " | n_sel=", df_sum$n_selected[i],
      "\n",
      sep = ""
    )
  }

  cat("\n")

  chosen <- choose_alpha(
    df_sum,
    auc_tol = auc_tol
  )

  pick_obj <- res_list[[
    which.min(
      abs(
        sapply(res_list, `[[`, "alpha") -
          chosen$alpha
      )
    )
  ]]

  print_chosen_model(
    pick_obj,
    chosen
  )

  cat("\n")

  all_rows[[outcome_var]] <- df_sum
  chosen_rows[[outcome_var]] <- chosen
  enet_objects[[outcome_var]] <- pick_obj
}

summary_all <- bind_rows(all_rows)

cat("\nSUMMARY (all alphas) — AUC is outer nested CV\n")
print(summary_all)

cat("\nCHOSEN (one model per outcome)\n")

chosen_df <- bind_rows(chosen_rows)
print(chosen_df)

# ============================================================
# Stored objects
# ============================================================
# enet_objects contains the chosen model for each outcome, including:
# - out-of-fold probabilities from the outer CV
# - observed outcome vector
# - final glmnet fit
#
# Example:
# obj <- enet_objects[["OUTCOME_1"]]
# roc_obj <- pROC::roc(obj$y, obj$prob_oof, quiet = TRUE)
# pROC::auc(roc_obj)
# ============================================================

# ============================================================
# Calibration plots from out-of-fold predictions
# ============================================================
# Purpose:
# - assess agreement between predicted probabilities and observed events
# - group participants into five quantile-based bins of predicted risk
# - calculate mean predicted probability and observed event rate
# - visualise grouped calibration against the identity line
# - save one calibration plot per outcome
#
# Notes:
# - calibration uses the outer-CV out-of-fold predictions
# - the same five-bin procedure is applied across outcomes
# - plots provide a descriptive assessment of grouped calibration
# - calibration intercepts and slopes are not estimated
# - no decision-curve analysis is performed
# ============================================================

# Create the output directory if necessary
dir.create(
  calibration_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ============================================================
# Calibration plot function
# ============================================================
# Input:
# - obj: selected model stored in enet_objects
# - outcome_name: name of the outcome
# - output_dir: directory for PDF output
# - n_bins: number of probability groups
#
# Output:
# - grouped calibration data
# - ggplot object
# - PDF calibration plot
# ============================================================

make_calibration_plot <- function(
    obj,
    outcome_name,
    output_dir = "Outputs",
    n_bins = 5
) {

  stopifnot(
    length(obj$y) == length(obj$prob_oof)
  )

  # Prepare observed outcomes and OOF predicted probabilities
  cal_data <- tibble(
    observed = as.numeric(obj$y),
    predicted = as.numeric(obj$prob_oof)
  ) %>%
    filter(
      is.finite(predicted),
      !is.na(observed)
    )

  if (nrow(cal_data) == 0) {

    warning(
      "No valid OOF predictions for outcome: ",
      outcome_name
    )

    return(NULL)
  }

  stopifnot(
    all(cal_data$observed %in% c(0, 1)),
    all(cal_data$predicted >= 0 & cal_data$predicted <= 1)
  )

  # Divide participants into quantile-based groups
  # and calculate grouped calibration statistics
  cal_bins <- cal_data %>%
    mutate(
      bin = ntile(
        predicted,
        min(n_bins, n())
      )
    ) %>%
    group_by(bin) %>%
    summarise(
      n = n(),
      mean_predicted = mean(predicted),
      observed_rate = mean(observed),
      .groups = "drop"
    )

  # Construct calibration plot
  p <- ggplot(
    cal_bins,
    aes(
      x = mean_predicted,
      y = observed_rate
    )
  ) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      color = "grey50"
    ) +
    geom_line() +
    geom_point(
      aes(size = n)
    ) +
    coord_equal(
      xlim = c(0, 1),
      ylim = c(0, 1)
    ) +
    labs(
      x = "Mean predicted probability",
      y = "Observed event rate",
      size = "Subjects"
    ) +
    scale_size_continuous(
      breaks = sort(unique(cal_bins$n)),
      labels = scales::label_number(accuracy = 1)
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 13),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    )

  # Create a filename from the outcome name
  file_label <- gsub(
    "[^A-Za-z0-9_-]",
    "_",
    outcome_name
  )

  output_file <- file.path(
    output_dir,
    paste0(
      "calibration_",
      file_label,
      "_",
      n_bins,
      "_bins.pdf"
    )
  )

  # Save PDF
  ggsave(
    filename = output_file,
    plot = p,
    width = 7,
    height = 7,
    units = "in"
  )

  cat(
    "Calibration plot saved: ",
    output_file,
    "\n",
    sep = ""
  )

  # Return the plot and underlying calibration data
  list(
    bins = cal_bins,
    plot = p,
    file = output_file
  )
}

# ============================================================
# Generate calibration plots for all selected models
# ============================================================
# The procedure is applied automatically to every outcome
# successfully fitted and stored in enet_objects.
#
# Example output:
# Outputs/calibration_Remission_5_bins.pdf
# Outputs/calibration_TR1_5_bins.pdf
# Outputs/calibration_TR2_5_bins.pdf
# ============================================================

calibration_results <- lapply(
  names(enet_objects),
  function(outcome_var) {

    make_calibration_plot(
      obj = enet_objects[[outcome_var]],
      outcome_name = outcome_var,
      output_dir = calibration_dir,
      n_bins = n_calibration_bins
    )
  }
)

names(calibration_results) <- names(enet_objects)

# ============================================================
# Calibration summary
# ============================================================
# calibration_results contains, for each outcome:
# - grouped predicted probabilities
# - observed event rates
# - number of participants per group
# - ggplot object and saved PDF path
# ============================================================

cat("\nCALIBRATION PLOTS\n")

for (outcome_var in names(calibration_results)) {

  result <- calibration_results[[outcome_var]]

  if (is.null(result)) next

  cat(
    "\nOutcome: ",
    outcome_var,
    "\n",
    sep = ""
  )

  print(result$bins)

  cat(
    "Saved to: ",
    result$file,
    "\n",
    sep = ""
  )
}

# ============================================================
# End of multivariable modelling and calibration workflow
# ============================================================
