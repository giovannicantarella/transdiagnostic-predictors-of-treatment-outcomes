# ============================================================
# GSRD | Elastic Net alpha-tuning with nested cross-validation
# - Predictors included here were preselected based on univariable associations
#   at p < 0.10
#   were carried forward to these multivariable models
# - OUTER CV: honest out-of-fold AUC for alpha selection
# - INNER CV: lambda.1se selected within each outer training fold
# - Gelman rescaling: center and scale continuous columns by 2*SD only
# - Ordinal predictors are handled with polynomial contrasts when requested
# - Missing predictors are reported once per outcome
# - Selection rule:
#     1) prefer models with n_selected > 0
#     2) highest outer AUC
#     3) within auc_tol, choose smallest n_selected
#     4) tie-break: higher AUC, then lower alpha
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(glmnet)
  library(pROC)
})

set.seed(1234)
options(contrasts = c("contr.treatment", "contr.poly"))

# -------------------------
# Precheck
# -------------------------
if (!exists("GSRD")) {
  stop("Object 'GSRD' not found in the environment. Load the dataset before running this script.")
}

# -------------------------
# Settings
# -------------------------
alpha_grid <- seq(0.10, 0.30, by = 0.05)

nfolds_outer <- 5
nfolds_inner <- 5

auc_tol <- 0.01
drop_vars <- character(0)

use_poly_for_ordinals <- TRUE
ordinal_level_order_3 <- c("Low", "Moderate", "High")

# -------------------------
# Predictors by outcome
# -------------------------
pred_vars_remission <- c(
  "MADRS.Retrospective",
  "Education_3cat_ord",
  "Age.at.Onset_clean",
  "SuicidalRisk_bin",
  "EDUC03_resid_z",
  "SCZ_resid_z",
  "Age",
  "Sex"
)

pred_vars_TR <- c(
  "MADRS.Retrospective",
  "SuicidalRisk_bin",
  "Age",
  "Sex"
)

outcome_map <- list(
  remission  = pred_vars_remission,
  Resistance = pred_vars_TR
)

outcome_map <- lapply(outcome_map, function(v) setdiff(unique(v), drop_vars))

cat("GSRD elastic net alpha-tuning (nested CV + Gelman rescaling, lambda.1se)\n")
cat("Alpha grid:", paste(alpha_grid, collapse = ", "), "\n")
cat("Outer folds:", nfolds_outer, " | Inner folds:", nfolds_inner, " | auc_tol:", auc_tol, "\n")
cat("Ordinal polynomial contrasts:", use_poly_for_ordinals, "\n\n")

# -------------------------
# Helpers
# -------------------------
to01 <- function(x) {
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x)) return(as.integer(x))
  x <- as.character(x)
  as.integer(x %in% c("1", "TRUE", "T", "true", "Yes", "YES", "case", "Case"))
}

make_strat_folds <- function(y, K) {
  y <- as.integer(y)
  if (length(unique(y)) < 2) stop("make_strat_folds: y must contain both classes.")
  
  idx1 <- which(y == 1)
  idx0 <- which(y == 0)
  
  foldid <- integer(length(y))
  foldid[idx1] <- sample(rep(seq_len(K), length.out = length(idx1)))
  foldid[idx0] <- sample(rep(seq_len(K), length.out = length(idx0)))
  
  foldid
}

choose_alpha <- function(df_summary, auc_tol = 0.01) {
  df_summary <- df_summary %>%
    mutate(
      n_selected = as.integer(n_selected),
      AUC = as.numeric(AUC)
    )
  
  df_pos <- df_summary %>% filter(n_selected > 0)
  
  if (nrow(df_pos) == 0) {
    cat("[WARN] All models have n_selected == 0. Falling back to highest AUC.\n")
    auc_best <- max(df_summary$AUC, na.rm = TRUE)
    keep <- df_summary %>% filter(AUC >= (auc_best - auc_tol))
    return(keep %>% arrange(desc(AUC), alpha) %>% slice(1))
  }
  
  auc_best <- max(df_pos$AUC, na.rm = TRUE)
  keep <- df_pos %>% filter(AUC >= (auc_best - auc_tol))
  
  keep %>% arrange(n_selected, desc(AUC), alpha) %>% slice(1)
}

apply_poly_contrasts <- function(df, predictors_in,
                                 level_order_3 = c("Low", "Moderate", "High")) {
  
  ord_vars <- intersect(grep("_3cat_ord$", predictors_in, value = TRUE), names(df))
  if (length(ord_vars) == 0) return(df)
  
  for (v in ord_vars) {
    if (is.numeric(df[[v]])) next
    
    x <- as.character(df[[v]])
    
    if (!is.null(level_order_3) && length(level_order_3) == 3) {
      vals <- unique(na.omit(x))
      if (length(vals) > 0 && all(vals %in% level_order_3)) {
        df[[v]] <- ordered(x, levels = level_order_3)
      } else {
        df[[v]] <- ordered(df[[v]])
      }
    } else {
      df[[v]] <- ordered(df[[v]])
    }
    
    k <- length(levels(df[[v]]))
    if (k >= 2) contrasts(df[[v]]) <- contr.poly(k)
  }
  
  df
}

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

gelman_fit_cont <- function(X, cont_idx) {
  X <- as.matrix(X)
  
  mu <- rep(0, ncol(X))
  sc <- rep(1, ncol(X))
  
  if (length(cont_idx) > 0) {
    mu[cont_idx] <- colMeans(X[, cont_idx, drop = FALSE], na.rm = TRUE)
    sdv <- apply(X[, cont_idx, drop = FALSE], 2, sd, na.rm = TRUE)
    sc[cont_idx] <- 2 * sdv
    sc[cont_idx][!is.finite(sc[cont_idx]) | sc[cont_idx] == 0] <- 1
  }
  
  list(center = mu, scale = sc)
}

gelman_apply <- function(X, sc) {
  X <- as.matrix(X)
  sweep(sweep(X, 2, sc$center, "-"), 2, sc$scale, "/")
}

nested_auc_alpha <- function(X_raw, y, alpha, foldid_outer, nfolds_inner, cont_idx) {
  n <- nrow(X_raw)
  stopifnot(length(y) == n, length(foldid_outer) == n)
  
  prob_oof <- rep(NA_real_, n)
  
  for (k in sort(unique(foldid_outer))) {
    idx_te <- which(foldid_outer == k)
    idx_tr <- which(foldid_outer != k)
    
    Xtr_raw <- X_raw[idx_tr, , drop = FALSE]
    Xte_raw <- X_raw[idx_te, , drop = FALSE]
    ytr <- y[idx_tr]
    
    sc <- gelman_fit_cont(Xtr_raw, cont_idx)
    Xtr <- gelman_apply(Xtr_raw, sc)
    Xte <- gelman_apply(Xte_raw, sc)
    
    storage.mode(Xtr) <- "double"
    storage.mode(Xte) <- "double"
    
    foldid_inner <- make_strat_folds(ytr, K = nfolds_inner)
    
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
      predict(fit_in, newx = Xte, type = "response")
    )
  }
  
  ok <- is.finite(prob_oof)
  if (!all(ok)) {
    warning("OOF predictions contain NA/Inf. AUC will be computed on finite predictions only.")
  }
  
  roc_obj <- pROC::roc(y[ok], prob_oof[ok], quiet = TRUE)
  auc_val <- as.numeric(pROC::auc(roc_obj))
  auc_ci  <- as.numeric(pROC::ci.auc(roc_obj))
  
  list(
    AUC = auc_val,
    AUC_low = auc_ci[1],
    AUC_high = auc_ci[3],
    prob_oof = prob_oof
  )
}

fit_enet_one_alpha_cc <- function(df_cc, outcome_var, predictors_in, alpha, foldid_outer,
                                  nfolds_inner, use_poly_for_ordinals, ordinal_level_order_3) {
  
  if (isTRUE(use_poly_for_ordinals)) {
    df_cc <- apply_poly_contrasts(df_cc, predictors_in, level_order_3 = ordinal_level_order_3)
  }
  
  y <- to01(df_cc[[outcome_var]])
  if (length(unique(y)) < 2) {
    stop("Outcome has no variation after complete-case filtering.")
  }
  
  n_cc <- nrow(df_cc)
  n_events <- sum(y == 1)
  n_nonevents <- sum(y == 0)
  
  X_raw <- model.matrix(reformulate(predictors_in), data = df_cc)[, -1, drop = FALSE]
  cont_idx <- get_continuous_cols(X_raw)
  
  if (length(foldid_outer) != n_cc) {
    stop("foldid_outer length does not match the number of complete-case observations.")
  }
  
  nested <- nested_auc_alpha(
    X_raw = X_raw,
    y = y,
    alpha = alpha,
    foldid_outer = foldid_outer,
    nfolds_inner = nfolds_inner,
    cont_idx = cont_idx
  )
  
  sc_full <- gelman_fit_cont(X_raw, cont_idx)
  X <- gelman_apply(X_raw, sc_full)
  storage.mode(X) <- "double"
  
  foldid_inner_full <- make_strat_folds(y, K = nfolds_inner)
  
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
  selected <- setdiff(rownames(b)[as.numeric(b[, 1]) != 0], "(Intercept)")
  
  list(
    alpha = alpha,
    n_total = NA_integer_,
    n_cc = n_cc,
    events = n_events,
    nonevents = n_nonevents,
    lambda_min = cvfit$lambda.min,
    lambda_1se = cvfit$lambda.1se,
    n_selected = length(selected),
    selected_terms = selected,
    AUC = nested$AUC,
    AUC_low = nested$AUC_low,
    AUC_high = nested$AUC_high,
    prob_oof = nested$prob_oof,
    y = y,
    cvfit = cvfit,
    fit = fit,
    X = X,
    predictors_used = predictors_in
  )
}

print_chosen_model <- function(pick_obj, outcome_var, chosen_row) {
  
  cat(
    "\nCHOSEN | ", outcome_var,
    " | alpha=", sprintf("%.2f", chosen_row$alpha),
    " | AUC(outer)=", sprintf("%.3f", chosen_row$AUC),
    " [", sprintf("%.3f", chosen_row$AUC_low), ", ", sprintf("%.3f", chosen_row$AUC_high), "]",
    " | n_sel=", chosen_row$n_selected,
    " | lambda.1se(full fit)=", sprintf("%.6f", pick_obj$lambda_1se),
    "\n", sep = ""
  )
  
  b <- as.matrix(coef(pick_obj$fit))
  b_df <- data.frame(
    term = rownames(b),
    beta = as.numeric(b[, 1]),
    stringsAsFactors = FALSE
  )
  b_df <- b_df[b_df$beta != 0, , drop = FALSE]
  b_df <- b_df[order(-abs(b_df$beta)), , drop = FALSE]
  
  cat("Non-zero coefficients (sorted by absolute beta):\n")
  if (nrow(b_df) == 0) {
    cat("(NONE)\n")
  } else {
    print(b_df, row.names = FALSE)
  }
  
  cat("Summary:\n")
  cat(" - N_total: ", pick_obj$n_total, "\n", sep = "")
  cat(" - N_cc: ", pick_obj$n_cc, "\n", sep = "")
  cat(" - Events: ", pick_obj$events, " | Nonevents: ", pick_obj$nonevents, "\n", sep = "")
  cat(" - Predictors used: ", length(pick_obj$predictors_used), "\n", sep = "")
  cat(" - Design matrix columns: ", ncol(pick_obj$X), "\n", sep = "")
  cat(" - lambda.min: ", sprintf("%.6f", pick_obj$lambda_min), "\n", sep = "")
  cat(" - lambda.1se: ", sprintf("%.6f", pick_obj$lambda_1se), "\n", sep = "")
  cat(" - n_nonzero (excluding intercept): ", pick_obj$n_selected, "\n", sep = "")
  
  if (length(pick_obj$selected_terms) > 0) {
    cat(" - Selected terms:\n")
    cat("   ", paste(pick_obj$selected_terms, collapse = ", "), "\n", sep = "")
  } else {
    cat(" - Selected terms: NONE\n")
  }
}

# -------------------------
# Run models
# -------------------------
all_rows <- list()
chosen_rows <- list()
enet_objects <- list()

missing_log_env <- new.env(parent = emptyenv())
missing_log_env$seen <- list()

for (outcome_var in names(outcome_map)) {
  
  cat("\n=====================================================\n")
  cat("OUTCOME:", outcome_var, "\n")
  cat("=====================================================\n")
  
  if (!(outcome_var %in% names(GSRD))) {
    cat(outcome_var, " | SKIP: outcome not found in GSRD.\n", sep = "")
    next
  }
  
  preds <- outcome_map[[outcome_var]]
  
  miss_key <- outcome_var
  missing_pred0 <- setdiff(unique(preds), names(GSRD))
  if (length(missing_pred0) > 0 && !isTRUE(missing_log_env$seen[[miss_key]])) {
    cat(
      outcome_var, " | Missing predictors dropped: ",
      paste(sort(unique(missing_pred0)), collapse = ", "),
      "\n\n", sep = ""
    )
    missing_log_env$seen[[miss_key]] <- TRUE
  }
  
  predictors_in <- intersect(preds, names(GSRD))
  if (length(predictors_in) == 0) {
    cat(outcome_var, " | SKIP: no predictors available after filtering.\n", sep = "")
    next
  }
  
  df_raw <- GSRD %>%
    dplyr::select(dplyr::all_of(c(outcome_var, predictors_in)))
  
  n_total <- nrow(df_raw)
  
  df_cc <- df_raw %>% tidyr::drop_na()
  if (nrow(df_cc) == 0) {
    cat(outcome_var, " | SKIP: no complete-case rows.\n", sep = "")
    next
  }
  
  if (isTRUE(use_poly_for_ordinals)) {
    df_cc <- apply_poly_contrasts(df_cc, predictors_in, level_order_3 = ordinal_level_order_3)
  }
  
  y_tmp <- to01(df_cc[[outcome_var]])
  if (length(unique(y_tmp)) < 2) {
    cat(outcome_var, " | SKIP: outcome has no variation after complete-case filtering.\n", sep = "")
    next
  }
  
  foldid_outer <- make_strat_folds(y_tmp, K = nfolds_outer)
  
  res_list <- lapply(alpha_grid, function(a) {
    tryCatch(
      {
        obj <- fit_enet_one_alpha_cc(
          df_cc = df_cc,
          outcome_var = outcome_var,
          predictors_in = predictors_in,
          alpha = a,
          foldid_outer = foldid_outer,
          nfolds_inner = nfolds_inner,
          use_poly_for_ordinals = FALSE,
          ordinal_level_order_3 = ordinal_level_order_3
        )
        obj$n_total <- n_total
        obj
      },
      error = function(e) {
        cat(outcome_var, " | alpha=", a, " | ERROR: ", conditionMessage(e), "\n", sep = "")
        NULL
      }
    )
  })
  
  res_list <- Filter(Negate(is.null), res_list)
  if (length(res_list) == 0) next
  
  df_sum <- bind_rows(lapply(res_list, function(x) {
    data.frame(
      outcome = outcome_var,
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
  }))
  
  for (i in seq_len(nrow(df_sum))) {
    cat(
      df_sum$outcome[i],
      " | alpha=", df_sum$alpha[i],
      " | N_cc=", df_sum$N_cc[i],
      " | events=", df_sum$events[i],
      " | AUC(outer)=", sprintf("%.3f", df_sum$AUC[i]),
      " | n_sel=", df_sum$n_selected[i],
      "\n", sep = ""
    )
  }
  cat("\n")
  
  chosen <- choose_alpha(df_sum, auc_tol = auc_tol)
  pick_obj <- res_list[[which.min(abs(sapply(res_list, `[[`, "alpha") - chosen$alpha))]]
  pick_obj$n_total <- n_total
  
  print_chosen_model(
    pick_obj = pick_obj,
    outcome_var = outcome_var,
    chosen_row = chosen
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

# enet_objects contains the chosen models with out-of-fold probabilities and outcomes
# Example:
# obj <- enet_objects[["remission"]]
# roc_obj <- pROC::roc(obj$y, obj$prob_oof, quiet = TRUE)
# pROC::auc(roc_obj)