############################################################
# Project: Transdiagnostic predictors of treatment outcomes
# Script: Preprocessing workflow for cohort-level analytic datasets
#
# Purpose:
# - assemble one baseline analytic dataset from cohort-specific raw files
# - derive harmonised sociodemographic and clinical predictors
# - merge polygenic scores and ancestry principal components
# - residualise polygenic scores on ancestry principal components
# - derive treatment outcome phenotypes when available
#
# Notes:
# - this script is intentionally generic and should be adapted to the
#   structure of each cohort before use
# - file names, visit labels, ID variables, and coding schemes may differ
#   across studies
# - the aim is not to force identical preprocessing across cohorts, but to
#   keep a common workflow and naming logic where possible
############################################################

library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(fs)

options(stringsAsFactors = FALSE)

# ============================================================
# Input paths
# ============================================================

raw_data_dir  <- "INSERT_PATH_TO_RAW_DATA"
prs_dir       <- "INSERT_PATH_TO_PRS"
pheno_dir     <- "INSERT_PATH_TO_PHENOTYPES"

# ============================================================
# Main dataset
# ============================================================

cohort_df <- read.csv("COHORT.csv", check.names = FALSE, stringsAsFactors = FALSE)

# ============================================================
# Helper functions
# ============================================================

canon_id <- function(x) {
  gsub("[^A-Z0-9]", "", toupper(as.character(x)))
}

parse_num <- function(x) {
  suppressWarnings(readr::parse_number(as.character(x)))
}

recode_yes_no <- function(x) {
  z <- trimws(as.character(x))
  case_when(
    z %in% c("1", "Yes", "YES", "yes", "Y", "y") ~ 1,
    z %in% c("0", "No", "NO", "no", "N", "n")    ~ 0,
    TRUE                                         ~ NA_real_
  )
}

read_ws <- function(path) {
  x <- readr::read_table(path, col_names = TRUE, progress = FALSE)
  x <- x %>% mutate(across(where(is.character), ~ str_squish(.)))
  if ("type_convert" %in% getNamespaceExports("readr")) {
    x <- tryCatch(readr::type_convert(x), error = function(e) x)
  }
  x
}

pick_baseline <- function(df) {
  if ("BASE1" %in% names(df)) {
    df %>% filter(BASE1 == 1)
  } else if ("VISIT" %in% names(df)) {
    df %>% filter(grepl("^Baseline$|^Entry$", VISIT, ignore.case = TRUE))
  } else if (all(c("SUBJECT_ID", "VISITID") %in% names(df))) {
    df %>%
      group_by(SUBJECT_ID) %>%
      slice_min(VISITID, with_ties = FALSE) %>%
      ungroup()
  } else if (all(c("ID", "VISITID") %in% names(df))) {
    df %>%
      group_by(ID) %>%
      slice_min(VISITID, with_ties = FALSE) %>%
      ungroup()
  } else {
    df
  }
}

# Standardise a PC table so that:
# - first column is the generic genetic ID
# - remaining columns are PC1, PC2, ...
# - duplicated ID columns are removed if present
normalize_pc_df <- function(x) {
  if (ncol(x) == 0) return(tibble())
  
  x <- x %>% mutate(across(where(is.character), ~ str_squish(.)))
  
  if (ncol(x) >= 2 && isTRUE(all(x[[1]] == x[[2]], na.rm = TRUE))) {
    x <- x[, -2, drop = FALSE]
  }
  
  names(x)[1] <- "GEN_ID"
  
  if (ncol(x) > 1) {
    names(x)[-1] <- paste0("PC", seq_len(ncol(x) - 1))
  }
  
  x <- x %>%
    mutate(across(starts_with("PC"), ~ suppressWarnings(as.numeric(.))))
  
  x$GEN_ID <- ifelse(
    grepl("_.+", x$GEN_ID),
    sub("_.+$", "", x$GEN_ID),
    x$GEN_ID
  )
  
  distinct(x, GEN_ID, .keep_all = TRUE)
}

read_pc_file <- function(path) {
  x <- tryCatch(
    readr::read_table(path, col_names = TRUE, progress = FALSE),
    error = function(e) NULL
  )
  
  bad_header <- !is.null(x) &&
    (any(grepl("^-?\\d", names(x))) || any(grepl("^[0-9]{2}-", names(x))))
  
  if (is.null(x) || bad_header) {
    x <- readr::read_table(path, col_names = FALSE, progress = FALSE)
  }
  
  normalize_pc_df(x)
}

read_cohort_pcs <- function(path = file.path(pheno_dir, "COHORT_PCs")) {
  if (fs::file_exists(path) && !fs::is_dir(path)) {
    return(read_pc_file(path))
  }
  
  if (fs::is_dir(path)) {
    files <- dir_ls(path, recurse = TRUE, type = "file")
    if (!length(files)) return(tibble())
    pcs_list <- lapply(files, read_pc_file)
    return(bind_rows(pcs_list))
  }
  
  tibble()
}

residualize_prs <- function(df, prs_name, pcs = c("PC1", "PC2", "PC3", "PC4")) {
  vars_needed <- c(prs_name, pcs)
  vars_needed <- vars_needed[vars_needed %in% names(df)]
  
  if (!prs_name %in% vars_needed) {
    df[[paste0(prs_name, "_resid")]]   <- NA_real_
    df[[paste0(prs_name, "_resid_z")]] <- NA_real_
    return(df)
  }
  
  tmp <- df[, vars_needed, drop = FALSE]
  ok  <- complete.cases(tmp)
  
  if (sum(ok) < 10) {
    df[[paste0(prs_name, "_resid")]]   <- NA_real_
    df[[paste0(prs_name, "_resid_z")]] <- NA_real_
    return(df)
  }
  
  pc_terms <- setdiff(vars_needed, prs_name)
  
  if (length(pc_terms) == 0) {
    resid_vals <- scale(tmp[[prs_name]][ok])
  } else {
    form <- as.formula(paste(prs_name, "~", paste(pc_terms, collapse = " + ")))
    fit  <- lm(form, data = tmp[ok, , drop = FALSE])
    resid_vals <- resid(fit)
  }
  
  df[[paste0(prs_name, "_resid")]]   <- NA_real_
  df[[paste0(prs_name, "_resid_z")]] <- NA_real_
  
  idx <- which(ok)
  df[[paste0(prs_name, "_resid")]][idx]   <- as.numeric(resid_vals)
  df[[paste0(prs_name, "_resid_z")]][idx] <- as.numeric(scale(resid_vals))
  
  df
}

# ============================================================
# Cohort-specific raw sources
# ============================================================

demographics_df <- read.csv(file.path(raw_data_dir, "DEMOGRAPHICS.csv"),
                            check.names = FALSE, stringsAsFactors = FALSE)

clinical_df <- read.csv(file.path(raw_data_dir, "CLINICAL.csv"),
                        check.names = FALSE, stringsAsFactors = FALSE)

function_df <- read.csv(file.path(raw_data_dir, "FUNCTIONING.csv"),
                        check.names = FALSE, stringsAsFactors = FALSE)

history_df <- read.csv(file.path(raw_data_dir, "HISTORY.csv"),
                       check.names = FALSE, stringsAsFactors = FALSE)

# ============================================================
# Baseline records
# ============================================================

clinical_bl <- pick_baseline(clinical_df)
function_bl <- pick_baseline(function_df)
history_bl  <- pick_baseline(history_df)

# ============================================================
# Core demographic variables
# ============================================================

cohort_df <- cohort_df %>%
  left_join(
    demographics_df %>%
      transmute(
        SUBJECT_ID     = as.character(SUBJECT_ID),
        Age_raw        = AGE,
        Sex_raw        = SEX,
        Employment_raw = EMPLOYMENT,
        Marital_raw    = MARITAL_STATUS,
        Education_raw  = EDUCATION,
        Income_raw     = INCOME,
        LivesAlone_raw = LIVE_ALONE,
        Race_raw       = RACE,
        Ethnicity_raw  = ETHNICITY
      ),
    by = c("ID" = "SUBJECT_ID")
  ) %>%
  mutate(
    Age = suppressWarnings(as.numeric(Age_raw)),
    
    Sex = case_when(
      Sex_raw %in% c("Male", "M", 1)   ~ "Male",
      Sex_raw %in% c("Female", "F", 2) ~ "Female",
      TRUE                             ~ NA_character_
    ),
    
    Employment = case_when(
      str_to_lower(trimws(as.character(Employment_raw))) %in%
        c("full time", "part time", "full-time", "part-time", "employed") ~ 1L,
      str_to_lower(trimws(as.character(Employment_raw))) %in%
        c("unemployed", "did not work", "not working") ~ 0L,
      TRUE ~ NA_integer_
    ),
    
    Married = case_when(
      str_to_lower(trimws(as.character(Marital_raw))) %in%
        c("married", "living as married", "cohabiting") ~ 1L,
      is.na(Marital_raw) ~ NA_integer_,
      TRUE               ~ 0L
    ),
    
    Lives_Alone = recode_yes_no(LivesAlone_raw),
    Education_years = suppressWarnings(as.numeric(Education_raw))
  ) %>%
  mutate(
    Education_3cat = case_when(
      Education_years <= 12 ~ 1L,
      Education_years <= 16 ~ 2L,
      Education_years > 16  ~ 3L,
      TRUE                  ~ NA_integer_
    ),
    Education_3cat_ord = factor(
      Education_3cat,
      levels = c(1, 2, 3),
      labels = c("Low", "Moderate", "High"),
      ordered = TRUE
    )
  )

# ============================================================
# Clinical severity and symptom dimensions
# ============================================================

cohort_df <- cohort_df %>%
  left_join(
    clinical_bl %>%
      transmute(
        SUBJECT_ID             = as.character(SUBJECT_ID),
        Severity_raw           = GLOBAL_SEVERITY,
        Anxiety_raw            = ANXIETY,
        Psychotic_raw          = PSYCHOTIC_SYMPTOMS,
        Suicidality_raw        = SUICIDALITY,
        SocialWithdrawal_raw   = SOCIAL_WITHDRAWAL,
        DepressiveMood_raw     = DEPRESSED_MOOD,
        PreviousEpisodes_raw   = PREVIOUS_EPISODES,
        Hospitalization_raw    = HOSPITALIZATION,
        Comorbidity_raw        = COMORBIDITY
      ),
    by = c("ID" = "SUBJECT_ID")
  ) %>%
  mutate(
    Severity = parse_num(Severity_raw),
    Anxiety  = parse_num(Anxiety_raw),
    Psychotic_Symptoms = parse_num(Psychotic_raw),
    Suicidality = parse_num(Suicidality_raw),
    Social_Withdrawal = parse_num(SocialWithdrawal_raw),
    Depressed_Mood = parse_num(DepressiveMood_raw),
    N_previous_episodes = parse_num(PreviousEpisodes_raw),
    Hospitalization = recode_yes_no(Hospitalization_raw),
    Comorbidity = recode_yes_no(Comorbidity_raw)
  ) %>%
  mutate(
    Severity_3cat = case_when(
      Severity %in% c(1, 2)    ~ 1L,
      Severity %in% c(3, 4)    ~ 2L,
      Severity %in% c(5, 6, 7) ~ 3L,
      TRUE                     ~ NA_integer_
    ),
    Severity_3cat_ord = factor(
      Severity_3cat,
      levels = c(1, 2, 3),
      labels = c("Low", "Moderate", "High"),
      ordered = TRUE
    ),
    
    Social_Withdrawal_3cat = case_when(
      is.na(Social_Withdrawal)   ~ NA_integer_,
      Social_Withdrawal %in% 0:2 ~ 1L,
      Social_Withdrawal == 3     ~ 2L,
      Social_Withdrawal >= 4     ~ 3L,
      TRUE                       ~ NA_integer_
    ),
    Social_Withdrawal_3cat_ord = factor(
      Social_Withdrawal_3cat,
      levels = c(1, 2, 3),
      labels = c("Low", "Moderate", "High"),
      ordered = TRUE
    ),
    
    N_prev_episodes_3cat = case_when(
      N_previous_episodes %in% 0:1 ~ 1L,
      N_previous_episodes %in% 2:3 ~ 2L,
      N_previous_episodes >= 4     ~ 3L,
      TRUE                         ~ NA_integer_
    ),
    N_prev_episodes_3cat_ord = factor(
      N_prev_episodes_3cat,
      levels = c(1, 2, 3),
      labels = c("Low", "Moderate", "High"),
      ordered = TRUE
    )
  )

# ============================================================
# Functioning
# ============================================================

cohort_df <- cohort_df %>%
  left_join(
    function_bl %>%
      transmute(
        SUBJECT_ID = as.character(SUBJECT_ID),
        Functioning = suppressWarnings(as.numeric(FUNCTIONING_TOTAL))
      ),
    by = c("ID" = "SUBJECT_ID")
  )

# ============================================================
# Optional cohort-specific variables
# ============================================================

# Example:
# cohort_df <- cohort_df %>%
#   mutate(
#     Age_at_onset = suppressWarnings(as.numeric(AGE_AT_ONSET)),
#     Family_History = recode_yes_no(FAMILY_HISTORY)
#   )

# ============================================================
# Polygenic risk scores
# ============================================================
# GEN_ID is used as a generic genetic identifier.
# If the source PRS file uses a cohort-specific name such as GEN_PUBLICID,
# rename it to GEN_ID immediately after reading.
# ============================================================

prs_df <- read_table(
  file.path(prs_dir, "all_COHORT_PRSs_standard"),
  col_types = cols(
    IID        = col_character(),
    GEN_ID     = col_character(),
    SUBJECT_ID = col_character(),
    .default   = col_double()
  )
)

# Optional compatibility step for datasets that still use GEN_PUBLICID
if ("GEN_PUBLICID" %in% names(prs_df) && !"GEN_ID" %in% names(prs_df)) {
  prs_df <- prs_df %>% rename(GEN_ID = GEN_PUBLICID)
}

prs_vars <- c(
  "ADHD05_PRS", "BODY04_PRS", "BP_PRS",
  "DEPR07_PRS", "EDUC03_PRS", "SCZ_PRS"
)

prs_vars_present <- prs_vars[prs_vars %in% names(prs_df)]

prs_keep <- prs_df %>%
  mutate(
    SUBJECT_ID = as.character(SUBJECT_ID),
    IID_key    = canon_id(sub("_[^_]+$", "", IID)),
    GEN_key    = canon_id(GEN_ID),
    SUBJ_key   = canon_id(SUBJECT_ID)
  )

cohort_df <- cohort_df %>%
  mutate(ID_key = canon_id(ID)) %>%
  left_join(
    prs_keep %>%
      transmute(
        ID_key = coalesce(SUBJ_key, IID_key, GEN_key),
        across(all_of(c("IID", "GEN_ID", prs_vars_present)))
      ) %>%
      distinct(ID_key, .keep_all = TRUE),
    by = "ID_key"
  ) %>%
  select(-ID_key)

# ============================================================
# Principal components
# ============================================================
# PCs are also standardised to use GEN_ID as the genetic identifier.
# The code tries to match PCs using either IID or GEN_ID.
# ============================================================

pcs <- read_cohort_pcs(file.path(pheno_dir, "COHORT_PCs"))

pcs_join <- pcs %>%
  mutate(PC_KEY = canon_id(GEN_ID)) %>%
  select(PC_KEY, starts_with("PC")) %>%
  distinct(PC_KEY, .keep_all = TRUE)

cohort_df <- cohort_df %>%
  mutate(
    K_IID = if ("IID" %in% names(.)) {
      canon_id(sub("_[^_]+$", "", IID))
    } else {
      NA_character_
    },
    K_GEN = if ("GEN_ID" %in% names(.)) {
      canon_id(GEN_ID)
    } else {
      NA_character_
    }
  )

n_match_iid <- inner_join(
  cohort_df %>% filter(!is.na(K_IID)) %>% select(K_IID),
  pcs_join %>% select(PC_KEY),
  by = c("K_IID" = "PC_KEY")
) %>% nrow()

n_match_gen <- inner_join(
  cohort_df %>% filter(!is.na(K_GEN)) %>% select(K_GEN),
  pcs_join %>% select(PC_KEY),
  by = c("K_GEN" = "PC_KEY")
) %>% nrow()

best_key <- if (n_match_iid >= n_match_gen) "K_IID" else "K_GEN"

pc_cols_exist <- grep("^PC\\d+$", names(cohort_df), value = TRUE)

cohort_df <- cohort_df %>%
  select(-any_of(pc_cols_exist)) %>%
  left_join(pcs_join, by = setNames("PC_KEY", best_key)) %>%
  select(-K_IID, -K_GEN)

# ============================================================
# PRS residualisation
# ============================================================
# Residualise PRSs on ancestry PCs only once, after both PRS and
# population PC data have been merged into the analytic dataset.
# ============================================================

pc_vars <- c("PC1", "PC2", "PC3", "PC4")

for (prs in prs_vars_present) {
  cohort_df <- residualize_prs(cohort_df, prs_name = prs, pcs = pc_vars)
}

# ============================================================
# Treatment outcome phenotypes
# ============================================================

outcome_df <- read_ws(file.path(pheno_dir, "COHORT_OUTCOMES.txt"))

cohort_df <- cohort_df %>%
  mutate(ID_chr = as.character(ID)) %>%
  left_join(
    outcome_df %>%
      mutate(SUBJECT_ID = as.character(SUBJECT_ID)) %>%
      select(any_of(c("SUBJECT_ID", "Remission", "TR1"))),
    by = c("ID_chr" = "SUBJECT_ID")
  ) %>%
  select(-ID_chr)

# ---------- Example: derive TR2 from repeated severity scores ----------
# Replace cgi_long and CGI_raw with the actual longitudinal file and variable.

# cgi_long <- read.csv(file.path(raw_data_dir, "CGI_LONG.csv"),
#                      check.names = FALSE, stringsAsFactors = FALSE) %>%
#   mutate(
#     SUBJECT_ID = as.character(SUBJECT_ID),
#     CGI_num = parse_num(CGI_raw)
#   )
#
# tr2_df <- cgi_long %>%
#   group_by(SUBJECT_ID) %>%
#   summarise(
#     n_vis = sum(!is.na(CGI_num)),
#     all_gt3 = ifelse(n_vis >= 2, all(CGI_num > 3, na.rm = TRUE), NA),
#     TR2 = case_when(
#       is.na(all_gt3) ~ NA_integer_,
#       all_gt3        ~ 1L,
#       TRUE           ~ 0L
#     ),
#     .groups = "drop"
#   )
#
# cohort_df <- cohort_df %>%
#   mutate(ID_chr = as.character(ID)) %>%
#   left_join(tr2_df, by = c("ID_chr" = "SUBJECT_ID")) %>%
#   select(-ID_chr)

# ============================================================
# Final checks and light recoding
# ============================================================

cohort_df <- cohort_df %>%
  distinct(ID, .keep_all = TRUE)

# ============================================================
# Save output
# ============================================================

write.csv(cohort_df, "cohort_preprocessed.csv", row.names = FALSE, na = "")
