############################################################
# Project: Transdiagnostic predictors of treatment outcomes
# Cohort: STEP-BD
# Script: Data preprocessing
############################################################

library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(fs)

options(stringsAsFactors = FALSE)

# -----------------------------
# Input folders
# -----------------------------
# Replace these with the local folders containing the STEP-BD source files

clinical_dir <- "path/to/STEP-BD/Clinician Rated Forms"
prs_dir      <- "path/to/STEP-BD/PRSs"
pheno_dir    <- "path/to/STEP-BD/phenotypes"

# -----------------------------
# Main dataset
# -----------------------------
STEP <- read.csv(
  "STEP.csv",
  header = TRUE,
  sep = ",",
  na.strings = c("", "NA"),
  stringsAsFactors = FALSE
)

# -----------------------------
# Helpers
# -----------------------------
canon_id <- function(x) {
  gsub("[^A-Z0-9]", "", toupper(as.character(x)))
}

parse_num <- function(x) {
  v <- if (is.numeric(x)) x else readr::parse_number(as.character(x))
  v[v %in% c(-55555, -3121, -3107, -961, -9, -8, -7, -6, -5, -4,
             -66, -666, -6666, -66666, -44444)] <- NA
  v
}

map_yesno <- function(x) {
  z <- str_trim(as.character(x))
  case_when(
    str_detect(z, "^Yes(\\b|\\s|\\()") ~ 1,
    str_detect(z, "^No(\\b|\\s|\\()")  ~ 0,
    TRUE ~ NA_real_
  )
}

pick_col <- function(nms, exact = character(), regex = NULL) {
  hit <- intersect(exact, nms)
  if (length(hit)) return(hit[1])
  if (!is.null(regex)) {
    g <- grep(regex, nms, value = TRUE, ignore.case = TRUE, perl = TRUE)
    if (length(g)) return(g[1])
  }
  NA_character_
}

read_ws <- function(path) {
  x <- readr::read_table(path, col_names = TRUE, progress = FALSE)
  x <- x %>% mutate(across(where(is.character), ~ str_squish(.)))
  if ("type_convert" %in% getNamespaceExports("readr")) {
    x <- tryCatch(readr::type_convert(x), error = function(e) x)
  }
  x
}

read_step_pcs <- function(path = file.path(pheno_dir, "STEP_PCs")) {
  first <- readr::read_lines(path, n_max = 1)
  toks  <- unlist(strsplit(first, "\\s+"))
  has_header <- any(grepl("[A-Za-z]", toks[-1], perl = TRUE))
  
  if (has_header) {
    df <- readr::read_table(path, col_names = TRUE, progress = FALSE)
  } else {
    df <- readr::read_table(path, col_names = FALSE, progress = FALSE)
  }
  
  df <- df %>% mutate(across(where(is.character), ~ str_squish(.)))
  
  if (!has_header) {
    colnames(df)[1] <- "GEN_PUBLICID"
    if (ncol(df) > 1) colnames(df)[-1] <- paste0("PC", seq_len(ncol(df) - 1))
  }
  
  if (ncol(df) >= 2 && isTRUE(all(df[[1]] == df[[2]], na.rm = TRUE))) {
    df <- df[, -2, drop = FALSE]
  }
  
  if (colnames(df)[1] != "GEN_PUBLICID") {
    names(df)[1] <- "GEN_PUBLICID"
  }
  
  if (ncol(df) > 1) {
    names(df)[-1] <- paste0("PC", seq_len(ncol(df) - 1))
  }
  
  df %>%
    mutate(
      across(starts_with("PC"), ~ suppressWarnings(as.numeric(.))),
      ID_canon = canon_id(GEN_PUBLICID)
    )
}

residualize_prs <- function(df, prs_name, pcs = c("PC1", "PC2", "PC3", "PC4")) {
  vars_needed <- c(prs_name, pcs)
  tmp <- df[, vars_needed, drop = FALSE]
  ok <- complete.cases(tmp)
  
  if (sum(ok) < 10) {
    df[[paste0(prs_name, "_resid")]]   <- NA_real_
    df[[paste0(prs_name, "_resid_z")]] <- NA_real_
    return(df)
  }
  
  form <- as.formula(paste(prs_name, "~", paste(pcs, collapse = " + ")))
  fit  <- lm(form, data = tmp[ok, , drop = FALSE])
  
  resid_vals <- resid(fit)
  
  df[[paste0(prs_name, "_resid")]]   <- NA_real_
  df[[paste0(prs_name, "_resid_z")]] <- NA_real_
  
  idx <- which(ok)
  df[[paste0(prs_name, "_resid")]][idx]   <- resid_vals
  df[[paste0(prs_name, "_resid_z")]][idx] <- as.numeric(scale(resid_vals))
  
  df
}

# -----------------------------
# Demographic variables
# -----------------------------
DEMO <- read.csv(
  file.path(clinical_dir, "df_gen.csv"),
  header = TRUE,
  sep = ",",
  na.strings = c("", "NA"),
  stringsAsFactors = FALSE
)

vars_needed <- c(
  "GEN_PUBLICID", "VISIT", "DAYSCONS", "AGE", "GENDER", "EMPLOY",
  "MARSTAT", "EDUCATE", "MENDAYS", "LIVEALON", "HINCOME",
  "RACE", "ETHNICTY", "ETHNSPEC"
)

missing_vars <- setdiff(vars_needed, names(DEMO))
if ("MENDAYS" %in% missing_vars) {
  DEMO$MENDAYS <- NA_integer_
  missing_vars <- setdiff(missing_vars, "MENDAYS")
}

if (length(missing_vars)) {
  message("Variables not found in DEMO: ", paste(missing_vars, collapse = ", "))
}

DEMO_base <- DEMO %>%
  filter(VISIT == "Entry") %>%
  mutate(abs_day = ifelse(is.na(DAYSCONS), Inf, abs(DAYSCONS))) %>%
  group_by(GEN_PUBLICID) %>%
  arrange(abs_day, DAYSCONS, .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  select(-abs_day) %>%
  arrange(GEN_PUBLICID) %>%
  mutate(
    across(where(is.character), ~ na_if(str_squish(.), "N/A")),
    across(where(is.character), ~ na_if(., "Unknown")),
    AGE     = ifelse(AGE < 0, NA_integer_, AGE),
    MENDAYS = ifelse(MENDAYS < 0, NA_integer_, MENDAYS)
  )

STEP <- DEMO_base %>%
  transmute(
    ID = GEN_PUBLICID,
    AGE_AT_RANDOMIZATION = AGE,
    GENDER = GENDER,
    CURRENT_EMPLOYMENT_STATUS = EMPLOY,
    CURRENT_MARITAL_STATUS = MARSTAT,
    YEARS_OF_EDUCATION_COMPLETED = EDUCATE,
    ONSET_MENSES_DAYS_SINCE_RANDOMIZATION = MENDAYS,
    LIVEALONE = LIVEALON,
    ANNUAL_HOUSEHOLD_INCOME = HINCOME,
    RACE = RACE,
    ETHNICITY_HISPANIC_LATINO = ETHNICTY,
    ETHNICITY_SPEC = ETHNSPEC
  ) %>%
  mutate(
    mar_str = str_squish(as.character(CURRENT_MARITAL_STATUS)),
    Married = case_when(
      mar_str %in% c("Married", "Living as Married") ~ 1L,
      mar_str %in% c("", "Unknown") | is.na(mar_str) ~ NA_integer_,
      TRUE ~ 0L
    ),
    emp_str = str_squish(as.character(CURRENT_EMPLOYMENT_STATUS)),
    Employment = case_when(
      emp_str %in% c("Full-time", "Part-time for pay") ~ 1L,
      emp_str %in% c("Retired", "Unknown", "") | is.na(emp_str) ~ NA_integer_,
      TRUE ~ 0L
    ),
    edu_str = str_squish(as.character(YEARS_OF_EDUCATION_COMPLETED)),
    Education7 = case_when(
      edu_str == "Seventh grade - ninth grade"           ~ 1L,
      edu_str == "Partial High School"                   ~ 2L,
      edu_str == "High School Diploma or GED"            ~ 3L,
      edu_str == "Technical School or Associates Degree" ~ 4L,
      edu_str == "Some college (at least one year)"      ~ 5L,
      edu_str == "College Diploma (Bachelors Degree)"    ~ 6L,
      edu_str == "Graduate or Professional Degree"       ~ 7L,
      TRUE ~ NA_integer_
    ),
    Education7_3cat = case_when(
      Education7 %in% 1:3 ~ 1L,
      Education7 %in% 4:5 ~ 2L,
      Education7 %in% 6:7 ~ 3L,
      TRUE ~ NA_integer_
    ),
    Education7_3cat_ord = factor(
      Education7_3cat,
      levels = c(1, 2, 3),
      labels = c("Low", "Moderate", "High"),
      ordered = TRUE
    ),
    LIVEALONE = case_when(
      LIVEALONE == "Yes" ~ 1L,
      LIVEALONE == "No"  ~ 0L,
      TRUE               ~ NA_integer_
    ),
    INCOME = case_when(
      str_detect(ANNUAL_HOUSEHOLD_INCOME, "^A") ~ 1,
      str_detect(ANNUAL_HOUSEHOLD_INCOME, "^B") ~ 2,
      str_detect(ANNUAL_HOUSEHOLD_INCOME, "^C") ~ 3,
      str_detect(ANNUAL_HOUSEHOLD_INCOME, "^D") ~ 4,
      str_detect(ANNUAL_HOUSEHOLD_INCOME, "^E") ~ 5,
      str_detect(ANNUAL_HOUSEHOLD_INCOME, "^F") ~ 6,
      str_detect(ANNUAL_HOUSEHOLD_INCOME, "^G") ~ 7,
      str_detect(ANNUAL_HOUSEHOLD_INCOME, "^H") ~ 8,
      str_detect(ANNUAL_HOUSEHOLD_INCOME, "^I") ~ 9,
      str_detect(ANNUAL_HOUSEHOLD_INCOME, "^J") ~ 10,
      ANNUAL_HOUSEHOLD_INCOME == "Refused" ~ NA_real_,
      TRUE ~ NA_real_
    ),
    INCOME_3cat = case_when(
      INCOME %in% 1:3  ~ 1L,
      INCOME %in% 4:6  ~ 2L,
      INCOME %in% 7:10 ~ 3L,
      TRUE             ~ NA_integer_
    ),
    INCOME_3cat_ord = factor(
      INCOME_3cat,
      levels = c(1, 2, 3),
      labels = c("Low", "Moderate", "High"),
      ordered = TRUE
    ),
    race_str = str_squish(as.character(RACE)),
    Ethnicity = case_when(
      race_str == "White or Caucasian"        ~ 1L,
      race_str == "Black or African American" ~ 2L,
      TRUE ~ NA_integer_
    )
  ) %>%
  select(-mar_str, -emp_str, -edu_str, -race_str)

# -----------------------------
# ADE-based variables
# -----------------------------
ade_raw <- read_csv(file.path(clinical_dir, "ade_gen.csv"), show_col_types = FALSE)

ade_base <- ade_raw %>%
  mutate(
    ADE_ID = GEN_PUBLICID,
    ADE_ID_canon = canon_id(GEN_PUBLICID),
    DAYSCONS_abs = if_else(!is.na(DAYSCONS), abs(DAYSCONS), Inf)
  ) %>%
  arrange(ADE_ID_canon, DAYSCONS_abs, DAYSCONS) %>%
  distinct(ADE_ID_canon, .keep_all = TRUE) %>%
  select(-DAYSCONS_abs)

hx_vars <- intersect(
  c("HXPAN", "HXGAD", "HXSOC", "HXOCD", "HXPTSD", "HXAN", "HXBUL", "HXPD", "HXHYP"),
  names(ade_base)
)

map_def_or_yes <- function(x) {
  z <- as.character(x)
  as.numeric(z %in% c("Definite", "Yes"))
}

ade_comorb <- ade_base %>%
  transmute(
    ADE_ID_canon,
    across(all_of(hx_vars), map_def_or_yes)
  ) %>%
  mutate(
    Comorbidity = as.numeric(rowSums(across(all_of(hx_vars))) > 0)
  ) %>%
  select(ADE_ID_canon, Comorbidity)

ade_suicide <- ade_base %>%
  transmute(
    ADE_ID_canon,
    Suicide_Attempt = if ("HXSUI" %in% names(ade_base)) map_yesno(HXSUI) else NA_real_
  )

ade_clinical <- ade_base %>%
  transmute(
    ADE_ID_canon,
    Psychotic_lifetime = if ("PSYCHOS2" %in% names(ade_base)) map_yesno(PSYCHOS2) else NA_real_,
    Lifetime_N_phases = if ("LIFETPH" %in% names(ade_base)) parse_num(LIFETPH) else NA_real_,
    Age_first_depressive_episode = if ("DAGE1ST" %in% names(ade_base)) parse_num(DAGE1ST) else NA_real_,
    ADE_Anxiety_Severity = if ("ANXSEV" %in% names(ade_base)) {
      case_when(
        parse_num(ANXSEV) < 0 ~ NA_real_,
        parse_num(ANXSEV) > 4 ~ NA_real_,
        TRUE ~ parse_num(ANXSEV)
      )
    } else {
      NA_real_
    }
  )

STEP <- STEP %>%
  mutate(ID_canon = canon_id(ID)) %>%
  left_join(ade_comorb,   by = c("ID_canon" = "ADE_ID_canon")) %>%
  left_join(ade_suicide,  by = c("ID_canon" = "ADE_ID_canon")) %>%
  left_join(ade_clinical, by = c("ID_canon" = "ADE_ID_canon")) %>%
  select(-ID_canon) %>%
  mutate(
    ADE_Anx_3cat = case_when(
      ADE_Anxiety_Severity %in% 0:1 ~ 1L,
      ADE_Anxiety_Severity == 2     ~ 2L,
      ADE_Anxiety_Severity %in% 3:4 ~ 3L,
      TRUE                          ~ NA_integer_
    ),
    ADE_Anx_3cat_ord = factor(
      ADE_Anx_3cat,
      levels = c(1, 2, 3),
      ordered = TRUE
    ),
    Lifetime_N_phases_1to7 = case_when(
      Lifetime_N_phases %in% 0:4 ~ as.integer(Lifetime_N_phases),
      Lifetime_N_phases >= 5  & Lifetime_N_phases <= 11 ~ 5L,
      Lifetime_N_phases >= 12 & Lifetime_N_phases <= 53 ~ 6L,
      Lifetime_N_phases > 53 ~ 7L,
      TRUE ~ NA_integer_
    ),
    Lifetime_N_phases_3cat = case_when(
      Lifetime_N_phases %in% 0:5 ~ 1L,
      Lifetime_N_phases >= 6 & Lifetime_N_phases <= 12 ~ 2L,
      Lifetime_N_phases > 12 ~ 3L,
      TRUE ~ NA_integer_
    ),
    Lifetime_N_phases_3cat_ord = factor(
      Lifetime_N_phases_3cat,
      levels = c(1, 2, 3),
      ordered = TRUE
    )
  )

# -----------------------------
# MANPHASE and family history
# -----------------------------
recode_manphase_4cat <- function(x) {
  x <- trimws(as.character(x))
  case_when(
    x %in% c("0", "1", "2", "3 - 4", "5 - 9") ~ 1,
    x == "10 - 20" ~ 2,
    x == "20 - 50" ~ 3,
    x == "Too many to count" ~ 4,
    x %in% c("Indeterminate", "Unknown") ~ NA_real_,
    TRUE ~ NA_real_
  )
}

ade_manphase <- ade_raw %>%
  mutate(abs_day = abs(DAYSCONS)) %>%
  arrange(GEN_PUBLICID, abs_day) %>%
  group_by(GEN_PUBLICID) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    ID = GEN_PUBLICID,
    MANPHASE_4cat_ord = factor(recode_manphase_4cat(MANPHASE), levels = 1:4, ordered = TRUE)
  )

ade1_raw <- read_csv(file.path(clinical_dir, "ade1_gen.csv"), show_col_types = FALSE)

map_bpifamhx <- function(x) {
  z <- str_trim(as.character(x))
  case_when(
    str_detect(z, "^10\\s*=") | str_detect(z, "^15\\s*=") | str_detect(z, "^20\\s*=") ~ 1,
    str_detect(z, "^0\\s*=")  | str_detect(z, "^2\\s*=")  | str_detect(z, "^5\\s*=")  ~ 0,
    z %in% c("N/A", "Unknown", "") | is.na(z) ~ NA_real_,
    TRUE ~ NA_real_
  )
}

FH_BPI <- ade1_raw %>%
  mutate(ADE_ID_canon = canon_id(GEN_PUBLICID)) %>%
  arrange(ADE_ID_canon) %>%
  distinct(ADE_ID_canon, .keep_all = TRUE) %>%
  transmute(
    ADE_ID_canon,
    family_history_BD = map_bpifamhx(BPIFAMHX)
  )

STEP <- STEP %>%
  left_join(ade_manphase, by = "ID") %>%
  mutate(ID_canon = canon_id(ID)) %>%
  left_join(FH_BPI, by = c("ID_canon" = "ADE_ID_canon")) %>%
  select(-ID_canon)

# -----------------------------
# CMF-based variables
# -----------------------------
cmf_raw <- read_csv(file.path(clinical_dir, "cmf_gen.csv"), show_col_types = FALSE)

id_cmf <- pick_col(names(cmf_raw), exact = c("GEN_PUBLICID", "STEP_PATIENT_ID", "PUBLICID", "ID"))
rcp    <- pick_col(names(cmf_raw), exact = c("DAYSRZ", "DAYS_SINCE_RANDOMIZATION", "RCP"),
                   regex = "DAYS?_SINCE_?RANDOM")
scp    <- pick_col(names(cmf_raw), exact = c("DAYSCONS", "DAYS_SINCE_CONSENT", "SCP"),
                   regex = "DAYS?_SINCE_?CONS")
yr     <- pick_col(names(cmf_raw), exact = c("CMFYEAR"), regex = "^CMF.*YEAR$")

cmf_base <- cmf_raw %>%
  mutate(
    CMF_ID = .data[[id_cmf]],
    t_rcp  = if (!is.na(rcp)) parse_num(.data[[rcp]]) else NA_real_,
    t_scp  = if (!is.na(scp)) parse_num(.data[[scp]]) else NA_real_,
    t_year = if (!is.na(yr))  parse_num(.data[[yr]])  else NA_real_,
    time_key = coalesce(t_rcp, t_scp, t_year),
    CMF_ID_canon = canon_id(CMF_ID)
  ) %>%
  arrange(CMF_ID_canon, time_key) %>%
  group_by(CMF_ID_canon) %>%
  slice(1L) %>%
  ungroup()

nms <- names(cmf_base)

col_mendays <- pick_col(nms, exact = c("MENDAYS", "ONSET_MENSES_DAYS_SINCE_RZ"),
                        regex = "MENDAYS|ONSET.*MENSES")
col_gaf     <- pick_col(nms, exact = c("GAFMONTH", "GAF_THIS_MONTH"),
                        regex = "GAF.*(MONTH|MTH)$")
col_subab   <- pick_col(nms, exact = c("SUBABUSE", "CURR_SUBSTANCE_ABUSE_OR_DEPEND"),
                        regex = "SUB.*ABUSE.*DEPEND")
col_abnanx  <- pick_col(nms, exact = c("ABNANX", "ABNORMAL_ANXIETY"),
                        regex = "ABN.*ANXI|ANXI.*ABN")
col_depresd <- pick_col(nms, exact = c("DEPRESD", "DEPRESSED_MOST_OF_THE_DAY"),
                        regex = "DEPRESS.*MOST.*DAY")
col_hallu   <- pick_col(nms, exact = c("SMSHALLU"), regex = "HALLU")
col_delus   <- pick_col(nms, exact = c("SMSDELUS"), regex = "DELUS")

cmf_vars <- cmf_base %>%
  mutate(
    raw_abnanx  = if (!is.na(col_abnanx)) parse_num(.data[[col_abnanx]]) else NA_real_,
    raw_depresd = if (!is.na(col_depresd)) parse_num(.data[[col_depresd]]) else NA_real_
  ) %>%
  transmute(
    CMF_ID_canon,
    CGI_baseline = suppressWarnings(as.numeric(CGI)),
    ONSET_MENSES_DAYS_SINCE_RANDOMIZATION =
      if (!is.na(col_mendays)) parse_num(.data[[col_mendays]]) else NA_real_,
    GAF_THIS_MONTH =
      if (!is.na(col_gaf)) parse_num(.data[[col_gaf]]) else NA_real_,
    CURR_SUBSTANCE_ABUSE_OR_DEPEND =
      if (!is.na(col_subab)) map_yesno(.data[[col_subab]]) else NA_real_,
    CMF_Abnormal_Anxiety_Score = case_when(
      is.na(raw_abnanx) ~ NA_real_,
      raw_abnanx < 0    ~ NA_real_,
      raw_abnanx > 100  ~ NA_real_,
      TRUE              ~ raw_abnanx
    ),
    CMF_Depressed_Most_Day_Score = case_when(
      is.na(raw_depresd) ~ NA_real_,
      raw_depresd < 0    ~ NA_real_,
      raw_depresd > 100  ~ NA_real_,
      TRUE               ~ raw_depresd
    ),
    CMF_Hallucinations_Severity =
      if (!is.na(col_hallu)) parse_num(.data[[col_hallu]]) else NA_real_,
    CMF_Delusions_Severity =
      if (!is.na(col_delus)) parse_num(.data[[col_delus]]) else NA_real_
  ) %>%
  mutate(
    Psychotic_symptoms_score = pmax(CMF_Hallucinations_Severity, CMF_Delusions_Severity, na.rm = TRUE),
    Psychotic_symptoms_score = ifelse(is.infinite(Psychotic_symptoms_score), NA_real_, Psychotic_symptoms_score),
    Psychotic_symptoms_bin = case_when(
      is.na(CMF_Hallucinations_Severity) & is.na(CMF_Delusions_Severity) ~ NA_real_,
      CMF_Hallucinations_Severity > 0 | CMF_Delusions_Severity > 0       ~ 1,
      CMF_Hallucinations_Severity == 0 & CMF_Delusions_Severity == 0     ~ 0,
      TRUE ~ NA_real_
    )
  )

STEP <- STEP %>%
  mutate(ID_canon = canon_id(ID)) %>%
  select(-any_of(c(
    "CGI_baseline",
    "ONSET_MENSES_DAYS_SINCE_RANDOMIZATION",
    "GAF_THIS_MONTH",
    "CURR_SUBSTANCE_ABUSE_OR_DEPEND",
    "CMF_Abnormal_Anxiety_Score",
    "CMF_Depressed_Most_Day_Score",
    "CMF_Hallucinations_Severity",
    "CMF_Delusions_Severity",
    "Psychotic_symptoms_score",
    "Psychotic_symptoms_bin"
  ))) %>%
  left_join(cmf_vars, by = c("ID_canon" = "CMF_ID_canon")) %>%
  select(-ID_canon) %>%
  mutate(
    CGI_baseline_3cat = case_when(
      CGI_baseline %in% c(1, 2) ~ 1L,
      CGI_baseline %in% c(3, 4) ~ 2L,
      CGI_baseline %in% c(5, 6) ~ 3L,
      TRUE ~ NA_integer_
    ),
    CGI_baseline_3cat_ord = factor(
      CGI_baseline_3cat,
      levels = c(1, 2, 3),
      labels = c("Low", "Moderate", "High"),
      ordered = TRUE
    ),
    CMF_Anx_3cat = case_when(
      is.na(CMF_Abnormal_Anxiety_Score) ~ NA_integer_,
      CMF_Abnormal_Anxiety_Score == 0   ~ 1L,
      CMF_Abnormal_Anxiety_Score < 50   ~ 2L,
      TRUE                              ~ 3L
    ),
    CMF_Anx_3cat_ord = factor(
      CMF_Anx_3cat,
      levels = c(1, 2, 3),
      labels = c("Low", "Moderate", "High"),
      ordered = TRUE
    )
  )

# -----------------------------
# MADRS baseline total
# -----------------------------
madrs_raw <- read_csv(file.path(clinical_dir, "madrs_gen.csv"), show_col_types = FALSE)

has_daycons <- "DAYSCONS" %in% names(madrs_raw)
has_visit   <- "VISIT" %in% names(madrs_raw)

madrs_base <- madrs_raw %>%
  mutate(
    MADRS_ID  = GEN_PUBLICID,
    t_scp     = if (has_daycons) parse_num(.data[["DAYSCONS"]]) else NA_real_,
    visit_chr = if (has_visit) as.character(.data[["VISIT"]]) else NA_character_
  ) %>%
  group_by(MADRS_ID) %>%
  slice({
    idx_entry <- which(!is.na(visit_chr) & visit_chr == "Entry")
    if (length(idx_entry) >= 1) idx_entry[1] else which.min(ifelse(is.na(t_scp), Inf, t_scp))
  }) %>%
  ungroup()

madrs_item_cols <- c(
  "MADRSSAD", "MADRSPES", "MADRSSUI", "MADRSSLP", "MADRSFEE",
  "MADRSLAS", "MADRSCON", "MADRSTEN", "MADRSAPP", "MADRSAPS"
)

madrs_base <- madrs_base %>%
  mutate(MADRS_ID_canon = canon_id(MADRS_ID)) %>%
  mutate(
    across(
      all_of(intersect(madrs_item_cols, names(.))),
      ~ parse_num(.x),
      .names = "num_{.col}"
    )
  ) %>%
  rowwise() %>%
  mutate(
    MADRS_Total_Baseline =
      if (any(is.na(c_across(starts_with("num_"))))) NA_real_
    else sum(c_across(starts_with("num_")))
  ) %>%
  ungroup() %>%
  select(MADRS_ID_canon, MADRS_Total_Baseline)

STEP <- STEP %>%
  mutate(ID_canon = canon_id(ID)) %>%
  left_join(madrs_base, by = c("ID_canon" = "MADRS_ID_canon")) %>%
  select(-ID_canon)

# -----------------------------
# Social variables from LRIFT
# -----------------------------
lrift_raw <- read_csv(file.path(clinical_dir, "lrift_gen.csv"), show_col_types = FALSE)

id_lrift <- intersect(c("GEN_PUBLICID", "STEP_PATIENT_ID", "PUBLICID", "ID"), names(lrift_raw))[1]

lr_refriend <- lrift_raw %>%
  filter(VISIT == "Entry") %>%
  transmute(
    ID_canon = canon_id(.data[[id_lrift]]),
    REFRIEND = as.character(REFRIEND)
  ) %>%
  arrange(ID_canon) %>%
  group_by(ID_canon) %>%
  slice(1L) %>%
  ungroup()

STEP <- STEP %>%
  mutate(ID_canon = canon_id(ID)) %>%
  select(-any_of(c("REFRIEND", "Relationship_with_friends"))) %>%
  left_join(lr_refriend, by = "ID_canon") %>%
  mutate(
    Relationship_with_friends = case_when(
      parse_num(REFRIEND) %in% 1:5 ~ as.numeric(parse_num(REFRIEND)),
      TRUE                         ~ NA_real_
    ),
    Rel_friends_3cat = case_when(
      Relationship_with_friends %in% 1:2 ~ 1L,
      Relationship_with_friends == 3     ~ 2L,
      Relationship_with_friends %in% 4:5 ~ 3L,
      TRUE                               ~ NA_integer_
    ),
    Rel_friends_3cat_ord = factor(
      Rel_friends_3cat,
      levels = c(1, 2, 3),
      labels = c("Low", "Moderate", "High"),
      ordered = TRUE
    )
  ) %>%
  select(-ID_canon)

# -----------------------------
# PRSs
# -----------------------------
PRS_STEP <- read_table(
  file.path(prs_dir, "all_STEP_PRSs_standard"),
  col_names = TRUE,
  trim_ws = TRUE,
  col_types = cols(
    .default     = col_double(),
    IID          = col_character(),
    GEN_PUBLICID = col_character()
  )
)

PRS_keep <- PRS_STEP %>%
  mutate(
    id1 = sub("_.+$", "", IID),
    id2 = sub("^.+_", "", IID),
    .mismatch = !is.na(id1) & !is.na(id2) & id1 != id2
  )

if (any(PRS_keep$.mismatch, na.rm = TRUE)) {
  message("Some IID values do not have matching left and right identifiers. Using the left part for the merge.")
}

PRS_keep <- PRS_keep %>%
  mutate(key_canon = canon_id(id1)) %>%
  distinct(key_canon, .keep_all = TRUE) %>%
  select(
    key_canon, IID, GEN_PUBLICID,
    ADHD05_PRS, BODY04_PRS, BP_PRS, DEPR07_PRS, EDUC03_PRS, SCZ_PRS
  )

STEP <- STEP %>%
  mutate(key_canon_step = canon_id(ID)) %>%
  left_join(PRS_keep, by = join_by(key_canon_step == key_canon)) %>%
  select(-key_canon_step)

# -----------------------------
# PCs
# -----------------------------
pcs <- read_step_pcs(file.path(pheno_dir, "STEP_PCs"))

STEP <- STEP %>%
  mutate(ID_canon = canon_id(ID)) %>%
  left_join(
    pcs %>% select(ID_canon, starts_with("PC")),
    by = "ID_canon"
  ) %>%
  select(-ID_canon)

# -----------------------------
# PRS residualisation
# -----------------------------
prs_vars <- c("ADHD05_PRS", "BODY04_PRS", "BP_PRS", "DEPR07_PRS", "EDUC03_PRS", "SCZ_PRS")
pc_vars  <- c("PC1", "PC2", "PC3", "PC4")

for (prs in prs_vars) {
  STEP <- residualize_prs(STEP, prs_name = prs, pcs = pc_vars)
}

# -----------------------------
# TR1 from phenotype file
# -----------------------------
tr_fu <- read_ws(file.path(pheno_dir, "STEP_TR_follow_up.txt"))

if (!"TR1" %in% names(tr_fu)) {
  stop("TR1 column not found in STEP_TR_follow_up.txt. Rename the corresponding column to TR1 before running the script.")
}

STEP <- STEP %>%
  mutate(ID_chr = as.character(ID)) %>%
  left_join(
    tr_fu %>%
      mutate(GEN_PUBLICID = as.character(GEN_PUBLICID)) %>%
      select(GEN_PUBLICID, TR1),
    by = c("ID_chr" = "GEN_PUBLICID")
  ) %>%
  select(-ID_chr)

# -----------------------------
# TR2 from longitudinal CGI
# -----------------------------
cmf_long <- cmf_raw %>%
  mutate(
    ID_canon = canon_id(GEN_PUBLICID),
    CGI_num  = suppressWarnings(as.numeric(CGI)),
    day      = suppressWarnings(as.numeric(DAYSCONS))
  ) %>%
  filter(!is.na(day)) %>%
  arrange(ID_canon, day) %>%
  group_by(ID_canon) %>%
  mutate(
    lead_day = lead(day),
    interval = lead_day - day
  )

TR2_df <- cmf_long %>%
  summarise(
    n_visits_cgi = sum(!is.na(CGI_num)),
    n_timepoints = n_distinct(day[!is.na(CGI_num)]),
    min_day      = min(day, na.rm = TRUE),
    max_day      = max(day, na.rm = TRUE),
    total_time   = sum(interval[interval > 0], na.rm = TRUE),
    bad_time     = sum(interval[interval > 0 & !is.na(CGI_num) & CGI_num > 3], na.rm = TRUE),
    TR2 = case_when(
      n_visits_cgi < 2 ~ NA_integer_,
      n_timepoints < 2 ~ NA_integer_,
      max_day == min_day ~ NA_integer_,
      total_time <= 0 ~ NA_integer_,
      (bad_time / total_time) > 0.5 ~ 1L,
      TRUE ~ 0L
    ),
    .groups = "drop"
  )

STEP <- STEP %>%
  mutate(ID_canon = canon_id(ID)) %>%
  left_join(TR2_df %>% select(ID_canon, TR2), by = "ID_canon") %>%
  select(-ID_canon)

# -----------------------------
# Final save
# -----------------------------
write.csv(STEP, "STEP_preprocessed.csv", row.names = FALSE)