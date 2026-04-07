############################################################
# Project: Transdiagnostic predictors of treatment outcomes
# Cohort: CATIE
# Script: Data preprocessing
############################################################

library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(data.table)
library(fs)

options(stringsAsFactors = FALSE)

# -----------------------------
# Paths
# -----------------------------
data_dir  <- "INSERT_PATH_TO_CATIE_RAW_DATA"
prs_dir   <- "INSERT_PATH_TO_CATIE_PRS"
pheno_dir <- "INSERT_PATH_TO_CATIE_PHENOTYPES"

# -----------------------------
# Main dataset
# -----------------------------
CATIE <- read.csv("CATIE.csv", check.names = FALSE)

# -----------------------------
# Helpers
# -----------------------------
canon_id <- function(x) {
  gsub("[^A-Z0-9]", "", toupper(as.character(x)))
}

parse_num <- function(x) {
  suppressWarnings(readr::parse_number(as.character(x)))
}

pick_baseline <- function(df) {
  if ("BASE1" %in% names(df)) {
    df %>% filter(BASE1 == 1)
  } else if ("VISIT" %in% names(df)) {
    df %>% filter(grepl("^Baseline$", VISIT, ignore.case = TRUE))
  } else if (all(c("NIMHID", "VISITID") %in% names(df))) {
    df %>% group_by(NIMHID) %>% slice_min(VISITID, with_ties = FALSE) %>% ungroup()
  } else {
    df
  }
}

read_ws <- function(path) {
  x <- readr::read_table(path, col_names = TRUE, progress = FALSE)
  x <- x %>% mutate(across(where(is.character), ~ str_squish(.)))
  if ("type_convert" %in% getNamespaceExports("readr")) {
    x <- tryCatch(readr::type_convert(x), error = function(e) x)
  }
  x
}

normalize_pc_df <- function(x) {
  if (ncol(x) == 0) return(tibble())
  
  x <- x %>% mutate(across(where(is.character), ~ str_squish(.)))
  
  if (ncol(x) >= 2 && isTRUE(all(x[[1]] == x[[2]], na.rm = TRUE))) {
    x <- x[, -2, drop = FALSE]
  }
  
  names(x)[1] <- "GEN_PUBLICID"
  
  if (ncol(x) > 1) {
    names(x)[-1] <- paste0("PC", seq_len(ncol(x) - 1))
  }
  
  x <- x %>%
    mutate(across(starts_with("PC"), ~ suppressWarnings(as.numeric(.))))
  
  x$GEN_PUBLICID <- ifelse(
    grepl("_.+", x$GEN_PUBLICID),
    sub("_.+$", "", x$GEN_PUBLICID),
    x$GEN_PUBLICID
  )
  
  distinct(x, GEN_PUBLICID, .keep_all = TRUE)
}

read_cat_pcs_file <- function(path) {
  x <- tryCatch(readr::read_table(path, col_names = TRUE, progress = FALSE), error = function(e) NULL)
  
  bad_header <- !is.null(x) &&
    (any(grepl("^-?\\d", names(x))) || any(grepl("^[0-9]{2}-", names(x))))
  
  if (is.null(x) || bad_header) {
    x <- readr::read_table(path, col_names = FALSE, progress = FALSE)
  }
  
  normalize_pc_df(x)
}

read_CATIE_PCs <- function(path = file.path(pheno_dir, "CATIE_PCs")) {
  if (fs::file_exists(path) && !fs::is_dir(path)) {
    return(read_cat_pcs_file(path))
  }
  
  if (fs::is_dir(path)) {
    files <- dir_ls(path, recurse = TRUE, type = "file")
    if (!length(files)) return(tibble())
    pcs_list <- lapply(files, read_cat_pcs_file)
    return(bind_rows(pcs_list))
  }
  
  tibble()
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
# Baseline demographic variables
# -----------------------------
DEMO <- read.csv(file.path(data_dir, "DEMO.csv"), check.names = FALSE)

CATIE <- CATIE %>%
  mutate(
    EMPLOY   = DEMO$EMPLOY,
    MARITAL  = DEMO$MARITAL,
    EDUC_YRS = DEMO$EDUC_YRS,
    RACE     = DEMO$RACE,
    ID       = DEMO$NIMHID,
    SIGOTHER = DEMO$SIGOTHER
  ) %>%
  mutate(
    SIGOTHER = case_when(
      SIGOTHER == "Yes" ~ 1,
      SIGOTHER == "No"  ~ 0,
      TRUE              ~ NA_real_
    ),
    Married = ifelse(MARITAL == 1, 1L, 0L),
    EMPLOY_norm = str_trim(str_to_lower(EMPLOY)),
    Employment = case_when(
      EMPLOY_norm %in% c("full time", "part time") ~ 1L,
      EMPLOY_norm %in% c("did not work", "unknown", "") ~ 0L,
      is.na(EMPLOY_norm) ~ NA_integer_,
      TRUE ~ 0L
    ),
    Ethnicity = ifelse(RACE %in% c(1, 2), RACE, NA)
  ) %>%
  select(-EMPLOY_norm)

# -----------------------------
# Income
# -----------------------------
SURFQ <- read.csv(file.path(data_dir, "SURFQ_2.csv"), check.names = FALSE)

SURFQ_B <- SURFQ %>%
  filter(VISIT == "Baseline") %>%
  select(NIMHID, INCOME) %>%
  distinct(NIMHID, .keep_all = TRUE)

CATIE <- CATIE %>%
  left_join(SURFQ_B, by = c("ID" = "NIMHID")) %>%
  mutate(
    INCOME_bin = case_when(
      INCOME %in% 1:2 ~ 0L,
      INCOME %in% 3:5 ~ 1L,
      TRUE            ~ NA_integer_
    ),
    INCOME_cat = factor(
      case_when(
        INCOME %in% c(1, 2)     ~ "Low",
        INCOME %in% c(3, 4, 5)  ~ "High",
        TRUE                    ~ NA_character_
      ),
      levels = c("Low", "High"),
      ordered = TRUE
    )
  )

# -----------------------------
# BMI
# -----------------------------
VITAL <- read.csv(file.path(data_dir, "VITAL.csv"), check.names = FALSE)

CATIE <- CATIE %>%
  left_join(
    VITAL %>% select(NIMHID, BMI_CAT) %>% distinct(NIMHID, .keep_all = TRUE),
    by = c("ID" = "NIMHID")
  ) %>%
  mutate(
    BMI_CAT = case_when(
      BMI_CAT %in% c(97, 98, 99) ~ NA_real_,
      TRUE                       ~ as.numeric(BMI_CAT)
    ),
    BMI_3cat_num = case_when(
      BMI_3cat == "normo_sotto" ~ 1L,
      BMI_3cat == "sovrappeso"  ~ 2L,
      BMI_3cat == "obeso"       ~ 3L,
      TRUE                      ~ NA_integer_
    ),
    BMI_3cat_ord = factor(
      BMI_3cat_num,
      levels = c(1, 2, 3),
      labels = c("Normal", "Overweight", "Obese"),
      ordered = TRUE
    )
  )

# -----------------------------
# Social withdrawal
# -----------------------------
FBASELN1 <- read.csv(file.path(data_dir, "FBASELN1.csv"), check.names = FALSE)

CATIE <- CATIE %>%
  left_join(
    FBASELN1 %>% select(NIMHID, FI17) %>% distinct(NIMHID, .keep_all = TRUE),
    by = c("ID" = "NIMHID")
  ) %>%
  rename(Social_Withdrawal = FI17) %>%
  mutate(
    Social_Withdrawal_num = case_when(
      Social_Withdrawal %in% 0:3 ~ as.numeric(Social_Withdrawal),
      Social_Withdrawal %in% c(97, 98, 99) ~ NA_real_,
      is.na(Social_Withdrawal) ~ NA_real_,
      TRUE ~ NA_real_
    ),
    Social_Withdrawal_ord = factor(
      Social_Withdrawal_num,
      levels = c(0, 1, 2, 3),
      ordered = TRUE
    )
  )

# Optional PANSS-based social withdrawal
PANSS <- read.csv(file.path(data_dir, "PANSS.csv"), check.names = FALSE)

panss_sw <- PANSS %>%
  filter(VISIT == "Baseline") %>%
  transmute(
    NIMHID,
    Social_Withdrawal_PANSS = parse_num(PANSSN04)
  )

CATIE <- CATIE %>%
  left_join(panss_sw, by = c("ID" = "NIMHID")) %>%
  mutate(
    Social_Withdrawal = coalesce(Social_Withdrawal_PANSS, Social_Withdrawal_num),
    Social_Withdrawal_3cat = case_when(
      is.na(Social_Withdrawal)   ~ NA_integer_,
      Social_Withdrawal %in% 1:2 ~ 1L,
      Social_Withdrawal == 3     ~ 2L,
      Social_Withdrawal %in% 4:7 ~ 3L,
      TRUE                       ~ NA_integer_
    ),
    Social_Withdrawal_3cat_ord = factor(
      Social_Withdrawal_3cat,
      levels = c(1, 2, 3),
      labels = c("Low", "Moderate", "High"),
      ordered = TRUE
    )
  ) %>%
  select(-Social_Withdrawal_PANSS)

# -----------------------------
# Family history and substance abuse
# -----------------------------
SCID_2 <- read.csv(file.path(data_dir, "SCID_2.csv"), check.names = FALSE)
SCID_1 <- read.csv(file.path(data_dir, "SCID_1.csv"), check.names = FALSE)

CATIE <- CATIE %>%
  left_join(
    SCID_2 %>% select(NIMHID, SCID22) %>% distinct(NIMHID, .keep_all = TRUE),
    by = c("ID" = "NIMHID")
  ) %>%
  rename(Family_History = SCID22) %>%
  mutate(
    Family_History = trimws(as.character(Family_History)),
    Family_History = ifelse(
      Family_History %in% c("Yes", "No"),
      as.integer(Family_History == "Yes"),
      NA_integer_
    )
  ) %>%
  left_join(
    SCID_1 %>% select(NIMHID, SCID11) %>% distinct(NIMHID, .keep_all = TRUE),
    by = c("ID" = "NIMHID")
  ) %>%
  rename(Substance_Abuse = SCID11) %>%
  mutate(
    Substance_Abuse = tolower(trimws(Substance_Abuse)),
    Substance_Abuse = ifelse(
      Substance_Abuse == "yes", 1,
      ifelse(Substance_Abuse == "no", 0, NA_real_)
    )
  )

# -----------------------------
# Comorbidity
# -----------------------------
SCREEN <- read.csv(file.path(data_dir, "SCREEN.csv"), check.names = FALSE)

CATIE <- CATIE %>%
  left_join(
    SCREEN %>% select(NIMHID, SCRN09K) %>% distinct(NIMHID, .keep_all = TRUE),
    by = c("ID" = "NIMHID")
  ) %>%
  rename(Comorbidity = SCRN09K)

# -----------------------------
# Previous episodes
# -----------------------------
DEMO_SCID19 <- DEMO %>%
  select(NIMHID, SCID19) %>%
  distinct(NIMHID, .keep_all = TRUE)

CATIE <- CATIE %>%
  left_join(DEMO_SCID19, by = c("ID" = "NIMHID")) %>%
  rename(N_previous_episodes = SCID19) %>%
  mutate(
    N_previous_episodes = trimws(as.character(N_previous_episodes)),
    N_previous_episodes = na_if(N_previous_episodes, ""),
    N_previous_episodes = suppressWarnings(as.numeric(N_previous_episodes)),
    N_previous_episodes_num = case_when(
      N_previous_episodes %in% c(97, 98, 99) ~ NA_real_,
      is.na(N_previous_episodes)             ~ NA_real_,
      TRUE                                   ~ as.numeric(N_previous_episodes)
    ),
    N_prev_episodes_3cat = case_when(
      N_previous_episodes %in% c(0, 1) ~ 1L,
      N_previous_episodes %in% c(2, 3) ~ 2L,
      N_previous_episodes == 4         ~ 3L,
      TRUE                             ~ NA_integer_
    ),
    N_prev_episodes_3cat_ord = factor(
      N_prev_episodes_3cat,
      levels = c(1, 2, 3),
      labels = c("Low", "Moderate", "High"),
      ordered = TRUE
    )
  )

# -----------------------------
# Previously screened
# -----------------------------
SCRN00_only <- SCREEN %>%
  select(NIMHID, SCRN00) %>%
  mutate(
    raw = as.character(SCRN00),
    is_label = grepl("^Dataset\\b", raw, ignore.case = TRUE),
    val = gsub("\\s+", " ", trimws(raw)),
    Previously_Screened = case_when(
      is_label ~ NA_real_,
      val == "" | is.na(raw) ~ 1,
      val %in% c("1", "Yes", "YES", "yes", "Y") ~ 1,
      val %in% c("0", "No", "NO", "no", "N")    ~ 0,
      grepl("^[01]$", val) ~ as.numeric(val),
      TRUE ~ NA_real_
    )
  ) %>%
  select(NIMHID, Previously_Screened) %>%
  distinct(NIMHID, .keep_all = TRUE)

CATIE <- CATIE %>%
  left_join(SCRN00_only, by = c("ID" = "NIMHID"))

# -----------------------------
# Hospitalization
# -----------------------------
SURFM_1 <- read.csv(file.path(data_dir, "SURFM_1.csv"), check.names = FALSE)

SURFM1_bl <- SURFM_1 %>%
  pick_baseline() %>%
  select(NIMHID, INPAT1) %>%
  mutate(
    raw = as.character(INPAT1),
    is_label = grepl("^Dataset\\b", raw, ignore.case = TRUE),
    val = gsub("\\s+", " ", trimws(raw)),
    Hospitalization = case_when(
      is_label ~ NA_real_,
      val %in% c("1", "Yes", "YES", "yes", "Y") ~ 1,
      val %in% c("0", "No", "NO", "no", "N")    ~ 0,
      grepl("^[01]$", val) ~ as.numeric(val),
      TRUE ~ NA_real_
    )
  ) %>%
  select(NIMHID, Hospitalization) %>%
  distinct(NIMHID, .keep_all = TRUE)

CATIE <- CATIE %>%
  left_join(SURFM1_bl, by = c("ID" = "NIMHID"))

# -----------------------------
# Symptom severity
# -----------------------------
CGIS <- read.csv(file.path(data_dir, "CGIS.csv"), check.names = FALSE)

CGIS_bl <- CGIS %>%
  filter(BASE1 == 1 | grepl("^Baseline$", VISIT, ignore.case = TRUE)) %>%
  transmute(
    NIMHID = as.character(NIMHID),
    CGI_Severity = suppressWarnings(
      as.numeric(sub("^\\s*([1-7]).*$", "\\1", as.character(CS16)))
    )
  ) %>%
  distinct(NIMHID, .keep_all = TRUE) %>%
  mutate(CGI_Severity = ifelse(CGI_Severity %in% c(97, 98, 99), NA, CGI_Severity))

CATIE <- CATIE %>%
  mutate(ID_chr = as.character(ID)) %>%
  left_join(CGIS_bl, by = c("ID_chr" = "NIMHID")) %>%
  select(-ID_chr) %>%
  mutate(
    CGI_Severity_3cat = case_when(
      CGI_Severity %in% c(1, 2)    ~ 1L,
      CGI_Severity %in% c(3, 4)    ~ 2L,
      CGI_Severity %in% c(5, 6, 7) ~ 3L,
      TRUE                         ~ NA_integer_
    ),
    CGI_Severity_3cat_ord = factor(
      CGI_Severity_3cat,
      levels = c(1, 2, 3),
      labels = c("Low", "Moderate", "High"),
      ordered = TRUE
    )
  )

# -----------------------------
# Functioning
# -----------------------------
QOL <- read.csv(file.path(data_dir, "QOL.csv"), check.names = FALSE)

QOL_bl <- QOL %>%
  filter(BASE1 == 1 | (!"BASE1" %in% names(.) & grepl("^Baseline$", VISIT, ignore.case = TRUE))) %>%
  select(NIMHID, QOL_TOT) %>%
  distinct(NIMHID, .keep_all = TRUE)

CATIE <- CATIE %>%
  left_join(QOL_bl, by = c("ID" = "NIMHID"))

if (any(grepl("^QOL_TOT\\.(x|y)$", names(CATIE)))) {
  CATIE <- CATIE %>%
    mutate(QOL_TOT = coalesce(QOL_TOT, QOL_TOT.x, QOL_TOT.y)) %>%
    select(-ends_with(".x"), -ends_with(".y"))
}

CATIE <- CATIE %>%
  rename(Functioning = QOL_TOT)

# -----------------------------
# Anxiety
# -----------------------------
PANSS_bl <- PANSS %>%
  filter(BASE1 == 1 | (!"BASE1" %in% names(.) & grepl("^Baseline$", VISIT, ignore.case = TRUE))) %>%
  select(NIMHID, PANSSG02) %>%
  distinct(NIMHID, .keep_all = TRUE)

CATIE <- CATIE %>%
  left_join(PANSS_bl, by = c("ID" = "NIMHID"))

if (any(grepl("^PANSSG02\\.(x|y)$", names(CATIE)))) {
  CATIE <- CATIE %>%
    mutate(PANSSG02 = coalesce(PANSSG02, PANSSG02.x, PANSSG02.y)) %>%
    select(-ends_with(".x"), -ends_with(".y"))
}

CATIE <- CATIE %>%
  rename(Anxiety = PANSSG02) %>%
  mutate(
    Anxiety = case_when(
      Anxiety %in% c(97, 98, 99) ~ NA_real_,
      TRUE                       ~ as.numeric(Anxiety)
    ),
    Anxiety_3cat_num = case_when(
      Anxiety %in% c(1, 2)    ~ 1L,
      Anxiety %in% c(3, 4)    ~ 2L,
      Anxiety %in% c(5, 6, 7) ~ 3L,
      TRUE                    ~ NA_integer_
    ),
    Anxiety_3cat_ord = factor(
      Anxiety_3cat_num,
      levels = c(1, 2, 3),
      labels = c("Low", "Moderate", "High"),
      ordered = TRUE
    )
  )

# -----------------------------
# Depressive symptoms
# -----------------------------
CLGRY <- read.csv(file.path(data_dir, "CLGRY.csv"), check.names = FALSE)

CLGRY_bl <- CLGRY %>%
  filter(BASE1 == 1 | (!"BASE1" %in% names(.) & VISIT == "Baseline")) %>%
  select(NIMHID, CALG1) %>%
  distinct(NIMHID, .keep_all = TRUE)

map_calgry <- c("Absent" = 0, "Mild" = 1, "Moderate" = 2, "Severe" = 3)

CLGRY_bl <- CLGRY_bl %>%
  mutate(
    CALG1_num = unname(map_calgry[CALG1]),
    CALG1_num = as.numeric(CALG1_num)
  )

CATIE <- CATIE %>%
  left_join(CLGRY_bl %>% select(NIMHID, CALG1_num), by = c("ID" = "NIMHID")) %>%
  rename(Calgary_Depressed_Mood = CALG1_num) %>%
  mutate(
    Calgary_Depressed_Mood_bin = ifelse(
      Calgary_Depressed_Mood == 0,
      0,
      ifelse(Calgary_Depressed_Mood %in% c(1, 2, 3), 1, NA)
    )
  )

# -----------------------------
# Psychotic symptoms
# -----------------------------
PANSS_b1 <- PANSS %>%
  filter(BASE1 == 1 | (!"BASE1" %in% names(.) & VISIT == "Baseline")) %>%
  select(NIMHID, B1_PPOS) %>%
  distinct(NIMHID, .keep_all = TRUE)

CATIE <- CATIE %>%
  left_join(PANSS_b1, by = c("ID" = "NIMHID")) %>%
  rename(PANSS_Psychotic = B1_PPOS)

# -----------------------------
# Suicidality and suicide attempt
# -----------------------------
CLGRY_su_bl <- CLGRY %>%
  filter(BASE1 == 1 | (!"BASE1" %in% names(.) & VISIT == "Baseline")) %>%
  select(NIMHID, CALG8) %>%
  distinct(NIMHID, .keep_all = TRUE) %>%
  mutate(
    v = trimws(as.character(CALG8)),
    n_from_str = suppressWarnings(as.numeric(sub("^\\s*([1-4]).*$", "\\1", v))),
    n_from_lab = unname(c("Absent" = 1, "Mild" = 2, "Moderate" = 3, "Severe" = 4)[v]),
    Suicidality = ifelse(!is.na(n_from_str) & n_from_str %in% 1:4, n_from_str, n_from_lab),
    Suicidality = as.numeric(Suicidality)
  ) %>%
  select(NIMHID, Suicidality)

CATIE <- CATIE %>%
  left_join(CLGRY_su_bl, by = c("ID" = "NIMHID")) %>%
  mutate(
    Suicidality_bin = case_when(
      Suicidality == 1     ~ 0L,
      Suicidality %in% 2:4 ~ 1L,
      TRUE                 ~ NA_integer_
    ),
    Suicidality_cat = factor(
      Suicidality_bin,
      levels = c(0, 1),
      labels = c("No", "Yes")
    )
  )

FBASELN3 <- read.csv(file.path(data_dir, "FBASELN3.csv"), check.names = FALSE)

FBASELN3_bl <- {
  df <- FBASELN3
  if ("BASE1" %in% names(df)) {
    df <- filter(df, BASE1 == 1)
  } else if ("VISIT" %in% names(df)) {
    df <- filter(df, grepl("^Baseline$", VISIT, ignore.case = TRUE))
  }
  
  df %>%
    select(NIMHID, FI78) %>%
    distinct(NIMHID, .keep_all = TRUE) %>%
    mutate(
      raw = trimws(as.character(FI78)),
      code = suppressWarnings(as.numeric(sub("^\\D*([0-9]+).*$", "\\1", raw))),
      code = case_when(
        !is.na(code)                 ~ code,
        tolower(raw) %in% c("yes", "y") ~ 1,
        tolower(raw) %in% c("no", "n")  ~ 0,
        TRUE ~ NA_real_
      ),
      Suicide_Attempt = ifelse(code %in% c(0, 1, 97, 98, 99), code, NA_real_)
    ) %>%
    select(NIMHID, Suicide_Attempt)
}

CATIE <- CATIE %>%
  left_join(FBASELN3_bl, by = c("ID" = "NIMHID"))

# -----------------------------
# Adverse events
# -----------------------------
AESPO <- read.csv(file.path(data_dir, "AESPO.csv"), check.names = FALSE)

AESPO_bl <- {
  df <- AESPO
  if ("BL1_VIS" %in% names(df)) {
    df <- df %>% filter(BL1_VIS == 1)
  } else {
    df <- df %>% filter(VISIT == "Baseline")
  }
  
  df %>%
    select(NIMHID, ANYGENAE, PTSEV) %>%
    distinct(NIMHID, .keep_all = TRUE) %>%
    mutate(
      any_raw = trimws(as.character(ANYGENAE)),
      any_num = suppressWarnings(as.numeric(sub("^\\s*([0-9]+).*$", "\\1", any_raw))),
      any_num = case_when(
        !is.na(any_num)  ~ any_num,
        any_raw == "Yes" ~ 1,
        any_raw == "No"  ~ 0,
        TRUE             ~ NA_real_
      ),
      Any_AE_Baseline = ifelse(any_num %in% c(0, 1), any_num, NA_real_),
      
      sev_raw = trimws(as.character(PTSEV)),
      sev_num = suppressWarnings(as.numeric(sub("^\\s*([0-9]+).*$", "\\1", sev_raw))),
      sev_num = case_when(
        !is.na(sev_num)             ~ sev_num,
        PTSEV == "Mild"             ~ 1,
        PTSEV == "Moderate"         ~ 2,
        PTSEV == "Severe"           ~ 3,
        PTSEV == "Life-threatening" ~ 4,
        TRUE                        ~ NA_real_
      ),
      AE_Severity_Baseline = ifelse(sev_num %in% 1:4, sev_num, NA_real_)
    ) %>%
    select(NIMHID, Any_AE_Baseline, AE_Severity_Baseline)
}

CATIE <- CATIE %>%
  left_join(AESPO_bl, by = c("ID" = "NIMHID"))

AESYS <- read.csv(file.path(data_dir, "AESYS.csv"), check.names = FALSE)

AESYS_bl <- AESYS %>%
  filter(VISIT == "Baseline") %>%
  select(NIMHID, PTAE) %>%
  distinct(NIMHID, .keep_all = TRUE) %>%
  mutate(
    raw = trimws(as.character(PTAE)),
    num = suppressWarnings(as.numeric(sub("^\\s*([0-9]+).*$", "\\1", raw))),
    num = case_when(
      !is.na(num)  ~ num,
      raw == "Yes" ~ 1,
      raw == "No"  ~ 0,
      TRUE         ~ NA_real_
    ),
    PTAE_baseline = ifelse(num %in% c(0, 1), num, NA_real_)
  ) %>%
  select(NIMHID, PTAE_baseline)

CATIE <- CATIE %>%
  left_join(AESYS_bl, by = c("ID" = "NIMHID"))

# -----------------------------
# PRSs
# -----------------------------
PRS_CATIE <- read_table(
  file.path(prs_dir, "all_CATIE_PRSs_standard"),
  col_types = cols(
    IID     = col_character(),
    SampID  = col_character(),
    NIMHID  = col_character(),
    .default = col_double()
  )
)

prs_vars <- c(
  "ADHD05_PRS", "BODY04_PRS", "BP_PRS",
  "DEPR07_PRS", "EDUC03_PRS", "SCZ_PRS",
  "IID", "SampID"
)

PRS_keep <- PRS_CATIE %>%
  mutate(NIMHID = as.character(NIMHID)) %>%
  select(NIMHID, all_of(prs_vars))

CATIE <- CATIE %>%
  mutate(NIMHID = as.character(ID)) %>%
  left_join(PRS_keep, by = "NIMHID") %>%
  select(-NIMHID)

# -----------------------------
# Principal components
# -----------------------------
pcs <- read_CATIE_PCs(file.path(pheno_dir, "CATIE_PCs"))

pcs_join <- pcs %>%
  mutate(PC_KEY = canon_id(GEN_PUBLICID)) %>%
  select(PC_KEY, starts_with("PC")) %>%
  distinct(PC_KEY, .keep_all = TRUE)

CATIE_keys <- CATIE %>%
  mutate(
    K_SAMPID = if ("SampID" %in% names(.)) canon_id(SampID) else NA_character_,
    K_IID    = if ("IID" %in% names(.)) canon_id(sub("_[^_]+$", "", IID)) else NA_character_
  )

n_match_sampid <- inner_join(
  CATIE_keys %>% filter(!is.na(K_SAMPID)) %>% select(K_SAMPID),
  pcs_join %>% select(PC_KEY),
  by = c("K_SAMPID" = "PC_KEY")
) %>% nrow()

n_match_iid <- inner_join(
  CATIE_keys %>% filter(!is.na(K_IID)) %>% select(K_IID),
  pcs_join %>% select(PC_KEY),
  by = c("K_IID" = "PC_KEY")
) %>% nrow()

best_key <- if (n_match_sampid >= n_match_iid) "K_SAMPID" else "K_IID"

pc_cols_exist <- grep("^PC\\d+$", names(CATIE_keys), value = TRUE)

CATIE <- CATIE_keys %>%
  select(-any_of(pc_cols_exist)) %>%
  left_join(pcs_join, by = setNames("PC_KEY", best_key)) %>%
  select(-K_SAMPID, -K_IID)

# -----------------------------
# PRS residualisation
# -----------------------------
prs_resid_vars <- c("ADHD05_PRS", "BODY04_PRS", "BP_PRS", "DEPR07_PRS", "EDUC03_PRS", "SCZ_PRS")
pc_vars <- c("PC1", "PC2", "PC3", "PC4")

for (prs in prs_resid_vars) {
  CATIE <- residualize_prs(CATIE, prs_name = prs, pcs = pc_vars)
}

# -----------------------------
# TR1 from phenotype file
# -----------------------------
TRS <- read_ws(file.path(pheno_dir, "CATIE_TRS.txt"))

# Replace 'TR1_original' with the actual resistance variable in CATIE_TRS.txt if needed
if (!"TR1" %in% names(TRS)) {
  stop("TR1 column not found in CATIE_TRS.txt. Rename the corresponding column to TR1 before running the script.")
}

CATIE <- CATIE %>%
  mutate(ID_chr = as.character(ID)) %>%
  left_join(
    TRS %>%
      mutate(NIMHID = as.character(NIMHID)) %>%
      select(NIMHID, TR1),
    by = c("ID_chr" = "NIMHID")
  ) %>%
  select(-ID_chr)

# -----------------------------
# TR2 from longitudinal CGI
# -----------------------------
CGIS_long <- CGIS %>%
  mutate(
    NIMHID  = as.character(NIMHID),
    CGI_num = suppressWarnings(as.numeric(sub("^\\s*([1-7]).*$", "\\1", as.character(CS16))))
  )

TR2_df <- CGIS_long %>%
  group_by(NIMHID) %>%
  summarise(
    n_vis = sum(!is.na(CGI_num)),
    all_gt3 = ifelse(n_vis >= 2, all(CGI_num > 3, na.rm = TRUE), NA),
    TR2 = case_when(
      is.na(all_gt3) ~ NA_integer_,
      all_gt3        ~ 1L,
      TRUE           ~ 0L
    ),
    .groups = "drop"
  )

CATIE <- CATIE %>%
  mutate(ID_chr = as.character(ID)) %>%
  left_join(TR2_df, by = c("ID_chr" = "NIMHID")) %>%
  select(-ID_chr)

# -----------------------------
# Final recoding
# -----------------------------
CATIE <- CATIE %>%
  mutate(
    MARITAL = case_when(
      MARITAL == "Married"       ~ 1,
      MARITAL == "Widowed"       ~ 2,
      MARITAL == "Divorced"      ~ 3,
      MARITAL == "Separated"     ~ 4,
      MARITAL == "Never married" ~ 5,
      TRUE ~ NA_real_
    ),
    RACE = case_when(
      RACE == "White alone" ~ 1,
      RACE == "Black or African-American alone" ~ 2,
      RACE == "American Indian or Alasks Native alone" ~ 3,
      RACE == "Asian alone" ~ 4,
      RACE == "Native Hawaiian or Other Pasific Islander alone" ~ 5,
      RACE == "Two or More Races" ~ 6,
      TRUE ~ NA_real_
    )
  )

# -----------------------------
# Save output
# -----------------------------
write.csv(CATIE, "CATIE_preprocessed.csv", row.names = FALSE)