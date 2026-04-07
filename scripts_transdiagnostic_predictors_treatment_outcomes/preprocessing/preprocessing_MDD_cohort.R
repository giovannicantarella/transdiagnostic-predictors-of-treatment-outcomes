############################################################
# Project: Transdiagnostic predictors of treatment outcomes
# Cohort: GSRD
# Script: Data preprocessing
############################################################

library(dplyr)
library(readr)
library(stringr)
library(tidyr)

options(stringsAsFactors = FALSE)

# -----------------------------
# Input files
# -----------------------------
# Replace these with the local paths to the GSRD source files

gsrd_file     <- "path/to/GSRD.csv"
features_file <- "path/to/TRD3_MDD_features.csv"

# -----------------------------
# Main datasets
# -----------------------------
GSRD <- read.csv(gsrd_file, stringsAsFactors = FALSE)
GSRD_MDD_feat <- read.csv(features_file, sep = ";", stringsAsFactors = FALSE)

# -----------------------------
# Helpers
# -----------------------------
make_resid_on_covars <- function(df, y, covars, min_cc = 10) {
  covars <- covars[covars %in% names(df)]
  dat <- df[, c(y, covars), drop = FALSE]
  cc <- stats::complete.cases(dat)
  
  out <- rep(NA_real_, nrow(df))
  if (sum(cc) < min_cc) return(out)
  
  fml <- stats::as.formula(paste0(y, " ~ ", paste(covars, collapse = " + ")))
  fit <- stats::lm(fml, data = df, subset = cc)
  out[cc] <- stats::residuals(fit)
  out
}

# -----------------------------
# Merge additional feature file
# -----------------------------
GSRD_MDD_feat <- GSRD_MDD_feat %>%
  mutate(GenID = paste0("TRD3", GenID))

GSRD <- GSRD %>%
  left_join(
    GSRD_MDD_feat %>%
      select(GenID, setdiff(names(GSRD_MDD_feat), names(GSRD))),
    by = "GenID"
  )

# -----------------------------
# Psychotic features
# -----------------------------
GSRD <- GSRD %>%
  rename(MINI_MDD_specifiers_raw = DiagNostic.Interview..M.I.N.I.._P3) %>%
  mutate(
    Psychotic_features_MINI = case_when(
      is.na(MINI_MDD_specifiers_raw) ~ NA_character_,
      str_detect(str_to_lower(MINI_MDD_specifiers_raw), "with psychotic features") ~ "Yes",
      TRUE ~ "No"
    ),
    Psychotic_features_MINI = factor(Psychotic_features_MINI, levels = c("No", "Yes")),
    Psychotic_MADRS = case_when(
      is.na(Pessimistic.thoughts.MADRS.retrospective) ~ NA_character_,
      Pessimistic.thoughts.MADRS.retrospective %in% c(5, 6) ~ "Yes",
      TRUE ~ "No"
    ),
    Psychotic_MADRS = factor(Psychotic_MADRS, levels = c("No", "Yes")),
    Psychotic_YMRS11 = case_when(
      is.na(Insight.Young.Current) ~ NA_character_,
      Insight.Young.Current == 0   ~ "No",
      Insight.Young.Current > 0    ~ "Yes"
    ),
    Psychotic_YMRS11 = factor(Psychotic_YMRS11, levels = c("No", "Yes"))
  )

# -----------------------------
# Basic cleaning
# -----------------------------
ids_bad_bmi <- c(254, 967, 1422, 140, 1356)

GSRD <- GSRD %>%
  mutate(
    bmi = ifelse(ID %in% ids_bad_bmi, NA, bmi),
    Age.at.Onset_clean = ifelse(Age.at.Onset <= 10, NA, Age.at.Onset)
  )

# -----------------------------
# Psychiatric comorbidity
# -----------------------------
GSRD <- GSRD %>%
  mutate(
    Psychiatric_Comorbidity = case_when(
      PanicDisorder == 1 | GAD == 1 | OCD == 1 | PTSD == 1 ~ 1,
      is.na(PanicDisorder) & is.na(GAD) & is.na(OCD) & is.na(PTSD) ~ NA_real_,
      TRUE ~ 0
    ),
    Psychiatric_Comorbidity = factor(
      Psychiatric_Comorbidity,
      levels = c(0, 1),
      labels = c("No", "Yes")
    )
  )

# -----------------------------
# Employment and functioning
# -----------------------------
GSRD <- GSRD %>%
  mutate(
    Employed = case_when(
      Work_Status %in% c("Full time", "Part time", "Other",
                         "Sheltered workshop", "Volunteer work", "Student") ~ 1,
      Work_Status %in% c("Unemployed") ~ 0,
      TRUE ~ NA_real_
    ),
    Employed = factor(
      Employed,
      levels = c(0, 1),
      labels = c("No", "Yes")
    ),
    Functioning = rowSums(
      select(., Sheehan.Work, Sheehan.Social, Sheehan.Family),
      na.rm = TRUE
    )
  )

# -----------------------------
# Education
# -----------------------------
GSRD <- GSRD %>%
  mutate(
    Education_3cat = case_when(
      Education %in% c(1, 2) ~ 1L,
      Education == 3         ~ 2L,
      Education == 4         ~ 3L,
      TRUE ~ NA_integer_
    ),
    Education_3cat_ord = factor(
      Education_3cat,
      levels = c(1, 2, 3),
      labels = c("Low", "Moderate", "High"),
      ordered = TRUE
    )
  )

# -----------------------------
# Sheehan impairment
# -----------------------------
GSRD <- GSRD %>%
  mutate(
    SheehanWork_3cat = case_when(
      Sheehan.Work %in% 0:4  ~ 1L,
      Sheehan.Work %in% 5:7  ~ 2L,
      Sheehan.Work %in% 8:10 ~ 3L,
      TRUE ~ NA_integer_
    ),
    SheehanSocial_3cat = case_when(
      Sheehan.Social %in% 0:4  ~ 1L,
      Sheehan.Social %in% 5:7  ~ 2L,
      Sheehan.Social %in% 8:10 ~ 3L,
      TRUE ~ NA_integer_
    ),
    SheehanFamily_3cat = case_when(
      Sheehan.Family %in% 0:4  ~ 1L,
      Sheehan.Family %in% 5:7  ~ 2L,
      Sheehan.Family %in% 8:10 ~ 3L,
      TRUE ~ NA_integer_
    ),
    SheehanWork_3cat_ord = factor(
      SheehanWork_3cat,
      levels = c(1, 2, 3),
      labels = c("High", "Moderate", "Low"),
      ordered = TRUE
    ),
    SheehanSocial_3cat_ord = factor(
      SheehanSocial_3cat,
      levels = c(1, 2, 3),
      labels = c("High", "Moderate", "Low"),
      ordered = TRUE
    ),
    SheehanFamily_3cat_ord = factor(
      SheehanFamily_3cat,
      levels = c(1, 2, 3),
      labels = c("High", "Moderate", "Low"),
      ordered = TRUE
    )
  )

# -----------------------------
# Suicidality
# -----------------------------
GSRD <- GSRD %>%
  mutate(
    SuicidalRisk_bin = if_else(Suicidal.Risk == 0, 0L, 1L),
    SuicidalRisk_bin = factor(
      SuicidalRisk_bin,
      levels = c(0, 1),
      labels = c("No", "Yes")
    )
  )

# -----------------------------
# Number of depressive episodes
# -----------------------------
GSRD <- GSRD %>%
  mutate(
    N_dep_episodes_3cat = case_when(
      Number.of.Depressive.Episodes %in% 0:2  ~ 1L,
      Number.of.Depressive.Episodes %in% 3:5  ~ 2L,
      Number.of.Depressive.Episodes %in% 6:10 ~ 3L,
      TRUE ~ NA_integer_
    ),
    N_dep_episodes_3cat_ord = factor(
      N_dep_episodes_3cat,
      levels = c(1, 2, 3),
      labels = c("Low", "Moderate", "High"),
      ordered = TRUE
    )
  )

# -----------------------------
# Inner tension
# -----------------------------
GSRD <- GSRD %>%
  mutate(
    InnerTension_3cat = case_when(
      Inner.tension.MADRS.Retrospective %in% 0:2 ~ 1L,
      Inner.tension.MADRS.Retrospective %in% 3:4 ~ 2L,
      Inner.tension.MADRS.Retrospective %in% 5:6 ~ 3L,
      TRUE ~ NA_integer_
    ),
    InnerTension_3cat_ord = factor(
      InnerTension_3cat,
      levels = c(1, 2, 3),
      labels = c("Low", "Moderate", "High"),
      ordered = TRUE
    )
  )

# -----------------------------
# Resistance
# -----------------------------
# This cohort includes a single resistance phenotype.
# Replace 'Resistance_original' with the corresponding
# variable name if the outcome is stored under a different name.

if ("Resistance" %in% names(GSRD)) {
  GSRD <- GSRD %>%
    mutate(Resistance = as.numeric(Resistance))
}

# -----------------------------
# PRS residualisation on PCs and site variables
# -----------------------------
pc_vars   <- grep("^PC[1-9]$", names(GSRD), value = TRUE)
site_vars <- grep("^site[0-9]+$", names(GSRD), value = TRUE)
covars    <- c(pc_vars, site_vars)

prs_vars <- c("ADHD05", "DEPR07", "EDUC03", "SCZ", "BP", "BODY04")
prs_vars <- prs_vars[prs_vars %in% names(GSRD)]

for (v in prs_vars) {
  if (!is.numeric(GSRD[[v]])) {
    suppressWarnings(GSRD[[v]] <- as.numeric(as.character(GSRD[[v]])))
  }
}

for (v in prs_vars) {
  GSRD[[paste0(v, "_resid")]] <- make_resid_on_covars(GSRD, y = v, covars = covars)
}

# -----------------------------
# Save output
# -----------------------------
write.csv(GSRD, "GSRD_preprocessed.csv", row.names = FALSE, na = "")