############################################################
# Project: Transdiagnostic predictors of treatment outcomes
# Cohort: STAR*D
# Script: Data preprocessing
############################################################

library(dplyr)
library(readr)
library(stringr)
library(tidyr)

options(stringsAsFactors = FALSE)

# -----------------------------
# Input folders
# -----------------------------
# Replace these with the local folders containing the STAR*D source files

enrollment_dir <- "path/to/STAR/Enrollment"
ivr_dir        <- "path/to/STAR/Level1_IVR"
prs_dir        <- "path/to/STAR/PRSs"
cov_dir        <- "path/to/STAR/phenotypes_cov"

# -----------------------------
# Main dataset
# -----------------------------
STAR <- read.csv("STAR.csv", stringsAsFactors = FALSE)

# -----------------------------
# Helpers
# -----------------------------
parse_num <- function(x) {
  suppressWarnings(readr::parse_number(as.character(x)))
}

residualize_prs_on_pcs <- function(df, prs_name, pcs, min_n = 50) {
  vars_needed <- c(prs_name, pcs)
  tmp <- df[, vars_needed, drop = FALSE]
  ok <- complete.cases(tmp)
  
  if (sum(ok) < min_n) {
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
# Demographics from DM
# -----------------------------
dm <- read.csv(file.path(enrollment_dir, "DM.csv"), stringsAsFactors = FALSE)

STAR <- dm %>%
  transmute(
    Patient.ID,
    Employment = Current.employment.status,
    Education = Years.of.schooling.completed,
    Married = Current.marital.status,
    N_relatives = Number.of.relatives.living.with.patient,
    N_friends = Number.of.friends.living.with.patient,
    Household_size = Total.number.of.persons.in.household,
    Living_conditions = case_when(
      !is.na(Total.number.of.persons.in.household) &
        Total.number.of.persons.in.household == 1 ~ "Lives_alone",
      !is.na(Total.number.of.persons.in.household) &
        Total.number.of.persons.in.household > 1 ~ "Lives_with_others",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    Education = ifelse(Education < 0, NA, Education),
    Education_3cat = case_when(
      Education <= 12 ~ "Low",
      Education <= 16 ~ "Medium",
      Education > 16  ~ "High",
      TRUE ~ NA_character_
    ),
    Education_3cat = factor(
      Education_3cat,
      levels = c("Low", "Medium", "High"),
      ordered = TRUE
    ),
    Employed = case_when(
      Employment %in% c("Full time", "Part time", "Self employed") ~ 1,
      Employment %in% c("Unempl, looking", "Unempl, not looking")  ~ 0,
      Employment %in% c("Retired, not working", "Unknown")         ~ NA_real_,
      TRUE ~ NA_real_
    )
  )

# -----------------------------
# Demographics from EL
# -----------------------------
EL <- read.csv(file.path(enrollment_dir, "EL.csv"), stringsAsFactors = FALSE)

EL_demo <- EL %>%
  transmute(
    Patient.ID,
    Age = Age.at.entry.to.study,
    Gender = Gender,
    Hispanic = `Is.the.patient.hispanic.`,
    Ethnicity = case_when(
      White == "Yes" ~ "White",
      `Black.or.African.American` == "Yes" ~ "Black",
      TRUE ~ NA_character_
    )
  )

STAR <- STAR %>%
  left_join(EL_demo, by = "Patient.ID")

# -----------------------------
# Clinical history from PHX
# -----------------------------
PHX <- read.csv(file.path(enrollment_dir, "PHX.csv"), stringsAsFactors = FALSE)

PHX_clin <- PHX %>%
  transmute(
    Patient.ID,
    Age_at_onset = ifelse(Age.at.onset.of.first.MDE >= 0, Age.at.onset.of.first.MDE, NA_integer_),
    N_previous_episodes = ifelse(Number.of.MDEs >= 0, Number.of.MDEs, NA_integer_),
    FH_depression = case_when(
      Family.hx.depression == "Yes" ~ 1,
      Family.hx.depression == "No"  ~ 0,
      TRUE ~ NA_real_
    ),
    FH_bipolar = case_when(
      Family.hx.bipolar.disorder == "Yes" ~ 1,
      Family.hx.bipolar.disorder == "No"  ~ 0,
      TRUE ~ NA_real_
    ),
    Family_History_mood = case_when(
      Family.hx.depression == "Yes" | Family.hx.bipolar.disorder == "Yes" ~ 1,
      Family.hx.depression == "No" & Family.hx.bipolar.disorder == "No"   ~ 0,
      TRUE ~ NA_real_
    ),
    FH_alcohol_abuse = case_when(
      Family.hx.alcohol.abuse == "Yes" ~ 1,
      Family.hx.alcohol.abuse == "No"  ~ 0,
      TRUE ~ NA_real_
    ),
    FH_drug_abuse = case_when(
      Family.hx.drug.abuse == "Yes" ~ 1,
      Family.hx.drug.abuse == "No"  ~ 0,
      TRUE ~ NA_real_
    ),
    FH_suicide = case_when(
      Family.hx.suicide == "Yes" ~ 1,
      Family.hx.suicide == "No"  ~ 0,
      TRUE ~ NA_real_
    )
  )

STAR <- STAR %>%
  left_join(PHX_clin, by = "Patient.ID") %>%
  mutate(
    Age_at_onset = ifelse(Age_at_onset <= 10, NA, Age_at_onset),
    N_previous_episodes_log = case_when(
      is.na(N_previous_episodes) ~ NA_real_,
      N_previous_episodes > 0    ~ log(N_previous_episodes),
      TRUE ~ NA_real_
    )
  )

# -----------------------------
# HRSD baseline variables
# -----------------------------
HRSD <- read.csv(file.path(enrollment_dir, "HRSD.csv"), stringsAsFactors = FALSE)

HRSD_clin <- HRSD %>%
  transmute(
    Patient.ID,
    HRSD_total = HRS.total.score..recorded.,
    HRSD_total_crc = HRS.total.score..CRC.,
    HRSD_suicide = HRS.Suicide,
    Suicidal_ideation_bin = case_when(
      is.na(HRS.Suicide) ~ NA_real_,
      HRS.Suicide == "Absent" ~ 0,
      TRUE ~ 1
    ),
    Functional_Impairment_3cat = case_when(
      HRS.Work.and.interests %in% c("", "-2") ~ NA_character_,
      HRS.Work.and.interests %in% c("No disturbance", "Mild") ~ "Low",
      HRS.Work.and.interests == "Moderate" ~ "Moderate",
      HRS.Work.and.interests %in% c("Moderately severe", "Severe") ~ "High",
      TRUE ~ NA_character_
    ),
    Anxiety_3cat = case_when(
      HRS.Psychic.anxiety %in% c("", "-2") ~ NA_character_,
      HRS.Psychic.anxiety %in% c("Absent", "Mild") ~ "Low",
      HRS.Psychic.anxiety == "Moderate" ~ "Moderate",
      HRS.Psychic.anxiety %in% c("Moderately severe", "Severe") ~ "High",
      TRUE ~ NA_character_
    ),
    Depressive_symptoms_3cat = case_when(
      HRS.Depressed.mood %in% c("", "-2") ~ NA_character_,
      HRS.Depressed.mood %in% c("Not depressed", "Mild") ~ "Low",
      HRS.Depressed.mood == "Moderate" ~ "Moderate",
      HRS.Depressed.mood %in% c("Moderately severe", "Severe") ~ "High",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    Functional_Impairment_3cat = factor(
      Functional_Impairment_3cat,
      levels = c("Low", "Moderate", "High"),
      ordered = TRUE
    ),
    Anxiety_3cat = factor(
      Anxiety_3cat,
      levels = c("Low", "Moderate", "High"),
      ordered = TRUE
    ),
    Depressive_symptoms_3cat = factor(
      Depressive_symptoms_3cat,
      levels = c("Low", "Moderate", "High"),
      ordered = TRUE
    )
  )

STAR <- STAR %>%
  mutate(Patient.ID = as.integer(Patient.ID)) %>%
  left_join(HRSD_clin %>% mutate(Patient.ID = as.integer(Patient.ID)), by = "Patient.ID")

# -----------------------------
# Anxiety and psychiatric comorbidity from PHX
# -----------------------------
PHX_hx <- PHX %>%
  rowwise() %>%
  mutate(
    HXPAN = {
      vals <- c_across(c(`AxI..Panic.w..agoraphobia`, `AxI..Panic.w.o.agoraphobia`))
      vals[vals %in% c("-2", "")] <- NA
      if (all(is.na(vals))) NA_real_ else if (any(vals == "Yes", na.rm = TRUE)) 1 else 0
    },
    HXGAD = {
      vals <- c_across(c(`AxI..Generalized.anxiety`))
      vals[vals %in% c("-2", "")] <- NA
      if (all(is.na(vals))) NA_real_ else if (any(vals == "Yes", na.rm = TRUE)) 1 else 0
    },
    HXSOC = {
      vals <- c_across(c(`AxI..Social.phobia`))
      vals[vals %in% c("-2", "")] <- NA
      if (all(is.na(vals))) NA_real_ else if (any(vals == "Yes", na.rm = TRUE)) 1 else 0
    },
    HXOCD = {
      vals <- c_across(c(`AxI..Obsessive.compulsive`))
      vals[vals %in% c("-2", "")] <- NA
      if (all(is.na(vals))) NA_real_ else if (any(vals == "Yes", na.rm = TRUE)) 1 else 0
    },
    HXPTSD = {
      vals <- c_across(c(`AxI..Posttraumatic.stress`))
      vals[vals %in% c("-2", "")] <- NA
      if (all(is.na(vals))) NA_real_ else if (any(vals == "Yes", na.rm = TRUE)) 1 else 0
    }
  ) %>%
  ungroup() %>%
  transmute(Patient.ID, HXPAN, HXGAD, HXSOC, HXOCD, HXPTSD)

STAR <- STAR %>%
  left_join(PHX_hx, by = "Patient.ID") %>%
  rowwise() %>%
  mutate(
    Anxiety = {
      vals <- c_across(c(HXPAN, HXGAD, HXSOC))
      if (any(vals == 1, na.rm = TRUE)) 1 else if (any(vals == 0, na.rm = TRUE)) 0 else NA_real_
    },
    Psychiatric_Comorbidity = {
      vals <- c_across(c(HXPAN, HXGAD, HXSOC, HXOCD, HXPTSD))
      if (any(vals == 1, na.rm = TRUE)) 1 else if (any(vals == 0, na.rm = TRUE)) 0 else NA_real_
    }
  ) %>%
  ungroup()

# -----------------------------
# Substance use problem from PDS
# -----------------------------
PDS <- read.csv(file.path(enrollment_dir, "PDS.csv"), stringsAsFactors = FALSE)

pds_sub_vars <- c(
  "Think.drink.too.much",
  "Think.had.alcohol.problem",
  "Think.about.cutting.down.on.drinking",
  "Think.using.drugs.too.much",
  "Think.had.a.drug.problem",
  "Think.about.cutting.down.on.drug.use"
)

PDS_substance <- PDS %>%
  rowwise() %>%
  mutate(
    Substance_Use_Problem = {
      vals <- c_across(all_of(pds_sub_vars))
      vals[vals %in% c("-2", "")] <- NA
      if (all(is.na(vals))) NA_real_ else if (any(vals == "Yes", na.rm = TRUE)) 1 else 0
    }
  ) %>%
  ungroup() %>%
  transmute(Patient.ID, Substance_Use_Problem)

STAR <- STAR %>%
  left_join(PDS_substance, by = "Patient.ID")

# -----------------------------
# WSAS baseline
# -----------------------------
IVR <- read.csv(file.path(ivr_dir, "IVR_assmt.csv"), stringsAsFactors = FALSE)

WSAS_baseline <- IVR %>%
  filter(
    Treatment.level == 1,
    CallType == "Base"
  ) %>%
  arrange(Patient.ID, Date.of.call) %>%
  group_by(Patient.ID) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    Patient.ID,
    WSAS_baseline = Sum.of.WSAS.items.1.5
  )

STAR <- STAR %>%
  left_join(WSAS_baseline, by = "Patient.ID")

# -----------------------------
# Psychotic symptoms from PDS
# -----------------------------
im_vars <- grep("^IM\\.", names(PDS), value = TRUE)

PDS_psy <- PDS %>%
  rowwise() %>%
  mutate(
    Psychotic_symptoms_bin = {
      vals <- c_across(all_of(im_vars))
      vals[vals %in% c("-2", "")] <- NA
      if (all(is.na(vals))) NA_real_ else if (any(vals == "Yes", na.rm = TRUE)) 1 else 0
    }
  ) %>%
  ungroup() %>%
  transmute(Patient.ID, Psychotic_symptoms_bin)

STAR <- STAR %>%
  left_join(PDS_psy, by = "Patient.ID")

# -----------------------------
# Social withdrawal from PDS
# -----------------------------
sw_vars <- c(
  "EM.Avoid.situ.b.c.afraid.embarrass.self",
  "EM.Avoid.situ.b.c.made.you.anxious",
  "FR.Avoid.situ.b.c.afraid.of.anx.attack",
  "FR.Avoid.sit.b.c.made.you.anxious"
)

PDS_sw <- PDS %>%
  rowwise() %>%
  mutate(
    Social_Withdrawal_bin = {
      vals <- c_across(all_of(sw_vars))
      vals[vals %in% c("-2", "")] <- NA
      if (all(is.na(vals))) NA_real_ else if (any(vals == "Yes", na.rm = TRUE)) 1 else 0
    }
  ) %>%
  ungroup() %>%
  transmute(Patient.ID, Social_Withdrawal_bin)

STAR <- STAR %>%
  left_join(PDS_sw, by = "Patient.ID")

# -----------------------------
# Genetic ID mapping
# -----------------------------
cellline2 <- read.table(
  file.path(cov_dir, "STAR_D_IDs_by_cellline_2.txt"),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
) %>%
  rename(genID = Lab_ID..SD.Name.)

STAR <- STAR %>%
  select(-any_of("genID")) %>%
  left_join(
    cellline2 %>% select(Patient.ID, genID),
    by = "Patient.ID"
  )

# -----------------------------
# PRSs
# -----------------------------
PRS <- read.table(
  file.path(prs_dir, "all_STAR_PRSs_standard_clinID"),
  header = TRUE,
  sep = "",
  stringsAsFactors = FALSE
)

STAR <- STAR %>%
  left_join(PRS, by = "Patient.ID")

# -----------------------------
# Covariates: centres + PCs
# -----------------------------
cov <- read.table(
  file.path(cov_dir, "cov_STARD"),
  header = FALSE,
  sep = "",
  stringsAsFactors = FALSE
)

colnames(cov) <- c(
  "genID_FID",
  "genID_IID",
  paste0("centre_", 1:11),
  "age",
  "baseline_severity",
  paste0("PC", 1:10)
)

cov_key <- cov %>%
  transmute(
    genID = genID_FID,
    across(starts_with("centre_")),
    age,
    baseline_severity,
    across(starts_with("PC"))
  ) %>%
  distinct(genID, .keep_all = TRUE)

STAR <- STAR %>%
  left_join(cov_key, by = "genID")

# -----------------------------
# PRS residualisation
# -----------------------------
prs_vars <- c(
  "ADHD05_PRS",
  "BODY04_PRS",
  "BP_PRS",
  "DEPR07_PRS",
  "EDUC03_PRS",
  "SCZ_PRS"
)

pc_vars <- paste0("PC", 1:10)

for (prs in prs_vars) {
  STAR <- residualize_prs_on_pcs(STAR, prs_name = prs, pcs = pc_vars, min_n = 50)
}

# -----------------------------
# Outcomes
# -----------------------------
trd_vs_resp_key <- read.table(
  file.path(cov_dir, "TRD_vs_resp"),
  header = FALSE,
  sep = "",
  stringsAsFactors = FALSE
) %>%
  transmute(
    genID = as.character(V1),
    TRD_vs_resp = as.integer(V3)
  ) %>%
  distinct(genID, .keep_all = TRUE)

trd_nonresp_vs_resp_key <- read.table(
  file.path(cov_dir, "TRD_non_resp_vs_resp"),
  header = FALSE,
  sep = "",
  stringsAsFactors = FALSE
) %>%
  transmute(
    genID = as.character(V1),
    TRD_non_resp_vs_resp = as.integer(V3)
  ) %>%
  distinct(genID, .keep_all = TRUE)

STAR <- STAR %>%
  left_join(trd_vs_resp_key, by = "genID") %>%
  left_join(trd_nonresp_vs_resp_key, by = "genID")

# Remission outcome
if ("QIDS_w12_rem" %in% names(STAR)) {
  STAR <- STAR %>%
    mutate(
      QIDS_w12_rem = case_when(
        QIDS_w12_rem == 1 ~ 1L,
        QIDS_w12_rem == 2 ~ 0L,
        TRUE ~ NA_integer_
      )
    )
}

# Resistance outcome
if ("TRD_vs_resp" %in% names(STAR)) {
  STAR <- STAR %>%
    mutate(
      TRD_vs_resp = case_when(
        TRD_vs_resp == 2 ~ 1L,
        TRD_vs_resp == 1 ~ 0L,
        TRUE ~ NA_integer_
      )
    )
}

# -----------------------------
# Save output
# -----------------------------
write.csv(STAR, "STAR_preprocessed.csv", row.names = FALSE)