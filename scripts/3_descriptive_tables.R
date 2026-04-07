

# Load packages ------
if (!require("pacman")) install.packages("pacman")

pacman::p_load(readr, dplyr, descr, gtsummary)


# Import data -----

# stud_pop = read_rds("/mnt/project/data/processed/2025_04_22_curated_dataset.RDS")
# hes = read_rds("/mnt/project/data/processed/HES.RDS")
# demographics = read_rds("/mnt/project/data/processed/2025_04_13_covariates.RDS")
# 
# stud_pop = stud_pop %>% 
#   dplyr::select(eid, age_first_pregnancy_final, gdm_final)
# 
# hes = hes %>% 
#   dplyr::select(eid, t2dm)
# 
# demographics = demographics %>% 
#   dplyr::select(eid, age_5y_group:htn_fam, `2734_0_0`) %>% 
#   dplyr::mutate(`2734_0_0` = if_else(`2734_0_0` %in% c(0, -3), NA_real_, `2734_0_0`)) %>% 
#   dplyr::rename(live_births = `2734_0_0`)
# 
# master_file = dplyr::left_join(stud_pop, demographics, by = "eid") %>% 
#   dplyr::left_join(hes, by = "eid")
# 
# rm(stud_pop,demographics, hes)
# 
# master_file = master_file %>% 
#   dplyr::mutate(diab_fam = case_when(is.na(diab_fam) ~ 9, 
#                                      TRUE ~ diab_fam))
# 
# master_file$age_first_pregnancy_final %>% descr::freq()
# master_file$ethnicity  %>% descr::freq()
# master_file$townsend_q5  %>% descr::freq()
# master_file$education_group  %>% descr::freq()
# master_file$diab_fam  %>% descr::freq()
# master_file$bmi_group  %>% descr::freq()
# master_file$menopause  %>% descr::freq()
# master_file$live_births %>%  descr::freq()
# master_file$t2dm  %>% descr::freq()
# master_file$smoking_status  %>% descr::freq()
# master_file$pa_rec  %>% descr::freq()
# master_file$fruit_veg_rec  %>% descr::freq()
# master_file$satfat_3  %>% descr::freq() # no es exactamente el mismo pero vale 
# master_file$alc_group  %>% descr::freq()

# Load packages ------
if (!require("pacman")) install.packages("pacman")

pacman::p_load(readr, dplyr,  survival, broom, ggplot2, forcats)

# metabolomics -----

metabolomics = read_rds("/mnt/project/data/processed/2025_12_15_metabolomic_dataset_limma_elasticnet.RDS")

colnames(metabolomics)

metabolomics = metabolomics %>% 
  dplyr::select(n_eid:ethnicity, gdm_final)

metabolomics = metabolomics %>% 
  dplyr::mutate(across(where(is.factor), droplevels))

signature = read_rds("/mnt/project/results/2_elasticnet/metabolic_signature_GDM.RDS")
signature = signature %>% 
  dplyr::select(eid, z_scaled) %>% 
  dplyr::rename(n_eid = eid)

hes = read_rds("/mnt/project/data/processed/HES.RDS")
demographics = read_csv("/mnt/project/data/raw/Demographics.csv")
covariates = read_rds("/mnt/project/data/processed/2025_12_10_covariates.RDS")


covariates = covariates %>%  
  dplyr::select(eid, `21022_0_0`, country, 
                pa_rec, satfat_3, alc_group, bmi_group, diab_fam, 
                education_group, smoking_status, townsend_q5, menopause, 
                fruit_veg_rec, medication_cholesterol
  ) %>% 
  dplyr::rename(age_at_recruitment = `21022_0_0`)

demographics = demographics %>% 
  dplyr::select(eid, `53-0.0`, `2734-0.0`) %>% 
  dplyr::mutate(`53-0.0` = as.Date(`53-0.0`), 
                `2734-0.0` = if_else(`2734-0.0` %in% c(0, -3), NA_real_, `2734-0.0`)) %>% 
  dplyr::rename(live_births = `2734-0.0`)

hes = hes %>% 
  dplyr::select(eid, date_t2dm, t2dm) %>% 
  dplyr::left_join(demographics, by = "eid")

metabolomics = metabolomics %>% 
  merge(signature, by = "n_eid") %>% 
  dplyr::rename(eid = n_eid) %>% 
  dplyr::left_join(hes, by = "eid") %>% 
  dplyr::left_join(covariates, by = "eid") %>% 
  dplyr::mutate(
    date_t2dm_last_update = case_when(
      t2dm == 1 ~ date_t2dm,
      t2dm == 0 & country == 1 ~ as.Date("2022-08-31"), 
      t2dm == 0 & country == 2 ~ as.Date("2022-05-31"), 
      t2dm == 0 & country == 3 ~ as.Date("2022-10-31")), 
    time_days = as.numeric(date_t2dm_last_update - `53-0.0`), 
    time_years = time_days/365.25, 
    age_out = age_at_recruitment + time_years)



# metabolomics = metabolomics %>%
#   dplyr::filter(menopause != 9)
# 
# metabolomics = metabolomics %>% 
#   dplyr::mutate(across(where(is.factor), droplevels))

dim(metabolomics)
# 197263     26
# metabolomics = metabolomics %>% 
#   dplyr::filter(time_years > 1)
# dim(metabolomics)
# 197199    26

# gdm = read_rds("/mnt/project/data/processed/2025_04_22_curated_dataset.RDS")
# 
# gdm = gdm %>% 
#   dplyr::select(eid, gdm_final) %>% 
#   dplyr::mutate(gdm_final = if_else(is.na(gdm_final), 0, gdm_final)) 

metabolomics = metabolomics %>% 
  # merge(gdm, by = "eid") %>% 
  # dplyr::filter(ethnicity != 9) %>%
  # dplyr::filter(education_group != 9) %>% 
  # dplyr::filter(smoking_status != 9) %>% 
  # dplyr::filter(townsend_q5 != 9) %>% 
  dplyr::mutate(gdm_final = case_when(gdm_final == "GDM" ~ 1, 
                                      gdm_final == "NoGDM" ~ 0), 
                gdm_final = as.factor(gdm_final), 
                diab_fam = as.numeric(diab_fam), 
                diab_fam = if_else(is.na(diab_fam), 9, diab_fam), 
                diab_fam = as.factor(diab_fam)
  ) %>% 
  # dplyr::filter(gdm_final == 1) %>% 
  droplevels()

# Time 
metabolomics = metabolomics %>% 
  dplyr::mutate(time_pregnant_baseline = age_at_recruitment - age_first_pregnancy_final)

master_file2 = metabolomics %>% 
  dplyr::select(gdm_final, 
                age_at_recruitment, 
                age_first_pregnancy_final, ethnicity, townsend_q5, 
                education_group, diab_fam, bmi_group, medication_cholesterol, 
                menopause, live_births, t2dm, smoking_status, pa_rec, 
                fruit_veg_rec, satfat_3, alc_group) #%>% 
  # dplyr::mutate(gdm_final = case_when(is.na(gdm_final) ~ 0, TRUE ~ gdm_final))


master_file2 %>% 
  # dplyr::select(-c(eid)) %>% 
  dplyr::mutate(medication_cholesterol = as.factor(medication_cholesterol)) %>% 
  gtsummary::tbl_summary(by = "gdm_final", 
                         type = all_continuous2() ~ "continuous2", 
                         statistic = all_continuous() ~ c("{mean} ({sd}) \  \ {median} ({p25}, {p75})"), 
                         digits = list(
                           all_continuous() ~ 2,
                           all_categorical() ~ c(0,2)
                         ), 
                         missing_text = "Unknown",
                         ) %>% 
  gtsummary::add_p(test = everything() ~ "chisq.test")


# metabolomics = read_rds("/mnt/project/data/processed/2025_04_22_metabolomics_dataset_elasticnet.RDS") %>% 
#   dplyr::select(eid)
# 
# proteomics = read_rds("/mnt/project/data/processed/2025_09_10_proteomic_dataste_limmma.RDS") %>% 
#   dplyr::select(eid)
# 
# master_file2 %>% 
#   dplyr::filter(eid %in% metabolomics$eid) %>% 
#   dplyr::select(-c(eid)) %>% 
#   gtsummary::tbl_summary(by = "gdm_final", 
#                          type = all_continuous2() ~ "continuous2", 
#                          statistic = all_continuous() ~ c("{mean} ({sd}) \  \ {median} ({p25}, {p75})"), 
#                          digits = list(
#                            all_continuous() ~ 2,
#                            all_categorical() ~ c(0,2)
#                          ), 
#                          missing_text = "Unknown",
#   ) %>% 
#   gtsummary::add_p(test = everything() ~ "chisq.test")
# 
# 
# master_file2 %>% 
#   dplyr::filter(eid %in% proteomics$eid) %>% 
#   dplyr::select(-c(eid)) %>% 
#   gtsummary::tbl_summary(by = "gdm_final", 
#                          type = all_continuous2() ~ "continuous2", 
#                          statistic = all_continuous() ~ c("{mean} ({sd}) \  \ {median} ({p25}, {p75})"), 
#                          digits = list(
#                            all_continuous() ~ 2,
#                            all_categorical() ~ c(0,2)
#                          ), 
#                          missing_text = "Unknown",
#   ) %>% 
#   gtsummary::add_p(test = everything() ~ "chisq.test")


# proteomics -----
proteomics = read_rds("/mnt/project/data/processed/2025_12_17_proteomic_dataset_limmma_elasticnet.RDS")

nrow(proteomics)
# 17895

colnames(proteomics)

proteomics = proteomics %>% 
  dplyr::select(n_eid:ethnicity, gdm_final)

proteomics = proteomics %>% 
  dplyr::mutate(across(where(is.factor), droplevels))

signature = read_rds("/mnt/project/results/2_elasticnet/proteomic_signature.RDS")
signature = signature %>% 
  dplyr::select(eid, z_scaled) %>% 
  dplyr::rename(n_eid = eid)

hes = read_rds("/mnt/project/data/processed/HES.RDS")
demographics = read_csv("/mnt/project/data/raw/Demographics.csv")
covariates = read_rds("/mnt/project/data/processed/2025_12_10_covariates.RDS")

covariates = covariates %>%  
  dplyr::select(eid, `21022_0_0`, country, 
                pa_rec, satfat_3, alc_group, bmi_group, diab_fam, 
                education_group, smoking_status, townsend_q5, menopause, 
                fruit_veg_rec, medication_cholesterol
  ) %>% 
  dplyr::rename(age_at_recruitment = `21022_0_0`)

demographics = demographics %>% 
  dplyr::select(eid, `53-0.0`, `2734-0.0`) %>% 
  dplyr::mutate(`53-0.0` = as.Date(`53-0.0`), 
                `2734-0.0` = if_else(`2734-0.0` %in% c(0, -3), NA_real_, `2734-0.0`)) %>% 
  dplyr::rename(live_births = `2734-0.0`)

hes = hes %>% 
  dplyr::select(eid, date_t2dm, t2dm) %>% 
  dplyr::left_join(demographics, by = "eid")

proteomics = proteomics %>% 
  merge(signature, by = "n_eid") %>% 
  dplyr::rename(eid = n_eid) %>% 
  dplyr::left_join(hes, by = "eid") %>% 
  dplyr::left_join(covariates, by = "eid") %>% 
  dplyr::mutate(
    date_t2dm_last_update = case_when(
      t2dm == 1 ~ date_t2dm,
      t2dm == 0 & country == 1 ~ as.Date("2022-08-31"), 
      t2dm == 0 & country == 2 ~ as.Date("2022-05-31"), 
      t2dm == 0 & country == 3 ~ as.Date("2022-10-31")), 
    time_days = as.numeric(date_t2dm_last_update - `53-0.0`), 
    time_years = time_days/365.25, 
    age_out = age_at_recruitment + time_years)

# proteomics$ethnicity %>% table(useNA = "ifany")
# proteomics$townsend_q5 %>% table(useNA = "ifany")

# proteomics = proteomics %>% 
#   dplyr::filter(ethnicity != 9) %>% 
#   dplyr::filter(townsend_q5 != 9) 

# proteomics = proteomics %>% 
#   dplyr::mutate(across(where(is.factor), droplevels))


dim(proteomics)
# [1] 17895    26
# proteomics = proteomics %>% 
#   dplyr::filter(time_years > 1)

dim(proteomics)
# [1] 17890    26

# gdm = read_rds("/mnt/project/data/processed/2025_04_22_curated_dataset.RDS")
# 
# gdm = gdm %>% 
#   dplyr::select(eid, gdm_final) %>% 
#   dplyr::mutate(gdm_final = if_else(is.na(gdm_final), 0, gdm_final)) 

proteomics = proteomics %>% 
  # merge(gdm, by = "eid") %>% 
  # dplyr::filter(ethnicity != 9) %>%
  # dplyr::filter(education_group != 3) %>% 
  # dplyr::filter(smoking_status != 9) %>% 
  # dplyr::filter(townsend_q5 != 9) %>% 
  dplyr::mutate(gdm_final = case_when(gdm_final == "GDM" ~ 1, 
                                      gdm_final == "NoGDM" ~ 0), 
                gdm_final = as.factor(gdm_final), 
                diab_fam = as.numeric(diab_fam), 
                diab_fam = if_else(is.na(diab_fam), 9, diab_fam), 
                diab_fam = as.factor(diab_fam)
  ) %>% 
  # dplyr::filter(gdm_final == 1) %>% 
  droplevels()

dim(proteomics)
# [1] 17841    26

# Time 
proteomics = proteomics %>% 
  dplyr::mutate(time_pregnant_baseline = age_at_recruitment - age_first_pregnancy_final)

master_file2 = proteomics %>% 
  dplyr::select(gdm_final, 
                age_at_recruitment, 
                age_first_pregnancy_final, ethnicity, townsend_q5, 
                education_group, diab_fam, bmi_group, medication_cholesterol, 
                menopause, live_births, t2dm, smoking_status, pa_rec, 
                fruit_veg_rec, satfat_3, alc_group) #%>% 
# dplyr::mutate(gdm_final = case_when(is.na(gdm_final) ~ 0, TRUE ~ gdm_final))


master_file2 %>% 
  # dplyr::select(-c(eid)) %>% 
  dplyr::mutate(medication_cholesterol = as.factor(medication_cholesterol)) %>% 
  gtsummary::tbl_summary(by = "gdm_final", 
                         type = all_continuous2() ~ "continuous2", 
                         statistic = all_continuous() ~ c("{mean} ({sd}) \  \ {median} ({p25}, {p75})"), 
                         digits = list(
                           all_continuous() ~ 2,
                           all_categorical() ~ c(0,2)
                         ), 
                         missing_text = "Unknown",
  ) %>% 
  gtsummary::add_p(test = everything() ~ "chisq.test")
