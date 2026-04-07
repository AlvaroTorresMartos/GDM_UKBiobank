
# https://github.com/kosukeimai/mediation
# https://scholar.harvard.edu/dtingley/scholar_software/mediation
library(mediation)
library(dplyr)
library(readr)

# T2DM -----
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
metabolomics = metabolomics %>% 
  dplyr::filter(time_years > 1)
dim(metabolomics)
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
  ) %>% 
  # dplyr::filter(gdm_final == 1) %>% 
  droplevels()

# Time 
metabolomics = metabolomics %>% 
  dplyr::mutate(time_pregnant_baseline = age_at_recruitment - age_first_pregnancy_final)

# Metabolomics ------

# metabolomics =  metabolomics %>%
#   dplyr::mutate(diab_fam = na_if(diab_fam, "9"),
#                 diab_fam = as.factor(diab_fam))

metabolomics = metabolomics %>% 
  dplyr::filter(!is.na(diab_fam)) %>% 
  dplyr::filter(!is.na(live_births))

metabolomics$diab_fam %>% table(useNA = "ifany")

med_model = lm(z_scaled ~ gdm_final + age_first_pregnancy_final + ethnicity, 
               data = metabolomics)
summary(med_model)

out_model = glm(t2dm ~ z_scaled + age_first_pregnancy_final + 
                  ethnicity + townsend_q5 + education_group + 
                  diab_fam + bmi_group + smoking_status + 
                  medication_cholesterol + 
                  pa_rec + fruit_veg_rec + satfat_3 + 
                  alc_group + menopause + live_births + gdm_final
                , data = metabolomics, family = "binomial") 

summary(out_model)



set.seed(123456789)
mediation = mediate(model.m = med_model, model.y = out_model, 
                    mediator = "z_scaled", treat = "gdm_final", 
                    sims = 1000)

# ACME = Average Causal Mediation Effects (ACME)
# ADE = Average Direct Effects
# summary(mediation)
# Causal Mediation Analysis 
# 
# Quasi-Bayesian Confidence Intervals
# 
# Estimate 95% CI Lower 95% CI Upper p-value    
# ACME (control)            0.00675      0.00556         0.01  <2e-16 ***
#   ACME (treated)            0.01628      0.01211         0.02  <2e-16 ***
#   ADE (control)             0.04980      0.03353         0.07  <2e-16 ***
#   ADE (treated)             0.05933      0.03962         0.08  <2e-16 ***
#   Total Effect              0.06608      0.04583         0.09  <2e-16 ***
#   Prop. Mediated (control)  0.10247      0.07330         0.15  <2e-16 ***
#   Prop. Mediated (treated)  0.24646      0.20332         0.30  <2e-16 ***
#   ACME (average)            0.01152      0.00894         0.01  <2e-16 ***
#   ADE (average)             0.05457      0.03652         0.08  <2e-16 ***
#   Prop. Mediated (average)  0.17446      0.13875         0.22  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Sample Size Used: 194310 
# 
# 
# Simulations: 1000 
 

# plot(mediation) # , xlim = c(0, 0.15))

# medsens(mediation)

# saveRDS(mediation, "/mnt/project/results/4_mediation/full_adjustment_mediation_model_metabolomics.RDS")

# EFFECT EXCLUDING BMI, PHYSICAL ACTIVITY AND DIET 
# med_model = lm(z_scaled ~ gdm_final + age_first_pregnancy_final + ethnicity, 
#                data = metabolomics)
# summary(med_model)
# 
# out_model = glm(t2dm ~ z_scaled + age_first_pregnancy_final + 
#                   ethnicity + townsend_q5 + education_group + 
#                   diab_fam + smoking_status + 
#                   alc_group + menopause + gdm_final
#                 , data = metabolomics, family = "binomial") # binomial("probit")
# 
# summary(out_model)
# 
# 
# set.seed(123456789)
# mediation = mediate(model.m = med_model, model.y = out_model, 
#                     mediator = "z_scaled", treat = "gdm_final", 
#                     sims = 1000)
# 
# # ACME = Average Causal Mediation Effects (ACME)
# # ADE = Average Direct Effects
# summary(mediation)
# 
# 
# plot(mediation) #, xlim = c(0, 0.2))
# 


# Proteomics ----

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
proteomics = proteomics %>% 
  dplyr::filter(time_years > 1)

dim(proteomics)
# [1] 17890    26

# gdm = read_rds("/mnt/project/data/processed/2025_04_22_curated_dataset.RDS")
# 
# gdm = gdm %>% 
#   dplyr::select(eid, gdm_final) %>% 
#   dplyr::mutate(gdm_final = if_else(is.na(gdm_final), 0, gdm_final)) 

proteomics = proteomics %>% 
  # merge(gdm, by = "eid") %>% 
  dplyr::filter(ethnicity != 9) %>%
  # dplyr::filter(education_group != 3) %>% 
  # dplyr::filter(smoking_status != 9) %>% 
  # dplyr::filter(townsend_q5 != 9) %>% 
  dplyr::mutate(gdm_final = case_when(gdm_final == "GDM" ~ 1, 
                                      gdm_final == "NoGDM" ~ 0), 
                gdm_final = as.factor(gdm_final), 
  ) %>% 
  # dplyr::filter(gdm_final == 1) %>% 
  droplevels()

dim(proteomics)
# [1] 17841    26

proteomics = proteomics %>% 
  dplyr::filter(!is.na(diab_fam)) %>% 
  dplyr::filter(!is.na(live_births))

# mediaiton
med_model = lm(z_scaled ~ gdm_final + age_first_pregnancy_final + ethnicity, 
               data = proteomics)
summary(med_model)

out_model = glm(t2dm ~ z_scaled + gdm_final 
                + townsend_q5 + education_group +
                  diab_fam + bmi_group + smoking_status +
                  medication_cholesterol + 
                  pa_rec + fruit_veg_rec + satfat_3 +
                  alc_group + menopause + live_births + age_first_pregnancy_final + ethnicity
                , data = proteomics, family = "binomial") 


summary(out_model)


set.seed(123456789)
mediation = mediate(model.m = med_model, model.y = out_model, 
                    mediator = "z_scaled", treat = "gdm_final", 
                    sims = 1000)

# ACME = Average Causal Mediation Effects (ACME)
# ADE = Average Direct Effects
summary(mediation)

# Causal Mediation Analysis 
# 
# Quasi-Bayesian Confidence Intervals
# 
# Estimate 95% CI Lower 95% CI Upper p-value  
# ACME (control)           -0.000150    -0.004044         0.00   0.914  
# ACME (treated)           -0.000523    -0.009713         0.01   0.914  
# ADE (control)             0.044090    -0.004344         0.12   0.078 .
# ADE (treated)             0.043716    -0.004599         0.11   0.078 .
# Total Effect              0.043566    -0.003263         0.11   0.086 .
# Prop. Mediated (control) -0.004514    -0.445631         0.26   0.912  
# Prop. Mediated (treated) -0.010553    -0.504268         0.31   0.912  
# ACME (average)           -0.000337    -0.006885         0.01   0.914  
# ADE (average)             0.043903    -0.004472         0.12   0.078 .
# Prop. Mediated (average) -0.007533    -0.476276         0.29   0.912  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Sample Size Used: 17576 
# 
# 
# Simulations: 1000

# saveRDS(mediation, "./upload/full_adjustment_mediation_model_proteomics.RDS")

# Upload the files to RAP UK Bibank----
## Locate all the files generated in the upload folder 

# setwd("./upload")

# system("dx ls data/processed")

# system("pwd")

system("dx upload /home/rstudio-server/upload/*")