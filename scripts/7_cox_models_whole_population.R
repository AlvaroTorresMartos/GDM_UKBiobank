
# Key notes -----
# https://stats.stackexchange.com/questions/317336/interpreting-r-coxph-cox-zph
# https://stats.stackexchange.com/questions/144923/extended-cox-model-and-cox-zph/238964#238964
# https://stats.stackexchange.com/questions/601592/violations-of-the-proportional-hazards-assumption-in-cox-how-to-use-strata-r

# Load packages ------
if (!require("pacman")) install.packages("pacman")

pacman::p_load(readr, dplyr,  survival, broom, ggplot2, forcats)

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
                diab_fam = as.numeric(diab_fam), 
                diab_fam = if_else(is.na(diab_fam), 9, diab_fam), 
                diab_fam = as.factor(diab_fam)
  ) %>% 
  # dplyr::filter(gdm_final == 1) %>% 
  droplevels()

# Time 
metabolomics = metabolomics %>% 
  dplyr::mutate(time_pregnant_baseline = age_at_recruitment - age_first_pregnancy_final)

surv_obj = with(metabolomics, Surv(time = age_at_recruitment,
                                   time2 = age_out,
                                   event = t2dm == 1))

# cox_metabolomics = coxph(surv_obj ~ z_scaled + 
#                            ethnicity + pa_rec + live_births +
#                            education_group + fruit_veg_rec + 
#                            medication_cholesterol +
#                            alc_group + bmi_group + 
#                            smoking_status + townsend_q5 + 
#                            menopause + diab_fam + 
#                            age_first_pregnancy_final + satfat_3, 
#                          data = metabolomics)

cox_metabolomics = coxph(surv_obj ~ z_scaled + 
                           ethnicity + pa_rec + live_births +
                           education_group + fruit_veg_rec + 
                           medication_cholesterol +
                           strata(alc_group) + strata(bmi_group) + 
                           strata(smoking_status) + strata(townsend_q5) + 
                           strata(menopause) + strata(diab_fam) + strata(satfat_3) +
                           strata(age_first_pregnancy_final), 
                         data = metabolomics)



# Test the Proportional Hazards Assumption of a Cox Regression                 
cox.zph(cox_metabolomics)
# chisq df    p
# z_scaled                0.7595  1 0.38
# ethnicity               5.2913  4 0.26
# pa_rec                  0.0852  1 0.77
# live_births             0.8893  1 0.35
# education_group         1.1343  1 0.29
# fruit_veg_rec           1.6282  1 0.20
# medication_cholesterol  0.0457  1 0.83
# GLOBAL                 10.8493 10 0.37
              
summary(cox_metabolomics)
#                         coef exp(coef)  se(coef)      z Pr(>|z|)    
# z_scaled                0.597728  1.817984  0.025826 23.144  < 2e-16 ***
                 
## Tidy general dataframe ----
cox_metabolomics = broom::tidy(cox_metabolomics, 
                               exponentiate = TRUE, conf.int = TRUE)

## Quantile dataframes and cox models ----

metabolomics = metabolomics %>% 
  dplyr::mutate(quantile_met = ntile(z_scaled, n = 5), 
                quantile_met = factor(quantile_met)
                )



surv_obj = with(metabolomics, Surv(time = age_at_recruitment,
                                   time2 = age_out,
                                   event = t2dm == 1))

cox_metabolomics2 = coxph(surv_obj ~ quantile_met + 
                           ethnicity + pa_rec + live_births +
                           education_group + fruit_veg_rec + 
                           medication_cholesterol +
                           strata(alc_group) + strata(bmi_group) + 
                           strata(smoking_status) + strata(townsend_q5) + 
                           strata(menopause) + strata(diab_fam) + 
                           strata(age_first_pregnancy_final) + strata(satfat_3), 
                         data = metabolomics)


cox.zph(cox_metabolomics2)
# chisq df    p
# quantile_met           1.88e+00  4 0.76
# ethnicity              4.89e+00  4 0.30
# pa_rec                 3.51e-02  1 0.85
# live_births            1.12e+00  1 0.29
# education_group        1.43e+00  1 0.23
# fruit_veg_rec          1.51e+00  1 0.22
# medication_cholesterol 1.99e-04  1 0.99
# GLOBAL                 1.11e+01 13 0.60

summary(cox_metabolomics2)
#                           coef  exp(coef)   se(coef)      z Pr(>|z|)    
# quantile_met2           0.2580952  1.2944621  0.0690029  3.740 0.000184 ***
#   quantile_met3           0.6333375  1.8838876  0.0677220  9.352  < 2e-16 ***
#   quantile_met4           0.9479250  2.5803500  0.0688953 13.759  < 2e-16 ***
#   quantile_met5           1.5104020  4.5285509  0.0741293 20.375  < 2e-16 ***



## Tidy general dataframe quantile----
cox_metabolomics2 = broom::tidy(cox_metabolomics2, 
                               exponentiate = TRUE, conf.int = TRUE)


## More information ----
info_metabolomics = metabolomics %>% 
  dplyr::summarise(
    count_total = n(), 
    count_t2dm = sum(t2dm == 1, na.rm = TRUE), 
    perc_t2dm = 100 * mean(t2dm == 1, na.rm = TRUE)
  )

info_metabolomics2 = metabolomics %>% 
  dplyr::group_by(quantile_met) %>% 
  dplyr::summarise(
    count_total = n(), 
    count_t2dm = sum(t2dm == 1, na.rm = TRUE), 
    perc_t2dm = 100 * mean(t2dm == 1, na.rm = TRUE)
  ) %>% 
  dplyr::select(-c(quantile_met))



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
                diab_fam = as.numeric(diab_fam), 
                diab_fam = if_else(is.na(diab_fam), 9, diab_fam), 
                diab_fam = as.factor(diab_fam)
  ) %>% 
  # dplyr::filter(gdm_final == 1) %>% 
  droplevels()

dim(proteomics)
# [1] 17841    26

surv_obj = with(proteomics, Surv(time = age_at_recruitment,
                                   time2 = age_out,
                                   event = t2dm == 1))

cox_proteomics = coxph(surv_obj ~ z_scaled + age_first_pregnancy_final + 
                         ethnicity + townsend_q5 + education_group  + 
                         medication_cholesterol + 
                         diab_fam + smoking_status + 
                         pa_rec + fruit_veg_rec + satfat_3 + 
                         alc_group + menopause +  live_births + 
                         strata(bmi_group)
                       , data = proteomics)

# cox_proteomics = coxph(surv_obj ~ z_scaled + age_first_pregnancy_final + 
#                ethnicity + townsend_q5 + education_group  + 
#                diab_fam + bmi_group + smoking_status + 
#                pa_rec + fruit_veg_rec + satfat_3 + 
#                alc_group + menopause +  live_births 
#             , data = proteomics)

cox.zph(cox_proteomics)
# chisq df    p
# z_scaled                  0.00577  1 0.94
# age_first_pregnancy_final 2.16069  1 0.14
# ethnicity                 1.85090  4 0.76
# townsend_q5               0.22949  1 0.63
# education_group           1.39455  1 0.24
# diab_fam                  0.64917  1 0.42
# smoking_status            1.07684  1 0.30
# pa_rec                    0.02006  1 0.89
# fruit_veg_rec             0.18875  1 0.66
# satfat_3                  0.08695  1 0.77
# alc_group                 0.19853  1 0.66
# menopause                 0.49919  1 0.48
# live_births               0.57258  1 0.45
# GLOBAL                    7.71980 16 0.96
summary(cox_proteomics)
#                           coef  exp(coef)   se(coef)      z Pr(>|z|)    
# z_scaled                   2.611e-01  1.298e+00  7.576e-02  3.447 0.000567 ***


## Tidy general dataframe ----
cox_proteomics = broom::tidy(cox_proteomics, 
                               exponentiate = TRUE, conf.int = TRUE)



## Quantile dataframes and cox models ----

proteomics = proteomics %>% 
  dplyr::mutate(quantile_met = ntile(z_scaled, n = 5), 
                quantile_met = factor(quantile_met)
  )

surv_obj = with(proteomics, Surv(time = age_at_recruitment,
                                 time2 = age_out,
                                 event = t2dm == 1))

cox_proteomics2 = coxph(surv_obj ~ quantile_met + age_first_pregnancy_final + 
                            ethnicity + townsend_q5 + education_group + 
                            medication_cholesterol +  
                            diab_fam + smoking_status + 
                            pa_rec + fruit_veg_rec + satfat_3 + 
                            alc_group + menopause +  live_births + 
                          strata(bmi_group)
                          , data = proteomics)

# cox_proteomics2 = coxph(surv_obj ~ quantile_met + age_first_pregnancy_final + 
#                           ethnicity + townsend_q5 + education_group + 
#                           diab_fam + bmi_group + smoking_status + 
#                           pa_rec + fruit_veg_rec + satfat_3 + 
#                           alc_group + menopause +  live_births
#                         , data = proteomics)

cox.zph(cox_proteomics2)
# baseline model 
# chisq df     p
# quantile_met               3.87715  4 0.423
# age_first_pregnancy_final  3.75969  1 0.053
# ethnicity                  0.80889  3 0.847
# townsend_q5                2.76700  4 0.598
# education_group            5.59769  4 0.231
# diab_fam                   1.51926  2 0.468
# bmi_group                 15.77784  6 0.015
# smoking_status             0.78267  3 0.854
# pa_rec                     0.00895  2 0.996
# fruit_veg_rec              0.18595  2 0.911
# satfat_3                   0.28582  2 0.867
# alc_group                  1.49397  4 0.828
# menopause                  0.47025  1 0.493
# live_births                1.75437  1 0.185
# GLOBAL                    32.79286 39 0.748
# strata model 
# chisq df    p
# quantile_met               3.406  4 0.49
# age_first_pregnancy_final  2.267  1 0.13
# ethnicity                  1.658  3 0.65
# townsend_q5                2.074  4 0.72
# education_group            4.666  4 0.32
# diab_fam                   0.766  2 0.68
# smoking_status             0.982  3 0.81
# pa_rec                     0.172  2 0.92
# fruit_veg_rec              0.316  2 0.85
# satfat_3                   0.186  2 0.91
# alc_group                  0.738  4 0.95
# menopause                  0.392  1 0.53
# live_births                0.987  1 0.32
# GLOBAL                    17.072 33 0.99
summary(cox_proteomics2)
# coef exp(coef)  se(coef)      z Pr(>|z|)    
# quantile_met2              0.116219  1.123241  0.123568  0.941  0.34695    
# quantile_met3              0.346544  1.414172  0.126870  2.731  0.00630 ** 
# quantile_met4              0.248341  1.281897  0.141306  1.757  0.07884 .  
# quantile_met5              0.492843  1.636963  0.157388  3.131  0.00174 ** 

## Tidy general dataframe quantile----
cox_proteomics2 = broom::tidy(cox_proteomics2, 
                                exponentiate = TRUE, conf.int = TRUE)


## More information -----

info_proteomics = proteomics %>% 
  dplyr::summarise(
    count_total = n(), 
    count_t2dm = sum(t2dm == 1, na.rm = TRUE), 
    perc_t2dm = 100 * mean(t2dm == 1, na.rm = TRUE)
  )

info_proteomics2 = proteomics %>% 
  dplyr::group_by(quantile_met) %>% 
  dplyr::summarise(
    count_total = n(), 
    count_t2dm = sum(t2dm == 1, na.rm = TRUE), 
    perc_t2dm = 100 * mean(t2dm == 1, na.rm = TRUE)
  ) %>% 
  dplyr::select(-c(quantile_met))

# tidy the cox models results -----

cox_metabolomics = cox_metabolomics %>% 
  dplyr::filter(term %in% "z_scaled") 

cox_metabolomics2 = cox_metabolomics2 %>% 
  dplyr::filter(term %in% c("quantile_met2", "quantile_met3", 
                            "quantile_met4", "quantile_met5")) 

cox_proteomics = cox_proteomics %>% 
  dplyr::filter(term %in% "z_scaled") 
  
cox_proteomics2 = cox_proteomics2 %>% 
  dplyr::filter(term %in% c("quantile_met2", "quantile_met3", 
                            "quantile_met4", "quantile_met5")) 

refs = data.frame("term" = "quantile_met1", 
                  "estimate" = 1, 
                  "std.error" = 0, 
                  "statistic" = 0,
                  "p.value" = "", 
                  "conf.low" = 0, 
                  "conf.high" = 0
                  )

cox_results = rbind(cox_metabolomics, refs, cox_metabolomics2, 
                    cox_proteomics, refs, cox_proteomics2)

info = rbind(info_metabolomics, info_metabolomics2, 
             info_proteomics, info_proteomics2)

cox_results = cbind(cox_results, info)


# saveRDS(cox_results, "./upload/2026_01_12_cox_models.RDS")


# Upload the files to RAP UK Bibank----
## Locate all the files generated in the upload folder 

# setwd("./upload")

# system("dx ls data/processed")

# system("pwd")

system("dx upload /home/rstudio-server/upload/*")

