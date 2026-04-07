


# Load packages ------
if (!require("pacman")) install.packages("pacman")

pacman::p_load(readr, dplyr,  forcats, ggplot2, waffle)

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


metabolomics = metabolomics %>% 
  dplyr::mutate(quantile_met = ntile(z_scaled, n = 5), 
                quantile_met = factor(quantile_met)
  )

metabolomics2 = metabolomics %>% 
  # merge(gdm, by = "eid") %>% 
  dplyr::mutate(outcome = case_when(gdm_final == 0 & t2dm == 0 ~ "No-GDM-No-T2DM", 
                                    gdm_final == 0 & t2dm == 1 ~ "No-GDM-T2DM", 
                                    gdm_final == 1 & t2dm == 0 ~ "GDM-No-T2DM", 
                                    gdm_final == 1 & t2dm == 1 ~ "GDM-T2DM"), 
                outcome = as.factor(outcome), 
                outcome = forcats::fct_relevel(outcome, "No-GDM-No-T2DM", "No-GDM-T2DM", "GDM-No-T2DM", "GDM-T2DM")
  ) %>% 
  droplevels()


data = metabolomics2 %>% 
  count(outcome, quantile_met) %>%
  group_by(outcome) %>% 
  mutate(
    total_outcome = sum(n),
    perc = round(100 * n / total_outcome)
  ) %>%
  ungroup()



data %>%
  dplyr::mutate(
    perc = case_when(
      outcome == "GDM-T2DM" & quantile_met == 5 ~ 52,
      TRUE ~ perc
    ), 
    Quintile = fct_recode(quantile_met, Q1 = "1", Q2 = "2", 
                          Q3 = "3", Q4 = "4", Q5 = "5"
                          )
  ) %>%
  ggplot(aes(fill = Quintile, values = perc)) +
  waffle::geom_waffle(color = "white", size = 1.5, n_rows = 4) +
  facet_wrap(~outcome, ncol = 1) +
  scale_fill_manual(
    values = c(
      "Q1" = "#fee5d9",  # rosita muy claro
      "Q2" = "#fcae91",  # rojo claro
      "Q3" = "#fb6a4a",  # rojo medio
      "Q4" = "#de2d26",  # rojo fuerte
      "Q5" = "#67000d"   # rojo muy oscuro (casi vino)
    )) +
  coord_equal() +
  theme_void() +
  theme(
    strip.text = element_text(size = 14), 
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11), 
    panel.spacing = unit(1, "lines"), 
    legend.position = "bottom"
             )

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

proteomics2 = proteomics %>% 
  dplyr::mutate(quantile_prot = ntile(z_scaled, n = 5), 
                quantile_prot = factor(quantile_prot)
  )

proteomics2 = proteomics2 %>% 
  dplyr::mutate(outcome = case_when(gdm_final == 0 & t2dm == 0 ~ "No-GDM-No-T2DM", 
                                    gdm_final == 0 & t2dm == 1 ~ "No-GDM-T2DM", 
                                    gdm_final == 1 & t2dm == 0 ~ "GDM-No-T2DM", 
                                    gdm_final == 1 & t2dm == 1 ~ "GDM-T2DM"), 
                outcome = as.factor(outcome), 
                outcome = forcats::fct_relevel(outcome, "No-GDM-No-T2DM", "No-GDM-T2DM", "GDM-No-T2DM", "GDM-T2DM")
  ) %>% 
  droplevels()


data = proteomics2 %>% 
  count(outcome, quantile_prot) %>%
  group_by(outcome) %>% 
  mutate(
    total_outcome = sum(n),
    perc = round(100 * n / total_outcome)
  ) %>%
  ungroup()



data %>%
  dplyr::mutate(
    perc = case_when(outcome == "No-GDM-No-T2DM" & quantile_prot == 5 ~ 20,
    TRUE ~ perc), 
    Quintile = fct_recode(quantile_prot, Q1 = "1", Q2 = "2", 
                          Q3 = "3", Q4 = "4", Q5 = "5"
    )
  ) %>%
  ggplot(aes(fill = Quintile, values = perc)) +
  waffle::geom_waffle(color = "white", size = 1.5, n_rows = 4) +
  facet_wrap(~outcome, ncol = 1) +
  scale_fill_manual(
    values = c(
      "Q1" = "#fee5d9",  # rosita muy claro
      "Q2" = "#fcae91",  # rojo claro
      "Q3" = "#fb6a4a",  # rojo medio
      "Q4" = "#de2d26",  # rojo fuerte
      "Q5" = "#67000d"   # rojo muy oscuro (casi vino)
    )) +
  coord_equal() +
  theme_void() +
  theme(
    strip.text = element_text(size = 14), 
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11), 
    panel.spacing = unit(1, "lines"), 
    legend.position = "bottom"
  )


