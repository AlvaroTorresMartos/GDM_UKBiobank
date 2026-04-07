
# Load packages ------
if (!require("pacman")) install.packages("pacman")

pacman::p_load(readr, dplyr, janitor, stringr, tidyr, tibble, Hmisc, limma)

# Import data -----
t1 = Sys.time()
proteomics = read_rds("/mnt/project/data/processed/2025_12_15_proteomics_imputated.RDS")

gdm = read_rds("/mnt/project/data/processed/2025_12_10_curated_dataset.RDS")

# check = read_rds("/mnt/project/data/processed/2025_demographics_with_hesoutcome.RDS")

covariates = read_rds("/mnt/project/data/processed/2025_12_10_covariates.RDS")

t2 = Sys.time()
t2-t1
# Time difference of 1.795043 mins
rm(t1, t2)

## Joining the datasets -----

gdm = gdm %>% dplyr::select(-c(ethnicity))

data = dplyr::right_join(gdm, proteomics, by = "n_eid") %>% 
  dplyr::mutate(gdm_final = replace_na(gdm_final, 0))

covariates = covariates %>% dplyr::select(eid, ethnicity) %>% 
  dplyr::rename(n_eid = eid)

data = data %>% dplyr::left_join(covariates, by = "n_eid")

missing_values_case = naniar::miss_case_summary(data)
missing_values_vars = naniar::miss_var_summary(data)

# nomissingvalues = missing_values_case %>% 
#   dplyr::filter(pct_miss < 99)
# 
# data = data %>% 
#   dplyr::slice(nomissingvalues$case)
# 
# missing_values_case = naniar::miss_case_summary(data)

# 8 participant who decide to remove their data in the UKBB dataset
data = data %>% na.omit()

rm(missing_values_case, missing_values_vars, gdm, proteomics, covariates)

nrow(data)
# 17895


data = data %>% 
  dplyr::mutate(
    ethnicity = as.factor(ethnicity), 
    gdm_final = if_else(gdm_final == 1, "GDM", "NoGDM"),
    gdm_final = as.factor(gdm_final)
  ) %>% 
  dplyr::relocate(gdm_final, .after = last_col()) %>% 
  dplyr::relocate(ethnicity, .after = age_first_pregnancy_final)

data$gdm_final %>% table(useNA = "ifany")
# GDM NoGDM 
# 62 17833 

## Scaling -----

proteins_label = data %>% 
  dplyr::select(A1BG:ZPR1) %>% 
  colnames()

data = data %>% 
  dplyr::mutate(across(all_of(proteins_label), scale)) %>% 
  dplyr::mutate(id = paste0("participant_", 1:nrow(data))) 

saveRDS(data, "./upload/2025_12_17_proteomic_dataset_limmma_elasticnet.RDS")

## Proteins (expression matrix) -----

data = readr::read_rds("/mnt/project/data/processed/2025_12_17_proteomic_dataset_limmma_elasticnet.RDS")

proteins = colnames(data)[4:(ncol(data) - 2)]

proteomics = data %>% 
  # dplyr::mutate(id = paste0("participant_", 1:nrow(data))) %>% 
  dplyr::select(id, all_of(proteins)) %>% 
  tidyr::pivot_longer(-id, names_to = "proteins", values_to = "Value") %>% 
  tidyr::pivot_wider(names_from = id, values_from = Value) %>% 
  tibble::column_to_rownames(var = "proteins")

colnames(proteomics)[1:5]
# [1] "participant_1" "participant_2" "participant_3" "participant_4" "participant_5"
row.names(proteomics)[1:5]
# "A1BG"    "AAMDC"   "AARSD1"  "ABCA2"   "ABHD14B"

## Covariates and outcome (design matrix) -----

covariate_outcome = data %>% 
  dplyr::mutate(id = paste0("participant_", 1:nrow(data))) %>% 
  dplyr::select(id, age_first_pregnancy_final, ethnicity, gdm_final) %>% 
  dplyr::mutate(across(c(ethnicity, gdm_final), as.factor)) %>% 
  remove_rownames() %>% 
  tibble::column_to_rownames(var = "id")

head(covariate_outcome)
#                  age_first_pregnancy_final ethnicity gdm_final
# participant_1                        25         1     NoGDM
# participant_2                        26         4     NoGDM
# participant_3                        26         1     NoGDM
# participant_4                        38         1     NoGDM
# participant_5                        41         2     NoGDM
# participant_6                        28         1     NoGDM


# Limma analysis -----

## Design matrix -----
design_matrix = model.matrix( ~ 0 + gdm_final + age_first_pregnancy_final + ethnicity, 
                              data = covariate_outcome)

colnames(design_matrix)

colnames(design_matrix) = c("yesGDM", "noGDM", "maternal_age", 
                            "ethnicity2", "ethnicity3", "ethnicity4", "ethnicity9")

## Contrast matrix -----

contrast_matrix = makeContrasts(
  yesGDM_vs_noGDM = yesGDM-noGDM, 
  levels=design_matrix)

## Fitting model -----

fit = lmFit(object = proteomics, design = design_matrix)

## Compute Contrasts from Linear Model Fit -----

contrasts_lm_fit = contrasts.fit(fit, contrast_matrix)

## Empirical Bayes Statistics for Differential Expression (Levels or Concentration in proteomics) -----

ebayes_lm_fit = eBayes(contrasts_lm_fit)

## Multiple Testing Across Genes (Metabolites) and Contrasts -----

test = decideTests(ebayes_lm_fit, method="separate", 
                   adjust.method = "fdr", p.value=0.05)

summary(test)
#             yesGDM_vs_noGDM
# Down                26
# NotSig            2896
# Up                  9


## Table of Top Metabolites from Linear Model Fit
results = topTable(ebayes_lm_fit, number=nrow(ebayes_lm_fit),  
                   coef="yesGDM_vs_noGDM",
                   adjust="fdr") 

results = results %>% 
  tibble::rownames_to_column(var = "Protein") %>% 
  dplyr::mutate(FC = ifelse(logFC >= 0, 2^logFC, -1 / (2^logFC))) %>% 
  dplyr::relocate(FC, .before = logFC)

# Export results  ------ 
write.csv(results, "/mnt/project/results/1_limma_analysis/2025_12_17_proteins_associated_GDM.csv",
          row.names = FALSE)

# Export results  ------ 
# write.csv(results, "./upload2025_12_15_metabolites_associated_GDM.csv",
#           row.names = FALSE)

# Upload the files to RAP UK Bibank----
## Locate all the files generated in the upload folder 

# setwd("./upload")

# system("dx ls data/processed")

# system("pwd")

system("dx upload /home/rstudio-server/upload/*")


