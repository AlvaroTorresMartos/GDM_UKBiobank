
# Key notes -----
# z=β0 +β1 ×x1 +β2 ×x2 (odd-log)
# e.g. z=−4+2+4=2

# p=1/1 + e^z
# e.g. e^2 = 0.13
# e.g. p = 0.88 (probability of class = 1)
# z=ln(1−p/p)

# Load packages ------
if (!require("pacman")) install.packages("pacman")

pacman::p_load(readr, dplyr,  fastDummies, tibble, glmnet, caret)

# Metabolomics -------
data = read_rds("/mnt/project/data/processed/2025_12_15_metabolomic_dataset_limma_elasticnet.RDS")

vars = data %>% dplyr::select(ethnicity) %>% colnames()

data = fastDummies::dummy_columns(data, select_columns = vars,
                                  remove_first_dummy = TRUE)

limma_results = read_csv("/mnt/project/results/1_limma_analysis/2025_12_15_metabolites_associated_GDM.csv") %>% 
  dplyr::filter(adj.P.Val < 0.05) %>% 
  dplyr::select(Metabolite) %>%  
  dplyr::pull()

data = data %>% 
  tibble::column_to_rownames(var = "n_eid") %>% 
  dplyr::select(age_first_pregnancy_final,
                              ethnicity_2:ethnicity_9,
                              all_of(limma_results), 
                              gdm_final) %>% 
  dplyr::mutate(gdm_final = case_when(gdm_final == "GDM" ~ 1, 
                                      gdm_final == "NoGDM" ~ 0), 
                gdm_final = as.factor(gdm_final), 
                age_first_pregnancy_final = scale(age_first_pregnancy_final)
  )

rm(vars, limma_results)

elasticnet_under = read_rds("/mnt/project/results/2_elasticnet/final_elasticnet_metabolomics.RDS")

metabolic_signature = predict(elasticnet_under, data[, 1:187], type = "prob")

metabolic_signature = metabolic_signature %>% 
  dplyr::select(`1`) %>% 
  dplyr::rename(p = `1`) 

# z=ln(1−p/p)

metabolic_signature = metabolic_signature %>% 
  dplyr::mutate(z = log( ((1-p)/p)), 
                z_scaled = scale(z), 
                z_scaled = as.numeric(z_scaled), 
                z_scaled = z_scaled * -1, 
                eid = row.names(data)) %>% 
  dplyr::relocate(eid, .before = "p")
                
metabolic_signature$z_scaled %>%  hist()

# saveRDS(metabolic_signature, "./upload/metabolic_signature_GDM.RDS")

# Proteomics -----

data = read_rds("/mnt/project/data/processed/2025_12_17_proteomic_dataset_limmma_elasticnet.RDS")

vars = data %>% dplyr::select(ethnicity) %>% colnames()

data = fastDummies::dummy_columns(data, select_columns = vars, 
                                  remove_first_dummy = TRUE)

limma_results = read_csv("/mnt/project/results/1_limma_analysis/2025_12_17_proteins_associated_GDM.csv") %>% 
  dplyr::filter(adj.P.Val < 0.05) %>% 
  dplyr::select(Protein) %>%  
  dplyr::pull()

data = data %>% 
  tibble::column_to_rownames(var = "n_eid") %>% 
  dplyr::select(age_first_pregnancy_final,
                              ethnicity_2:ethnicity_9,
                              all_of(limma_results), 
                              gdm_final) %>% 
  dplyr::mutate(gdm_final = case_when(gdm_final == "GDM" ~ 1, 
                                      gdm_final == "NoGDM" ~ 0), 
                gdm_final = as.factor(gdm_final), 
                age_first_pregnancy_final = scale(age_first_pregnancy_final)
  )

rm(vars, limma_results)

elasticnet_under = read_rds("/mnt/project/results/2_elasticnet/final_elasticnet_proteomics.RDS")

proteomic_signature = predict(elasticnet_under, data[, 1:40], type = "prob")

proteomic_signature = proteomic_signature %>% 
  dplyr::select(`1`) %>% 
  dplyr::rename(p = `1`) 

# z=ln(1−p/p)

proteomic_signature = proteomic_signature %>% 
  dplyr::mutate(z = log( ((1-p)/p)), 
                # z = case_when(is.infinite(z) ~ -35, 
                              # TRUE ~ z), 
                z_scaled = scale(z), 
                z_scaled = as.numeric(z_scaled), 
                z_scaled = z_scaled * -1, 
                eid = row.names(data)) %>% 
  dplyr::relocate(eid, .before = "p")

proteomic_signature$z_scaled %>%  hist()


# saveRDS(proteomic_signature, "./upload/proteomic_signature.RDS")


# Upload the files to RAP UK Bibank----
## Locate all the files generated in the upload folder 

# setwd("./upload")

# system("dx ls data/processed")

# system("pwd")

system("dx upload /home/rstudio-server/upload/*")