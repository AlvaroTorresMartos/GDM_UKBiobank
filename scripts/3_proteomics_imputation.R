
# Load packages ------
if (!require("pacman")) install.packages("pacman")

pacman::p_load(readr, dplyr, janitor, stringr, tidyr, tibble, Hmisc, limma, 
               missForest, missRanger, doParallel)

# Import data -----
t1 = Sys.time()
proteomics = read_csv("/mnt/project/data/raw/proteomics.csv")

gdm = read_rds("/mnt/project/data/processed/2025_12_10_curated_dataset.RDS")

t2 = Sys.time()
t2-t1
# Time difference of 1.672683 mins
rm(t1, t2)


# Preprocessing -----
## Change the format of column names ----

proteomics = proteomics %>% 
  dplyr::rename_with(.col = 2:last_col(), .fn = toupper)


## Evaluate the missingness in proteomics data -----

missing_values_vars = naniar::miss_var_summary(proteomics, order = FALSE)
vars_with_acceptable_nas = missing_values_vars %>% 
  dplyr::filter(pct_miss < 50)

proteomics = proteomics %>% 
  dplyr::select(all_of(vars_with_acceptable_nas$variable))

missing_values_vars = naniar::miss_var_summary(proteomics)
missing_values_case = naniar::miss_case_summary(proteomics)

case_with_acceptable_nas = missing_values_case %>% 
  dplyr::filter(pct_miss < 30)

proteomics = proteomics %>% 
  dplyr::slice(case_with_acceptable_nas$case)

missing_values_vars = naniar::miss_var_summary(proteomics)
missing_values_case = naniar::miss_case_summary(proteomics)

rm(vars_with_acceptable_nas, case_with_acceptable_nas, 
   missing_values_case, missing_values_vars)



# Import data -----
t1 = Sys.time()
metabolomics = read_csv("/mnt/project/data/raw/Biomarkers.csv")

gdm = read_rds("/mnt/project/data/processed/2025_12_10_curated_dataset.RDS")

# covariates = read_rds("/mnt/project/data/processed/2025_12_10_covariates.RDS")

names = readr::read_tsv("/mnt/project/data/raw/fieldsum.tsv")

t2 = Sys.time()
t2-t1
# Time difference of 1.622121 mins
rm(t1, t2)

## Select only women with gdm outcome -----

proteomics = proteomics %>% 
  dplyr::filter(eid %in% gdm$n_eid)

## Imputation of missing values using missForest -----

ids = proteomics$eid

proteomics_preimputation = proteomics %>% 
  dplyr::select(2:last_col()) 

proteomics_preimputation = proteomics_preimputation %>% 
  dplyr::mutate(across(everything(), as.numeric)) %>% 
  as.data.frame()

t1 = Sys.time()
missranger = missRanger(
  proteomics_preimputation,
  num.trees = 50,
  mtry = NULL,
  maxiter = 2,
  pmm.k = 2, 
  seed = 123456789, 
  num.threads = 40
)
t2 = Sys.time()
t2-t1

proteomics_imputated = cbind(n_eid = ids, 
                               missranger)

# system("pkill -9 -u alvaro.torres R")

## Exporting ------
saveRDS(proteomics_imputated, "./upload/2025_12_15_proteomics_imputated.RDS")

# Upload the files to RAP UK Bibank----
## Locate all the files generated in the upload folder 

# setwd("./upload")

# system("dx ls data/processed")

# system("pwd")

system("dx upload /home/rstudio-server/upload/*")

