
# Load packages ------
if (!require("pacman")) install.packages("pacman")

pacman::p_load(readr, dplyr, janitor, stringr, tidyr, tibble, Hmisc, limma, 
               missForest, missRanger, doParallel)

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



# Preprocessing -----
## Change the format of column names ----
metabolomics = metabolomics %>%
  dplyr:: rename_with(~ paste0("n_", .x) %>%
                        gsub("-", "_", .) %>%
                        gsub("\\.", "_", .))

names = names %>%
  dplyr::mutate(field_id = paste0("n_", field_id, "_0_0"),
                title = janitor::make_clean_names(title)
  )

metabolites = names$field_id[c(1176:1177, 2313:2561)]
metabolites_label = names$title[c(1176:1177, 2313:2561)]
metabolites_label = c(metabolites_label)

## Select the id and metabolites and rename -----
metabolomics = metabolomics %>%
  dplyr::select("n_eid", all_of(metabolites)) %>%
  dplyr::rename_with(~ metabolites_label, all_of(metabolites)) 

## Evaluate the missingness in metabolomics data and filter -----
missing_values_case = naniar::miss_case_summary(metabolomics)
missing_values_vars = naniar::miss_var_summary(metabolomics)

missing_values_case$pct_miss %>% table()

no_nas = missing_values_case %>%
  dplyr::filter(pct_miss < 20)

nrow(metabolomics)
# 501936
metabolomics = metabolomics %>% dplyr::slice(no_nas$case)
nrow(metabolomics)
# 488083

missing_values_case = naniar::miss_case_summary(metabolomics)
missing_values_vars = naniar::miss_var_summary(metabolomics)

rm(no_nas, missing_values_case, missing_values_vars,
   names, metabolites)

## Select only women with gdm outcome -----

metabolomics = metabolomics %>% 
  dplyr::filter(n_eid %in% gdm$n_eid)

## Imputation of missing values using missForest -----

ids = metabolomics$n_eid

metabolomics_preimputation = metabolomics %>% 
  dplyr::select(glucose_lactate:triglycerides_to_total_lipids_in_small_hdl_percentage) 

metabolomics_preimputation = metabolomics_preimputation %>% 
  dplyr::mutate(across(everything(), as.numeric)) %>% 
  as.data.frame()

t1 = Sys.time()
missranger = missRanger(
  metabolomics_preimputation,
  num.trees = 50,
  mtry = floor(sqrt(ncol(metabolomics_preimputation))),
  maxiter = 10,
  pmm.k = 5, 
  seed = 123456789, 
  num.threads = 40
)
t2 = Sys.time()
t2-t1
# Time difference of 1.376609 hours

metabolomics_imputated = as.data.frame(cbind(n_eid = ids, 
                               missranger))

# system("pkill -9 -u alvaro.torres R")

# sapply(metabolomics_preimputation, class)

# We discarded missForest because is slower than missRanger -----
# n_cores <- parallel::detectCores() - 10
# cl = makeCluster(n_cores)
# registerDoParallel(cl)
# 
# t1 = Sys.time()
# set.seed(123456789)
# forest = missForest(metabolomics_preimputation, 
#                     verbose = TRUE,
#                     ntree = 50, 
#                     maxiter = 2,
#                     parallelize = c('forests'))
# t2 = Sys.time()
# t2-t1
# stopCluster(cl)
# system("pkill -9 -u alvaro.torres R")




## Exporting ------
saveRDS(metabolomics_imputated, "./upload/2025_12_12_metabolomics_imputated.RDS")

# Upload the files to RAP UK Bibank----
## Locate all the files generated in the upload folder 

# setwd("./upload")

# system("dx ls data/processed")

# system("pwd")

system("dx upload /home/rstudio-server/upload/*")

