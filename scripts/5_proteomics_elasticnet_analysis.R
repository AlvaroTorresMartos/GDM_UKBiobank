

# Load packages ------
if (!require("pacman")) install.packages("pacman")

pacman::p_load(readr, dplyr, janitor, stringr, naniar, imputeLCMD, broom, 
               tidyr, tibble, Hmisc, glmnet, descr, fastDummies, ggplot2, 
               caret, MLmetrics, foreign, parallel, doParallel)

# # Import data -----
# t1 = Sys.time()
# proteomics = read.csv("/mnt/project/data/raw/proteomics.csv")
# 
# gdm = read_rds("/mnt/project/data/processed/2025_04_22_curated_dataset.RDS")
# 
# # check = read_rds("/mnt/project/data/processed/2025_demographics_with_hesoutcome.RDS")
# 
# covariates = read_rds("/mnt/project/data/processed/2025_04_13_covariates.RDS")
# 
# hes = read_rds("/mnt/project/data/processed/HES.RDS")
# 
# t2 = Sys.time()
# t2-t1
# # Time difference of 1.672683 mins
# rm(t1, t2)
# 
# 
# # Preprocessing -----
# ## Change the format of column names ----
# 
# proteomics = proteomics %>% 
#   dplyr::rename_with(.col = 2:last_col(), .fn = toupper)
# 
# 
# ## Evaluate the missingness in proteomics data -----
# 
# missing_values_vars = naniar::miss_var_summary(proteomics, order = FALSE)
# vars_with_acceptable_nas = missing_values_vars %>% 
#   dplyr::filter(pct_miss < 50)
# 
# proteomics = proteomics %>% 
#   dplyr::select(all_of(vars_with_acceptable_nas$variable))
# 
# missing_values_vars = naniar::miss_var_summary(proteomics)
# missing_values_case = naniar::miss_case_summary(proteomics)
# 
# case_with_acceptable_nas = missing_values_case %>% 
#   dplyr::filter(pct_miss < 30)
# 
# proteomics = proteomics %>% 
#   dplyr::slice(case_with_acceptable_nas$case)
# 
# missing_values_vars = naniar::miss_var_summary(proteomics)
# missing_values_case = naniar::miss_case_summary(proteomics)
# 
# rm(vars_with_acceptable_nas, case_with_acceptable_nas, 
#    missing_values_case, missing_values_vars)
# 
# ## Imputation -----
# qrilc_object = imputeLCMD::impute.QRILC(dataSet.mvs = proteomics)
# proteomics_imp = qrilc_object[[1]] 
# 
# proteomics$CTSS %>% hist()
# proteomics_imp$CTSS %>% hist()
# 
# proteomics$CST1 %>% hist()
# proteomics_imp$CST1 %>% hist()
# 
# proteomics$TACSTD2 %>% hist()
# proteomics_imp$TACSTD2 %>% hist()
# 
# rm(qrilc_object, proteomics)
# 
# ## Joining the datasets -----
# 
# gdm = gdm %>% dplyr::select(-c(ethnicity))
# 
# data = dplyr::left_join(gdm, proteomics_imp, by = "eid") %>% 
#   dplyr::mutate(gdm_final = replace_na(gdm_final, 0))
# 
# missing_values_case = naniar::miss_case_summary(data)
# missing_values_vars = naniar::miss_var_summary(data)
# 
# nomissingvalues = missing_values_case %>% 
#   dplyr::filter(pct_miss < 99)
# 
# data = data %>% 
#   dplyr::slice(nomissingvalues$case)
# 
# missing_values_case = naniar::miss_case_summary(data)
# 
# rm(missing_values_case, missing_values_vars, nomissingvalues, gdm, 
#    proteomics_imp)
# 
# 
# covariates = covariates %>% 
#   dplyr::select(eid, ethnicity, townsend_q5, education_group, diab_fam, 
#                 bmi_group, smoking_status, pa_rec, fruit_veg_rec, satfat_3, 
#                 alc_group, menopause)
# 
# data = covariates %>% dplyr::right_join(data, by = "eid")
# 
# 
# 
# ## Filter only the women with GDM (N = 405)
# # data = data %>% dplyr::filter(gdm_final == 1)
# 
# ## Join the T2DM HES 
# hes = hes %>%  dplyr::select(eid, t2dm)
# 
# data = data %>% dplyr::left_join(hes, by = "eid")
# 
# data$t2dm %>% freq()
# # Diagnosis of Type 2 Diabetes Mellitus 
# # Frequency Percent
# # 0         17501  96.525
# # 1           630   3.475
# # Total     18131 100.000
# 
# data$ethnicity %>% freq()
# data$townsend_q5 %>% freq()
# data$education_group %>% freq()
# data$diab_fam %>% freq()
# data$bmi_group %>% freq()
# data$smoking_status %>% freq()
# data$pa_rec %>% freq()
# data$fruit_veg_rec %>% freq()
# data$satfat_3 %>% freq()
# data$alc_group %>% freq()
# data$menopause %>% freq()
# data$t2dm %>% freq()
# 
# data = data %>%
#   dplyr::mutate(
#     ethnicity = as.factor(ethnicity),
#     # country = as.factor(country),
#     townsend_q5 = as.factor(townsend_q5),
#     education_group = as.factor(education_group),
#     diab_fam = case_when(is.na(diab_fam) ~ 9, 
#                          TRUE ~ diab_fam), 
#     diab_fam = as.factor(diab_fam),
#     bmi_group = as.factor(bmi_group),
#     smoking_status = as.factor(smoking_status),
#     pa_rec = as.factor(pa_rec),
#     fruit_veg_rec = as.factor(fruit_veg_rec),
#     satfat_3 = as.factor(satfat_3),
#     alc_group = as.factor(alc_group),
#     menopause = as.factor(menopause),
#     t2dm = case_when(is.na(t2dm)  ~ 2, # si no funciona bien eliminar este linea
#                      TRUE ~ t2dm),
#     t2dm = as.factor(t2dm),
#     gdm_final = as.factor(gdm_final)) %>%
#   dplyr::filter(satfat_3 != 9) %>% 
#   dplyr::filter(menopause != 9) 
# # dplyr::select(-c(gdm_final)) %>%
# # dplyr::filter(ethnicity != 9) %>%
# # dplyr::filter(townsend_q5 != 9) %>%
# # dplyr::filter(education_group != 9) %>%
# # dplyr::filter(bmi_group != 0) %>%
# # dplyr::filter(smoking_status != 9) %>%
# # tidyr::drop_na(diab_fam)
# 
# sum(is.na(data))
# # [1] 0
# 
# 
# data$gdm_final %>% freq()
# #         Frequency  Percent
# # 0         17829  99.6034
# # 1            71   0.3966
# # Total     17900 100.0000
# 
# 
# # saveRDS(data, "/mnt/project/data/processed/2025_04_22_proteomics_dataset_elasticnet.RDS")

## Dummies ------

data = read_rds("/mnt/project/data/processed/2025_12_17_proteomic_dataset_limmma_elasticnet.RDS")

vars = data %>% dplyr::select(ethnicity) %>% colnames()

data = fastDummies::dummy_columns(data, select_columns = vars, 
                                  remove_first_dummy = TRUE)

limma_results = read_csv("/mnt/project/results/1_limma_analysis/2025_12_17_proteins_associated_GDM.csv") %>% 
  dplyr::filter(adj.P.Val < 0.05) %>% 
  dplyr::select(Protein) %>%  
  dplyr::pull()

data = data %>% dplyr::select(age_first_pregnancy_final,
                              ethnicity_2:ethnicity_9,
                              all_of(limma_results), 
                              gdm_final) %>% 
  dplyr::mutate(gdm_final = case_when(gdm_final == "GDM" ~ 1, 
                                      gdm_final == "NoGDM" ~ 0), 
                gdm_final = as.factor(gdm_final), 
                age_first_pregnancy_final = scale(age_first_pregnancy_final)
  )

rm(vars, limma_results)

# Training Predictive modelling (elasticnet) -----

## function for classification metrics including GMean -----
classification_metrics <- function(data, lev = NULL, model = NULL) {
  # Confusion matrix
  cm = confusionMatrix(data = data$pred, reference = data$obs, positive = lev[2])
  
  ##  Accuracy 
  acc_val = cm[["overall"]][["Accuracy"]]
  ## Kappa 
  kap_val = cm[["overall"]][["Kappa"]]
  ## F1
  f1_val = cm[["byClass"]][["F1"]]
  ##  Sensitivity, Specificity
  sensitivity_val = cm$byClass["Sensitivity"] %>% as.numeric()
  specificity_val = cm$byClass["Specificity"] %>% as.numeric()
  ##  Positive Predictive Value and Negative Predictive Value
  pos_pred_val = cm$byClass["Pos Pred Value"] %>% as.numeric()
  neg_pred_val = cm$byClass["Neg Pred Value"] %>% as.numeric()
  ## G-Mean (Geometric Mean)
  gmean_val = sqrt(sensitivity_val * specificity_val)
  
  # Return all metrics in a named vector
  output <- c(
    "GMean" = gmean_val,                          # Geometric Mean
    "Accuracy" = acc_val,                         # Accuracy
    "Kappa" = kap_val,                            # Kappa 
    "F1" = f1_val,                                # F1-Score
    "Sensitivity" = sensitivity_val,              # Sensitivity 
    "Specificity" = specificity_val,              # Specificity
    "Positive_Predictive_Value" = pos_pred_val,   # Positive Predictive Value 
    "Negative_Predictive_Value" = neg_pred_val   # Negative Predictive Value
  )
  return(output)
}

## Set-ups -----
set.seed(123456789)
k1 = caret::createFolds(data$gdm_final, k = 5,
                        list = TRUE, returnTrain =  TRUE)
k2 = caret::createFolds(data$gdm_final, k = 5,
                        list = TRUE, returnTrain =  TRUE)
k3 = caret::createFolds(data$gdm_final, k = 5,
                        list = TRUE, returnTrain =  TRUE)
k4 = caret::createFolds(data$gdm_final, k = 5,
                        list = TRUE, returnTrain =  TRUE)
k5 = caret::createFolds(data$gdm_final, k = 5,
                        list = TRUE, returnTrain =  TRUE)
folds = c(k1, k2, k3, k4, k5)

## Experimental design ------
set.seed(123456789)
control = trainControl(method = "repeatedcv", 
                       number = 5,
                       repeats = 5,
                       savePredictions = "all",
                       summaryFunction = classification_metrics,
                       sampling = "down",
                       index = folds,
                       allowParallel = TRUE
)

### Parallel process and run the elasticnet model 
parallel::detectCores()
# [1] 80
library(doParallel)
getDoParWorkers()
library(doParallel)
cl <- makePSOCKcluster(65)
registerDoParallel(cl)
# https://stackoverflow.com/questions/41117127/r-parallel-programming-error-in-task-1-failed-could-not-find-function
clusterCall(cl, function() library(magrittr))
t1 = Sys.time()
set.seed(123456789)
elasticnet_under = train(x = data[, c(1:40)], y = data$gdm_final, 
                         method = "glmnet",   
                         family = "binomial",
                         penalty.factor = c(rep(0, 5), rep(1, 35)),
                         metric = "GMean",
                         maximize =  TRUE,
                         trControl = control, 
                         tuneGrid = expand.grid(alpha = 0.5, #seq(0, 1, by = 0.1), 
                                                lambda = seq(0, 1, by = 0.01))
                         
)
registerDoSEQ()
unregister <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
stopCluster(cl)
t2 = Sys.time()
t2 - t1
# Time difference of 10.3718 secs
results = elasticnet_under$results

## Save the model as R object ----
# saveRDS(elasticnet_under, "./upload/elasticnet_proteomics.RDS")


# Final Predictive Model ------

elasticnet_under = read_rds("./upload/elasticnet_proteomics.RDS")

set.seed(123456789)
control = trainControl(method = "none",
                       sampling = "down"
)
t1 = Sys.time()
set.seed(123456789)
elasticnet_under = train(x = data[, c(1:40)], y = data$gdm_final, 
                         method = "glmnet",
                         family = "binomial",
                         penalty.factor = c(rep(0, 5), rep(1, 35)),
                         trControl = control,
                         tuneGrid = elasticnet_under$bestTune
)
t2 = Sys.time()
t2 - t1
# Time difference of 0.1271021 secs

## Save the model as R object ----
# saveRDS(elasticnet_under, "./upload/final_elasticnet_proteomics.RDS")

elasticnet_under = read_rds("./upload/final_elasticnet_proteomics.RDS")

pred_under = predict(elasticnet_under, data[, 1:40])
results = broom::tidy(confusionMatrix(data = pred_under, reference = data$gdm_final, 
                                      positive = "1"))
confusionMatrix(data = pred_under, reference = data$gdm_final, positive = "1")

coef_elasticnet = as.matrix(coef(elasticnet_under$finalModel, elasticnet_under$bestTune$lambda))
coef_elasticnet = coef_elasticnet %>% as.data.frame() %>%
  tibble::rownames_to_column(var = "Coefficient")

covariates = data %>% dplyr::select(age_first_pregnancy_final:ethnicity_9) %>% 
  colnames()

coef_elasticnet_under = coef_elasticnet %>% 
  slice(which(coef_elasticnet$Coefficient %in% c(covariates,"(Intercept)")  == FALSE)) %>% 
  filter(s1 != 0)
  
# saveRDS(coef_elasticnet_under, "./upload/coef_proteomic_signature.RDS")

# Coefficient plot -----



coef_elasticnet_under %>% 
  dplyr::rename(value = s1) %>% 
  # dplyr::filter(abs(value) > 0.5) %>% 
  dplyr::mutate(direction = ifelse(value > 0, "+", "-")) %>% 
  ggplot(aes(x = value, y = reorder(Coefficient, value), fill = direction)) + 
  geom_col() + 
  labs(x = "Regularized regression coefficients", y = "") + 
  scale_fill_manual(values = c("#A7C7E7", "#F7A1A1")) +
  # scale_fill_manual(values = c("#0072B2", "#D55E00")) +
  theme_bw() + 
  theme(legend.position = "none", 
        plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 15),
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 16)) #+ 
  # labs(title = "Proteomic signature associated with GDM")#, 
#subtitle = "GDM (Y) ~ Metabolites (X) + counfounders")#, 
#caption = "Metrics in whole dataset (repeated k-fold cross validation, k = 5 with 5 repeats): accuracy 0.64 (0.64), sensitivity: 0.62 (0.58) and specificity: 0.64 (0.62).")

# Upload the files to RAP UK Bibank----
## Locate all the files generated in the upload folder 

# setwd("./upload")

# system("dx ls data/processed")

# system("pwd")

system("dx upload /home/rstudio-server/upload/*")


