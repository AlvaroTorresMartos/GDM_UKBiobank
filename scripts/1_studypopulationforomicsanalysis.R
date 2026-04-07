
# Load packages ------
if (!require("pacman")) install.packages("pacman")

pacman::p_load(readr, dplyr, stringr, tidyr, lubridate)

# Import data -----
t1 = Sys.time()
demo = readr::read_csv("/mnt/project/data/raw/Demographics.csv")

hes = readr::read_csv("/mnt/project/data/raw/Health_Outcomes_v2.csv")

bio = readr::read_csv("/mnt/project/data/raw/old/Biomarkers.csv")

t2 = Sys.time()
t2-t1
# Time difference of 1.381225 mins

# Preprocessing -----
## Change the format of column names ----

demo = demo %>% 
 dplyr:: rename_with(~ paste0("n_", .x) %>%
                gsub("-", "_", .) %>%
                gsub("\\.", "_", .))



hes = hes %>% 
  dplyr:: rename_with(~ paste0("n_", .x) %>%
                        gsub("-", "_", .) %>%
                        gsub("\\.", "_", .))

bio = bio %>% 
  dplyr:: rename_with(~ paste0("n_", .x) %>%
                        gsub("-", "_", .) %>%
                        gsub("\\.", "_", .))

## Merge the datasets and remove the old ones ------

hes = hes %>% dplyr::select(-c(n_31_0_0, n_21022_0_0))

data = merge(demo, hes, by = "n_eid") %>% 
  merge(bio, by = "n_eid")

# rm(demo, hes)

### Remove the duplicated columns (after merging) ----
data = data %>% 
  dplyr::select(-ends_with(".y")) %>% 
  dplyr::rename_with(~ gsub(".x", "", .))

data$n_31_0_0 %>% table(useNA = "ifany")
# 0      1 
# 273155 228973 

## Filter only female participants  -----
# Sex (female coded by 0)
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=31
data = data %>% 
  dplyr::filter(n_31_0_0 == 0)
nrow(data)
# [1] 273155

## Define women with GDM -----
# GDM_Q
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=4041
# 4041: Gestational diabetes only
# This is variable 4041: GDM only - i.e. people who previously said they have 
# a diabetes diagnosis from a doctor and then answered whether their diabetes was 
# ONLY during pregnancy. Possible answers were:
# * - Yes (1) - we assume this is true GDM
# * - No (0) - i.e. not only during pregnancy but also at other points in their life. 
# This means either they already had diabetes before pregnancy, or they had GDM, 
# and then went on to develop diabetes later in life - for this analysis, we will 
# exclude any known pre-pregnancy diabetes, but we will assume all the rest means 
# no GDM
# * - Not applicable (-2) - we assume this means no GDM
# * - Unsure (-1) - code as missing
# * - Prefer not to answer (-3) - code as missing
# GDM_Int
# *This is variable 20002, and the answer 1221 corresponds to GDM.
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20002
# 20002: Non-cancer illness code, self-reported
# HES_GDM
# https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=2002
# https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=19&nl=1
# n_41270_0*
# *GDM is most probably code O244 - diabetes arising in pregnancy. 
# This doesn't necessarily mean true GDM, but the codes below are pre-existing 
# or unspecified diabetes, so I assumed O24.4 could be true GDM, plus O249 (unspecified diabetes in pregnancy)
# *O240 is pre-existing type 1 diabetes
# *O241 is pre-existing type 2 diabetes
# *O242 is pre-existing malnutrition-related diabetes
# *O243 is pre-existing unspecified diabetes
data2 = data %>% 
  dplyr::mutate(gdm_q = if_else(n_4041_0_0 == 1, 1, NA_real_), 
                gdm_int = if_else(rowSums(select(., starts_with("n_20002_0_")) == 1221, na.rm = TRUE) > 0, 1, NA_real_),
                hes_gdm = if_else(
                  rowSums(across(c(n_41270_0_0, n_41271_0_0, n_41202_0_0, 
                                   n_41203_0_0, n_41204_0_0, n_41205_0_0),
                                 ~ str_detect(., "O244|O249"))) > 0,
                  1, NA_real_), 
                gdm_final = dplyr::coalesce(gdm_q, gdm_int, hes_gdm)
                ) 




data2$gdm_q %>% table(useNA = "ifany")  
#     1   <NA> 
#   1070 272085   

data2$gdm_int %>% table(useNA = "ifany")  
#  1   <NA> 
# 232 272923 
  
data2$hes_gdm %>% table(useNA = "ifany") 
# <NA> 
#   273155 

data2$gdm_final %>% table(useNA = "ifany") 
# 1   <NA> 
#   1170 271985   

# hes_code = data2 %>% 
#   dplyr::select(eid, n_41270_0_0) %>% 
#   tidyr::separate(n_41270_0_0, 
#                   into = paste0("code_", 1:259), 
#                   sep = "\\|", 
#                   fill = "right")
# 
# hes_dates = data2 %>% 
#   dplyr::select(eid, starts_with("n_41280_0"))
# 
# dates_gdm = data.frame(eid = data2$eid, 
#                        date_244 = rep(NA_character_, nrow(hes_code)),
#                        date_249 = rep(NA_character_, nrow(hes_code)))
# 
# t1 = Sys.time()
# for (j in 1:nrow(hes_code)) {
#   print(j)
#   for (i in 2:ncol(hes_code)) {
#     if (!is.na(hes_code[j, i])) {
#       if (hes_code[j, i] == "O244") {
#         dates_gdm$date_244[j] = hes_dates[j, i]
#       }
#       if (hes_code[j, i] == "O249") {
#         dates_gdm$date_249[j] = hes_dates[j, i]
#       }
#     }
#   }
# }
# t2 = Sys.time()
# t2 - t1
# Time difference of 14.19295 mins
# saveRDS(dates_gdm, "/mnt/project/data/processed/2025_01_17_dates_gdm.RDS")
dates_gdm = readRDS("/mnt/project/data/processed/2025_01_17_dates_gdm.RDS")

dates_gdm = dates_gdm %>% dplyr::filter(eid %in% data2$n_eid)

# 21022: Age at recruitment
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=21022
data2 = data2 %>% 
  dplyr::mutate(
    dob = lubridate::make_date(year = n_34_0_0, month = n_52_0_0 , 
                               day = 01), 
    date_244 = as.Date(dates_gdm$date_244, format = "%Y-%m-%d"),
    date_249 = as.Date(dates_gdm$date_249, format = "%Y-%m-%d"),
    date_gdm = coalesce(date_244, date_249),
    age_gdm_hes = round(as.numeric((date_gdm - dob))/365.25, digits = 0),
    gdm_hes_before_baseline = ifelse(n_21022_0_0 > age_gdm_hes, TRUE, FALSE), 
  )

data2$gdm_hes_before_baseline %>% table(useNA = "ifany")
# FALSE   TRUE   <NA> 
#   29    226 272903 

### 1) Filtering by GDM after baseline (N = 29) ------
nrow(data2)
# 273155
data2 = data2 %>% 
  dplyr::filter(gdm_hes_before_baseline == TRUE | is.na(gdm_hes_before_baseline))
nrow(data2)
# 273126

## Calculate covariates and refining the number of participants included in our dataset ------ 
# Variables in touchscreen: Female-specific factors
# https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=100069
# Variables in Hospital inpatient: Summary Maternity 
# https://biobank.ndph.ox.ac.uk/showcase/label.cgi?id=2003

# Touchscreen
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=2734
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=2744
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=2754
# 2734: Number of live births
# 2744: Birth weight of first child
# 2754: Age at first live birth
# 2764: Age at last live birth
# 3872: Age of primiparous women at birth of child

# Hospital inpatient (HES)
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=41286 
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=41284
# 41286: Mother age on date of delivery
# 41284: Baby birth weight

data3 = data2 %>% 
  dplyr::mutate(across(starts_with("n_41284_0_"), as.numeric)) %>% 
  # dplyr::mutate(across(starts_with("n_41284_0_"), 
  #                      ~ if_else(. == 9999, NA, .))) %>%  
  dplyr::mutate(live_births = if_else(n_2734_0_0 %in% c(0, -3), NA_real_, n_2734_0_0), 
                child_weight = if_else(n_2744_0_0 %in% c(-1, -2, -3), NA_real_, n_2744_0_0),
                matage_extra = if_else(n_3872_0_0 %in% c(-3, -4), NA_real_, n_3872_0_0),
                matage_first = if_else(n_2754_0_0 %in% c(-3, -4), NA_real_, n_2754_0_0), 
                matage_last = if_else(n_2764_0_0 %in% c(-3, -4), NA_real_, n_2764_0_0), 
                matage = coalesce(!!!select(., starts_with("n_41286_0_"))), 
                baby_weight = coalesce(!!!select(., starts_with("n_41284_0_"))), 
                births_wrong = if_else(is.na(live_births) & is.na(child_weight) &
                                       is.na(matage_first) & is.na(matage_last) &
                                       is.na(matage_extra) & is.na(matage) & 
                                       is.na(baby_weight), 1, 0)
                )


 
data3 %>% 
  dplyr::filter(gdm_final == 1) %>% 
  select(births_wrong) %>% table(useNA = "ifany")
# births_wrong
#   0    1 
# 1251   11

### Generating the final age of the first pregnancy -----
data3 = data3 %>% 
  dplyr::mutate(age_first_pregnancy_final = pmin(matage_extra, matage_first, 
                                                 matage, na.rm = TRUE))


### 2) Filtering by age at first birth and lack of information about pregnancy (birth wrong)  ----- 
## Remove 52038 because: 
# - one women because she was 10 years old when she had the first pregnancy
# - 11 due to lack of information about pregnancy (birth wrong == 1), 
# one of them is the women with 10 years
# - 52026 because we have missing values in the age of the first pregnancy
# (52015 because the 11 birth wrong have not pregnancy information)
data3$age_first_pregnancy_final %>% summary()
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   10.00   22.00   25.00   25.94   29.00   54.00   52026 


nrow(data3)
# 273126

data3 = data3 %>% 
  dplyr::filter(!is.na(age_first_pregnancy_final)) %>% 
  dplyr::filter(age_first_pregnancy_final != 10) 
nrow(data3)
# 221099

data3$gdm_final %>% table(useNA = "ifany")
# 1   <NA> 
#   1142 219957 

## Decision made in the Moscho's paper 
# *Although the people who were included as GDM 1 above should be true GDM, 
# there may be some mistakes. Thus, need to make sure anyone coded as GDM 1, 
# did not have any pre-existing diabetes.
# *Ideally, I need to know the age/date of GDM diagnosis, then the age/date 
# of any other diabetes diagnosis, and if other diabetes preceded pregnancy 
# with GDM, exclude.
# *Explore if I can find age of GDM:
# *Therefore, decided as rule of thumb to compare for all participants, 
# the age for first pregnancy with the age of other diabetes, and exclude
# if age other diabetes<=age first pregnancy, noting the limitation of this.

## Define women with diabetes from all sources -----
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=2443
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=6153
# Questionnaire 
# *Related fields include:
#   2443 is a problematic question because most people with GDM answered Yes 
#   by the GDM (nor ordinary diabetes)
#   * 2443: Has a doctor ever told you that you have diabetes? 
#   (1=Yes, 0=No, -3=Prefer not to answer, -1=Don't know:then interview)
# 	* 6153: Do you regularly take any medication - option for insulin 
# 1	Cholesterol lowering medication
# 2	Blood pressure medication
# 3	Insulin
# 4	Hormone replacement therapy
# 5	Oral contraceptive pill or minipill
# -7	None of the above
# -1	Do not know
# -3	Prefer not to answer
# 	* 2492: Taking other prescription medication - if yes, then interview

data4 = data3 %>% 
  dplyr::mutate(diab_diag_q = case_when(n_2443_0_0 %in% c(-1, -3) ~ NA_real_, 
                                      n_2443_0_0 == 0 ~ NA_real_,
                                      TRUE ~ n_2443_0_0), 
                diag_med_q = if_else(n_6153_0_0 == 3, 1, NA_real_),
                medication_cholesterol = if_else(n_6153_0_0 == 1, 1, NA_real_),  
                diabetes_q = coalesce(diab_diag_q, diag_med_q))

data4$diab_diag_q %>% table(useNA = "ifany")
# 1   <NA> 
#   8624 212475 
data4$diag_med_q %>% table(useNA = "ifany")
# 1   <NA> 
#   216 220883
data4$diabetes_q %>% table(useNA = "ifany")
# 1   <NA> 
#   8643 212456 

table(data4$gdm_final == data4$diab_diag_q, useNA = "ifany")
# TRUE   <NA> 
#   1090 220009 
table(data4$gdm_final == data4$diag_med_q, useNA = "ifany")
# TRUE   <NA> 
#   17 221082  

### 3) Filtering by diabetes medication questionaries (N = 216) ----- 
nrow(data4)
# 221099
data4 = data4 %>% 
  dplyr::filter(is.na(diag_med_q))
nrow(data4)
# 220883

## * Using HbA1c measured at UKBB assessment centre (30750)- diabetes if HbA1c >=48mmol/L - mean of all timepoints

data4 = data4 %>% dplyr::left_join(bio, by = "n_eid")
data4 = data4 %>% 
  dplyr::select(-ends_with(".y")) %>% 
  dplyr::rename_with(~ gsub(".x", "", .))

data4 = data4 %>% 
  dplyr::mutate(diab_hba1c = if_else(n_30750_0_0 >= 48, 1, NA_real_))
  
data4$diab_hba1c %>% table(useNA = "ifany")
#    1   <NA> 
#   5340 215543
table(data4$gdm_final == data4$diab_hba1c, useNA = "ifany")
# TRUE   <NA> 
#   266 220617 

### 4) Filtering by levels of Hba1c (N = 5340) ----- 
nrow(data4)
# 220883
data4 = data4 %>% 
  dplyr::filter(is.na(diab_hba1c))
nrow(data4)
# 215543



## Interview: 
# 20002:Non-cancer illness code, self-reported
# 1223 is T2DM, 1222 is T1DM, and 1220 is generic for diabetes
  
data4 = data4 %>% 
  dplyr::mutate(
     gendiab_int = if_else(rowSums(select(., starts_with("n_20002_0_")) == 1220, 
                                   na.rm = TRUE) > 0, 1, NA_real_),
     t1dm_diab_int = if_else(rowSums(select(., starts_with("n_20002_0_")) == 1222, 
                                     na.rm = TRUE) > 0, 1, NA_real_),
     t2dm_diab_int = if_else(rowSums(select(., starts_with("n_20002_0_")) == 1223, 
                                     na.rm = TRUE) > 0, 1, NA_real_),
     diabetes_int = dplyr::coalesce(gendiab_int, t1dm_diab_int, t2dm_diab_int))


data4$gendiab_int %>% table(useNA = "ifany")
#     1   <NA> 
#   2946 212597
table(data4$gdm_final == data4$gendiab_int, useNA = "ifany")
# TRUE   <NA> 
#   29 215514 
data4$t1dm_diab_int %>% table(useNA = "ifany")
# 1   <NA> 
#   30 215513
table(data4$gdm_final == data4$t1dm_diab_int, useNA = "ifany")
# TRUE   <NA> 
#   9 215534 
data4$t2dm_diab_int %>% table(useNA = "ifany")
# 1   <NA> 
#   505 215038 
table(data4$gdm_final == data4$t2dm_diab_int, useNA = "ifany")
# TRUE   <NA> 
#   78 215465
data4$diabetes_int %>% table(useNA = "ifany")
# 1   <NA> 
#   3474 212069 

### 5) Filtering by diabetes in interviews (N = 3474) ----- 
nrow(data4)
# 215543
data4 = data4 %>% 
  dplyr::filter(is.na(diabetes_int))
nrow(data4)
# 212069

## HES DIABETES (hes_otherdiabetes)

# data4 = data4 %>% 
#   dplyr::mutate(
#     hes_otherdiabetes = if_else(
#       rowSums(across(c(n_41270_0_0, n_41271_0_0, n_41202_0_0, 
#                        n_41203_0_0, n_41204_0_0, n_41205_0_0), 
#                      ~ str_detect(., "O240|O241|O242|O243"))) > 0,
#       1, NA_real_))
# 
# data4$hes_otherdiabetes %>% table(useNA = "ifany")
#   1   <NA> 
#   10 212062 

### 6) Filtering by diabetes in HES (N = 10) ----- 
# nrow(data4)
# # 212072
# data4 = data4 %>% 
#   dplyr::filter(is.na(hes_otherdiabetes))
# nrow(data4)
# 212062

### 7) Filtering by CVD in questionaries (N = 6950) ----- 
# 6150: Vascular/heart problems diagnosed by doctor
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=6150
# 1	Heart attack
# 2	Angina
# 3	Stroke
# 4	High blood pressure
# -7	None of the above
# -3	Prefer not to answer
# data4 %>%  
#   dplyr::filter(gdm_final == 1) %>% 
#   dplyr::select(n_6150_0_0) %>% table(useNA = "ifany")
# -3    -7     1   1|2 1|2|4   1|3   1|4     2   2|4     3   3|4     4 
#  2   582     3     1     1     1     2     4     7     4     1   202 

data4 = data4 %>% 
  dplyr::mutate(
    cvd_exclusion = case_when(n_6150_0_0 %in% c("-7", "-3") ~ NA_character_, 
                              n_6150_0_0 == "4" ~ NA_character_,
                              n_6150_0_0 == "" ~ NA_character_,
                              TRUE ~ n_6150_0_0))

table(data4$cvd_exclusion, useNA = "ifany")
# 1        1|2   1|2|3 1|2|3|4   1|2|4     1|3   1|3|4     1|4     2     2|3   2|3|4     2|4       3     3|4    <NA> 
# 547     228      19      60     350      24      33     373    1532      58     106    1579    1098     943  205112

nrow(data4)
# 212069
data4 = data4 %>% 
  dplyr::filter(is.na(cvd_exclusion))
nrow(data4)
# 205119


## Calculate covarite: ethnicity  -----
# Ethnic background
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=21000

data4 = data4 %>% 
  dplyr::mutate(    
    n_21000_0_0 = if_else(n_21000_0_0 %in% c("-3", "-1"), NA_real_, n_21000_0_0), 
    n_21000_0_0 = case_when(
      n_21000_0_0 %in% c("1", "1001", "1002", "1003") ~ 1,
      n_21000_0_0 %in% c("4", "4001", "4002", "4003") ~ 2, 
      n_21000_0_0 %in% c("3", "5", "3001", "3002", "3003", "3004") ~ 3,
      n_21000_0_0 %in% c("2", "6", "2001", "2002", "2003", "2004") ~ 4, 
      is.na(n_21000_0_0) ~ 4,  
      TRUE ~ n_21000_0_0)) %>% 
  dplyr::rename(ethnicity = n_21000_0_0)

data4 %>% 
  dplyr::filter(gdm_final == 1) %>% 
  dplyr::select(ethnicity) %>% 
  table(useNA = "ifany")

# ethnicity
# 1   2   3   4 
# 642  18  40  19

table(data4$ethnicity, useNA = "ifany")
# 1      2      3      4 
# 194539   3164   3825   3591 

# HES excluding ----
hes = readRDS("/mnt/project/data/processed/2025_demographics_with_hesoutcome.RDS")

hes1 = hes %>% dplyr::filter(t2dm_inc == 1 | is.na(t2dm_inc))

hes2 = hes %>% 
  dplyr::mutate(cvd_inc = case_when(
  chd_inc == 2  ~ 2, 
  chf_inc == 2 ~ 2, 
  stroke_inc == 2 ~ 2,
  TRUE ~ 1)) %>% 
  dplyr::filter(cvd_inc == 1 | is.na(cvd_inc))

# CVD HES (N = 2062)

nrow(data4)
# 205119
data5 = data4 %>% dplyr::filter(n_eid %in% hes1$eid)
nrow(data5)
# 204959

nrow(data5)
# 204959
data5 = data5 %>% dplyr::filter(n_eid %in% hes2$eid)
nrow(data5)
# 203057

# Sample size  -----
table(data5$gdm_final, useNA = "ifany")
#     1   <NA> 
#   706 202351 


# Export study population for omics analysis (N = 205112) ----

data_def = data5 %>%
  dplyr::select(n_eid, age_first_pregnancy_final, ethnicity, gdm_final)

# saveRDS(data_def, "./upload/2025_07_09_curated_dataset.RDS")

# Upload the files to RAP UK Bibank----
## Locate all the files generated in the upload folder 

# setwd("./upload")

# system("dx ls data/processed")

# system("pwd")

system("dx upload /home/rstudio-server/upload/*")

