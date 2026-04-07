
# Load packages ------
if (!require("pacman")) install.packages("pacman")

pacman::p_load(readr, dplyr, readr, descr, stringr, tidyr)

# Import data -----
t1 = Sys.time()
demo = readr::read_csv("/mnt/project/data/raw/Demographics.csv")

hes = readr::read_csv("/mnt/project/data/raw/Health_Outcomes_v2.csv")

bio = readr::read_csv("/mnt/project/data/raw/Biomarkers.csv")

diet = readr::read_csv("/mnt/project/data/raw/Diet.csv")

hes = hes %>% dplyr::select(-c("21022-0.0"))

diet = diet %>% dplyr::select(-c("21022-0.0", "31-0.0"))

t2 = Sys.time()
t2-t1
# Time difference of 1.563938 mins

# Merge the datasets and remove the old ones ------
data = merge(demo, hes, by = "eid") %>% 
  merge(bio, by = "eid") %>% 
  merge(diet, by = "eid")

rm(bio, demo, hes, diet)

# Preprocessing -----
## Change the format of column names ----

data = data %>% 
  dplyr:: rename_with(~ gsub("-", "_", .) %>%
                        gsub("\\.", "_", .))

# Covariates ----


## Age groups ----
# 21022 Age at recruitment
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=21022

data2 = data %>% 
  dplyr::mutate(age_5y_group = case_when(
    `21022_0_0` >= 39 &  `21022_0_0` < 45 ~ 1, 
    `21022_0_0` >= 45 &  `21022_0_0` < 50 ~ 2, 
    `21022_0_0` >= 50 &  `21022_0_0` < 55 ~ 3, 
    `21022_0_0` >= 55 &  `21022_0_0` < 60 ~ 4, 
    `21022_0_0` >= 60 &  `21022_0_0` < 65 ~ 5, 
    `21022_0_0` >= 65 &  `21022_0_0` < 70 ~ 6, 
    `21022_0_0` >= 70  ~ 7, 
  ))

data2$age_5y_group %>% freq()
# Frequency   Percent Valid Percent
# 1         51727 1.030e+01       10.3015
# 2         65988 1.314e+01       13.1416
# 3         76260 1.519e+01       15.1872
# 4         90760 1.807e+01       18.0749
# 5        121417 2.418e+01       24.1803
# 6         93558 1.863e+01       18.6322
# 7          2422 4.823e-01        0.4823
# NA's          2 3.983e-04              
# Total    502134 1.000e+02      100.0000


## BMI group  ------
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=21001
# 21001: Body mass index (BMI)

# Calculate the row mean for BMI across the different time points
data2 = data2 %>%
  dplyr::rename(bmi = `21001_0_0`)

# Check the missing values for the new BMI variable
# sum(is.na(data2$`21001_0_0`)) # cogemos esta

data2 = data2 %>%
  dplyr::mutate(bmi_group = case_when(
    bmi < 18.5 ~ 0,
    bmi >= 18.5 & bmi < 25 ~ 1,
    bmi >= 25 & bmi < 30 ~ 2,
    bmi >= 30 & bmi < 35 ~ 3,
    bmi >= 35 & bmi < 40 ~ 4,
    bmi >= 40 ~ 5, 
    is.na(bmi) ~ 9
  ))

data2$bmi_group %>% freq()
# Frequency  Percent
# 0          2625   0.5228
# 1        162261  32.3143
# 2        211973  42.2144
# 3         87501  17.4258
# 4         24973   4.9734
# 5          9698   1.9314
# 9          3103   0.6180
# Total    502134 100.0000

# Socio-demographic covariates -----

## Ethnicity -----
data2 = data2 %>%
  dplyr::rename(ethnicity0 = `21000_0_0`, ethnicity1 = `21000_1_0`,
                ethnicity2 = `21000_2_0`, ethnicity3 = `21000_3_0`) %>% 
  dplyr::mutate(across(starts_with("ethnicity"), ~ case_when(
    . %in% c("1", "1001", "1002", "1003") ~ 1,
    . %in% c("4", "4001", "4002", "4003") ~ 2, 
    . %in% c("3", "3001", "3002", "3003", "3004", "5") ~ 3,
    . %in% c("2", "2001",  "2002",  "2003", "2004", "6") ~ 4,
    . %in% c("-3", "-1") ~ 9,  
    TRUE ~ .))) %>% 
  dplyr::mutate(ethnicity = coalesce(ethnicity0, ethnicity1, 
                                     ethnicity2, ethnicity3), 
                ethnicity = if_else(is.na(ethnicity), 9, ethnicity)
  )

data2$ethnicity %>% freq()
# Frequency  Percent
# 1        472381  94.0747
# 2          8048   1.6028
# 3         11443   2.2789
# 4          7502   1.4940
# 9          2760   0.5497
# Total    502134 100.0000


## Townsend index of deprivation -----
data2 = data2 %>%
  dplyr::mutate(townsend = `22189_0_0`, 
                # Create quintiles for Townsend (5 quantiles)
                townsend_q5 = ntile(townsend, 5),
                # Recode missing values to 9 (replace NA with 9)
                townsend_q5 = replace_na(townsend_q5, 9L))

data2$townsend_q5 %>% freq()
# Frequency  Percent
# 1        100302  19.9751
# 2        100302  19.9751
# 3        100302  19.9751
# 4        100302  19.9751
# 5        100302  19.9751
# 9           624   0.1243
# Total    502134 100.0000

# library(ggplot2)
# data2 %>% ggplot(aes(x = factor(townsend_q5), y = townsend)) + geom_boxplot()

## Education groups ------
data2 = data2 %>%
dplyr::mutate(`6138_0_0` = coalesce(`6138_0_0`, `6138_1_0`, `6138_2_0`, `6138_3_0`)) %>% 
  dplyr::rename(edu_group0 = `6138_0_0`) %>% 
  dplyr::mutate(across(starts_with("edu_group0"), 
                       ~ str_split(.x, "\\|"))) %>%
  unnest_wider(where(is.list), names_sep = "_") %>% 
  dplyr::mutate(education_group =  case_when(
    edu_group0_1 %in% c(1, 6) | edu_group0_2 %in% c(1, 6) | 
      edu_group0_3 %in% c(1, 6) | edu_group0_4 %in% c(1, 6) | 
      edu_group0_5 %in% c(1, 6) | edu_group0_6 %in% c(1, 6) ~ 1, 
    edu_group0_1 %in% c(2, 3, 4) | edu_group0_2 %in% c(2, 3, 4) | 
      edu_group0_3 %in% c(2, 3, 4) | edu_group0_4 %in% c(2, 3, 4) | 
      edu_group0_5 %in% c(2, 3, 4) | edu_group0_6 %in% c(2, 3, 4) ~ 2, 
    edu_group0_1 %in% c(5) | edu_group0_2 %in% c(5) | 
      edu_group0_3 %in% c(5) | edu_group0_4 %in% c(5) | 
      edu_group0_5 %in% c(5) | edu_group0_6 %in% c(5) ~ 3, 
    edu_group0_1 %in% c(-7) | edu_group0_2 %in% c(-7) | 
      edu_group0_3 %in% c(-7) | edu_group0_4 %in% c(-7) | 
      edu_group0_5 %in% c(-7) | edu_group0_6 %in% c(-7) ~ 4, 
    edu_group0_1 %in% c(-3) | edu_group0_2 %in% c(-3) | 
      edu_group0_3 %in% c(-3) | edu_group0_4 %in% c(-3) | 
      edu_group0_5 %in% c(-3) | edu_group0_6 %in% c(-3) ~ 9, 
    is.na(edu_group0_1) |  is.na(edu_group0_2) |  is.na(edu_group0_3) | 
      is.na(edu_group0_4) | is.na(edu_group0_5) |  is.na(edu_group0_6) ~ 9
  ))



data2$education_group %>% freq()
#          Frequency Percent
# 1        234105  46.622
# 2        144287  28.735
# 3         29450   5.865
# 4         85310  16.989
# 9          8982   1.789
# Total    502134 100.000

## Region ----
# 1 "East Midlands" 
# 2 "London" 
# 3 "North East" 
# 4 "North West" 
# 5 "Scotland" 
# 6 "South East" 
# 7 "South West"
# 8 "Wales" 
# 9 "West Midlands" 
# 10 "York

data2 = data2 %>% 
  dplyr::mutate(region = case_when(
    `54_0_0` %in% c(11013) ~ 1, 
    `54_0_0` %in% c(11012, 11018, 11020) ~ 2, 
    `54_0_0` %in% c(11009, 11017) ~ 3, 
    `54_0_0` %in% c(10003, 11001, 11008, 11016) ~ 4, 
    `54_0_0` %in% c(11004, 11005) ~ 5, 
    `54_0_0` %in% c(11002, 11007) ~ 6, 
    `54_0_0` %in% c(11011) ~ 7, 
    `54_0_0` %in% c(11003, 11022, 11023) ~ 8, 
    `54_0_0` %in% c(11006, 11021) ~ 9, 
    `54_0_0` %in% c(11010, 11014) ~ 10))

data2$region %>% freq()
# Frequency Percent
# 1         33853   6.742
# 2         68753  13.692
# 3         58253  11.601
# 4         78799  15.693
# 5         35830   7.136
# 6         43429   8.649
# 7         42977   8.559
# 8         20803   4.143
# 9         44898   8.941
# 10        74539  14.844
# Total    502134 100.000

## Country (alternative to region) -----
data2 = data2 %>% 
  dplyr::mutate(country = case_when(
    region == 5 ~ 1, # Scotland
    region == 8 ~ 2, # Wales
    TRUE ~ 3)) # England

data2$country %>% freq()
# Frequency Percent
# 1         35830   7.136
# 2         20803   4.143
# 3        445501  88.722
# Total    502134 100.000



# Lifestyle covariates -----
## A) Smoking status -----
# - I will use any info I have, giving priority to more recent data
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20116
# 20116: Smoking status
# -3	Prefer not to answer
# 0	Never
# 1	Previous
# 2	Current

# Recode smoking variables and create the final smoking status variable
data2 = data2 %>%
  # Recode and create smoke_* variables
  mutate(smoking_status = case_when(
    `20116_0_0` == -3 ~ 9,
    `20116_0_0` == 0 ~ 0,  # Never
    `20116_0_0` == 1 ~ 1,  # Previous
    `20116_0_0` == 2 ~ 2,  # Current
    TRUE ~ 9))

# Show the distribution of the final smoking status
data2$smoking_status %>% freq()
# Frequency  Percent
# 0        273328  54.4333
# 1        172921  34.4372
# 2         52937  10.5424
# 9          2948   0.5871
# Total    502134 100.0000


## B) Physical acticity -----
# - I will use any info I have, taking the average of all timepoints
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=904
# 904 : No. days/week of vigorous PA
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=884
# 884 : No. days/week of moderate PA
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=864
# 864 : No. days/week walked 10+ mins
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=914
# 914 : Duration of vigorous PA
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=914
# 894 : Duration of moderate PA 
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=914
# 874 : Duration of walks

# Calculate MET minutes
# Walking MET-minutes/week = 3.3 *  walking minutes * walking days
# Moderate MET-minutes/week = 4.0 * moderate minutes *  moderate days
# Vigorous MET-minutes/week = 8.0 *  vigorous minutes *   vigorous days

# Step 1 - Recode walking days, duration, moderate activity, and vigorous activity
data2 = data2 %>%
  # Recoding walking days (-3, -1 to NA; -2 to 0)
  dplyr::mutate(
    walkdays = case_when(
      `864_0_0` %in% c(-3, -1) ~ NA_real_,
      `864_0_0` == -2 ~ 0,
      TRUE ~ `864_0_0`),
    # Recoding walking duration (-3, -1 to NA and setting to 0 if walking days is 0)
    walkdur = case_when(
      `874_0_0` %in% c(-3, -1) ~ NA_real_,
      `864_0_0` == 0 ~ 0,
      TRUE ~ `874_0_0`), 
    # Recoding moderate activity days (-3, -1 to NA)
    modactdays = case_when(
      `884_0_0` %in% c(-3, -1) ~ NA_real_, 
      TRUE ~ `884_0_0`),
    # Recoding moderate activity duration (-3, -1 to NA and setting to 0 if moderate activity days is 0)
    modactdur = case_when(
      `894_0_0` %in% c(-3, -1) ~ NA_real_,
      `884_0_0` == 0 ~ 0,
      TRUE ~ `894_0_0`), 
    # Recoding vigorous activity days (-3, -1 to NA)
    vigactdays = case_when(
      `904_0_0` %in% c(-3, -1) ~ NA_real_,
      TRUE ~ `904_0_0`),
    # Recoding vigorous activity duration (-3, -1 to NA and setting to 0 if vigorous activity days is 0)
    vigactdur = case_when(
      `914_0_0` %in% c(-3, -1) ~ NA_real_,
      `904_0_0` == 0 ~ 0,
      TRUE ~ `914_0_0`))

# Step 2 - Calculate row means for all timepoints
data2 = data2 %>%
  # Step 3 - Obtain mins/week by multiplying duration and number of days
  dplyr::mutate(minsperweek_light = walkdur * walkdays,
                minsperweek_mod = modactdur * modactdays,
                minsperweek_vig = vigactdur * vigactdays) %>%
  # Step 4 - Convert mins/week into mins/day
  dplyr::mutate(minsperday_light = minsperweek_light / 7,
                minsperday_mod = minsperweek_mod / 7,
                minsperday_vig = minsperweek_vig / 7) %>%
  
  # Step 5 - Calculate METs per category of exercise
  dplyr::mutate(MET_light = minsperday_light * 3.3,
                MET_mod = minsperday_mod * 4.0,
                MET_vig = minsperday_vig * 8.0) %>%
  
  # Step 6 - Convert MET-mins/day to MET-hours/week
  dplyr::mutate(MET_light_hrwk = (MET_light * 7) / 60,
                MET_mod_hrwk = (MET_mod * 7) / 60,
                MET_vig_hrwk = (MET_vig * 7) / 60, 
                #  IPAQ Recommendations: remove activities less than 10 minutes
                walks_duration_minrecode = if_else(walkdur <= 9, 0, 1), 
                mpa_duration_minrecode = if_else(modactdur <= 9, 0, 1),
                vpa_duration_minrecode = if_else(vigactdur <= 9, 0, 1), 
                MET_light_hrwk = if_else(walks_duration_minrecode == 0, 0, MET_light_hrwk),
                MET_mod_hrwk = if_else(mpa_duration_minrecode == 0, 0, MET_mod_hrwk),
                MET_vig_hrwk = if_else(vpa_duration_minrecode == 0, 0, MET_vig_hrwk)) 
# Step 7 - Calculate total MET-hours/week
data2 = data2 %>%
  dplyr::mutate(total_MET_hrwk = if_else(
    is.na(MET_light_hrwk) & is.na(MET_mod_hrwk) & is.na(MET_vig_hrwk), NA_real_,
    rowSums(select(., MET_light_hrwk, MET_mod_hrwk, MET_vig_hrwk), na.rm = TRUE))) %>%
  # Step 8 - Recode total MET-hours/week into meaningful categories
  dplyr::mutate(pa_group = case_when(
    is.na(total_MET_hrwk) ~ 9, 
    total_MET_hrwk < 10 ~ 1,
    total_MET_hrwk >= 10 & total_MET_hrwk < 50 ~ 2,
    total_MET_hrwk >= 50 ~ 3)) %>%
  # Step 9 - Categorize PA based on recommendations
  dplyr::mutate(pa_rec = case_when(
    minsperweek_vig >= 75 ~ 1,
    minsperweek_mod >= 150 & minsperweek_vig < 75 ~ 1,
    is.na(total_MET_hrwk) ~ 9,
    TRUE ~ 0))

# Summary of the final variables for validation
data2$pa_group %>% freq()
# Frequency Percent
# 1        114561  22.815
# 2        232195  46.242
# 3        134530  26.792
# 9         20848   4.152
# Total    502134 100.000
data2$pa_rec %>%  freq()
#        Frequency Percent
# 0        249608  49.709
# 1        231678  46.139
# 9         20848   4.152
# Total    502134 100.000


## C) Alcohol -----
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=1568
# 1568: Average weekly red wine intake 
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=1588
# 1588: Average weekly beer plus cider intake
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=1578
# 1578: Average weekly champagne plus white wine intake
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=1598
# 1598: Average weekly spirits intake
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=1608
# 1608: Average weekly fortified wine intake
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=5364
# 5364: Average weekly intake of other alcoholic drinks

# Step 1 - Recode alcohol intake variables for weekly data (-3, -1 to NA_real_)
data2 = data2 %>%
  dplyr::mutate(`1568_0_0` =  if_else(`1568_0_0` %in% c(-3, -1), NA_real_, `1568_0_0`),
                `1588_0_0` = if_else(`1588_0_0` %in% c(-3, -1), NA_real_,  `1588_0_0`),
                `1578_0_0` =  if_else(`1578_0_0` %in% c(-3, -1), NA_real_, `1578_0_0`),
                `1598_0_0` = if_else(`1598_0_0` %in% c(-3, -1), NA_real_, `1598_0_0`),
                `1608_0_0` =  if_else(`1608_0_0` %in% c(-3, -1), NA_real_, `1608_0_0`),
                `5364_0_0` = if_else(`5364_0_0` %in% c(-3, -1), NA_real_, `5364_0_0`), 
                `4407_0_0` =  if_else(`4407_0_0` %in% c(-3, -1), NA_real_, `4407_0_0`),
                `4429_0_0` =  if_else(`4429_0_0` %in% c(-3, -1), NA_real_, `4429_0_0`),
                `4418_0_0` =  if_else(`4418_0_0` %in% c(-3, -1), NA_real_, `4418_0_0`),
                `4440_0_0` =  if_else(`4440_0_0` %in% c(-3, -1), NA_real_, `4440_0_0`),
                `4451_0_0` =  if_else(`4451_0_0` %in% c(-3, -1), NA_real_, `4451_0_0`),
                `4462_0_0` =  if_else(`4462_0_0` %in% c(-3, -1), NA_real_, `4462_0_0`))

# Step 2 - Multiply alcohol weekly units by conversion factors
data2 = data2 %>%
  dplyr::mutate(rw_u_wk = `1568_0_0` * 2,
                bc_u_wk =`1588_0_0` * 2.5,
                cw_u_wk =`1578_0_0` * 2,
                sp_u_wk =`1598_0_0` * 1,
                fw_u_wk =`1608_0_0` * 2, 
                other_u_wk =`5364_0_0` * 1) 
# Step 3 - Multiply alcohol monthly units by conversion factors and divide by 4.3 to convert to weekly
data2 = data2 %>%
  dplyr::mutate(rw_u_mth = (`4407_0_0` * 2) / 4.3,
                bc_u_mth = (`4429_0_0` * 2.5) / 4.3,
                cw_u_mth = (`4418_0_0` * 2) / 4.3,
                sp_u_mth = (`4440_0_0` * 1) / 4.3, 
                fw_u_mth = (`4451_0_0` * 2) / 4.3,
                other_u_mth = (`4462_0_0` * 1) / 4.3)

# Step 5 - Generate total alcohol intake variable
data2 = data2 %>%
  dplyr::mutate(rw_u_wk_def = coalesce(rw_u_wk, rw_u_mth), 
                bc_u_wk_def = coalesce(bc_u_wk, bc_u_mth), 
                cw_u_wk_def = coalesce(cw_u_wk, cw_u_mth), 
                sp_u_wk_def = coalesce(sp_u_wk, sp_u_mth), 
                fw_u_wk_def = coalesce(fw_u_wk, fw_u_mth), 
                other_u_wk_def = coalesce(other_u_wk, other_u_mth))

data2 = data2 %>%
  dplyr::mutate(alcohol_intake = if_else(
    is.na(rw_u_wk_def) & is.na(bc_u_wk_def) & is.na(cw_u_wk_def) & 
      is.na(sp_u_wk_def) & is.na(fw_u_wk_def) & is.na(other_u_wk_def) , NA_real_,
    rowSums(select(., rw_u_wk_def, bc_u_wk_def, cw_u_wk_def, sp_u_wk_def, 
                   fw_u_wk_def, other_u_wk_def), na.rm = TRUE)))

# Step 6 - Categorize alcohol intake into groups (none, occasional, moderate, heavy, unknown)
data2 = data2 %>%
  dplyr::mutate(alc_group = case_when(
    alcohol_intake == 0 ~ 0,
    alcohol_intake <= 1 ~ 1,
    alcohol_intake > 1 & alcohol_intake <= 14 ~ 2,
    alcohol_intake > 14 ~ 3,
    is.na(alcohol_intake) ~ 9))

# Final validation: Summary of alcohol intake groups
data2$alc_group %>% freq()
#        Frequency Percent
# 0         10622   2.115
# 1          9764   1.945
# 2        164228  32.706
# 3        202194  40.267
# 9        115326  22.967
# Total    502134 100.000
data2$alcohol_intake %>% summary()

# Dietary exposures from FFQ -----

## A) Fruits -----
# Fruit/veg variables have specific coding: -10 (less than one),
# -3 (prefer not to answer), -1 (do not know)
# 1 portion of cooked (1289)/raw veg (1299) = 3 tbsp/day
# 1 portion of fresh fruit (1309) = 1 piece
# Portion sizes based on: 
#   https://www.nhs.uk/live-well/eat-well/5-a-day-portion-sizes/#5-a-day-fruit-portions

# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=1289
# 1289: Cooked vegetable intake
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=1299
# 1200: Salad / raw vegetable intake
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=1309
# 1309: Fresh fruit intake

# Step 1 - Recoding cooked vegetable intake, raw vegetable intake and 
# fruit intake
data2 = data2 %>%
  dplyr::mutate(cooked_veg  = case_when(
    `1289_0_0` == -10 ~ 0.5,  # recode -10 to 0.5
    `1289_0_0`  %in% c(-3, -1) ~ NA_real_,  # recode -3 and -1 to NA
    TRUE ~  `1289_0_0`), 
    raw_veg  = case_when(
      `1299_0_0` == -10 ~ 0.5,  # recode -10 to 0.5
      `1299_0_0`  %in% c(-3, -1) ~ NA_real_,  # recode -3 and -1 to NA
      TRUE ~  `1299_0_0`), 
    fruit  = case_when(
      `1309_0_0` == -10 ~ 0.5,  # recode -10 to 0.5
      `1309_0_0`  %in% c(-3, -1) ~ NA_real_,  # recode -3 and -1 to NA
      TRUE ~  `1309_0_0`),
  )

# Step 2 - Generate fruit-veg portions per day
data2 = data2 %>%
  dplyr::mutate(fruit_veg = rowSums(select(., cooked_veg, raw_veg, fruit)), 
                fruit_veg = round(fruit_veg, digits = 0),
                fruit_veg_3 = ntile(fruit_veg, 3), 
                fruit_veg_3 = as.integer(fruit_veg_3),
                fruit_veg_3 = dplyr::case_when(
                  is.na(fruit) & is.na(raw_veg) & is.na(cooked_veg) ~ 9L,
                  fruit_veg == 0 ~ 0L,
                  is.na(fruit_veg_3) ~ 9L,
                  TRUE ~ fruit_veg_3
                ))


# Step 3 - Indicate whether meeting the fruit-veg recommendations
data2 = data2 %>%
  dplyr::mutate(fruit_veg_rec = case_when(
    fruit_veg >= 5 ~ 1,  # meeting recommendation
    fruit_veg < 5 ~ 0,   # not meeting recommendation
    TRUE ~ 9))  # missing values recoded as 9

# Summary of the final variables for validation
# data2$fruit_veg %>% freq()
data2$fruit_veg_3 %>% freq()
#        Frequency  Percent
# 0          3223   0.6419
# 1        160244  31.9126
# 2        163467  32.5545
# 3        163467  32.5545
# 9         11733   2.3366
# Total    502134 100.0000
data2$fruit_veg_rec %>% freq()
# Frequency Percent
# 0        111742  22.253
# 1        378659  75.410
# 9         11733   2.337
# Total    502134 100.000


## B)  Saturated fats -----
# Step 1 - Recoding intake categories for beef, lamb, pork, and cheese
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=1369
# 1369: Beef intake
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=1379
# 1379: Lamb/mutton intake
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=1389
# 1389: Pork intake
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=1408
# 1408: Cheese intake
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=1428
# 1428: Spread type


data2 = data2 %>%
  dplyr::mutate(beef = case_when(
    `1369_0_0` %in% c(-3, -1) ~ NA_real_,  # recode missing values
    `1369_0_0` == 0 ~ 0,  # none
    `1369_0_0` == 1 ~ 1,  # little
    `1369_0_0` %in% c(2, 3) ~ 2,  # moderate
    `1369_0_0` %in% c(4, 5) ~ 3,  # high
    TRUE ~ `1369_0_0`), 
    lamb = case_when(
      `1379_0_0` %in% c(-3, -1) ~ NA_real_,  # recode missing values
      `1379_0_0` == 0 ~ 0,  # none
      `1379_0_0` == 1 ~ 1,  # little
      `1379_0_0` %in% c(2, 3) ~ 2,  # moderate
      `1379_0_0` %in% c(4, 5) ~ 3,  # high
      TRUE ~ `1379_0_0`), 
    pork = case_when(
      `1389_0_0` %in% c(-3, -1) ~ NA_real_,  # recode missing values
      `1389_0_0` == 0 ~ 0,  # none
      `1389_0_0` == 1 ~ 1,  # little
      `1389_0_0` %in% c(2, 3) ~ 2,  # moderate
      `1389_0_0` %in% c(4, 5) ~ 3,  # high
      TRUE ~ `1389_0_0`), 
    cheese = case_when(
      `1408_0_0` %in% c(-3, -1) ~ NA_real_,  # recode missing values
      `1408_0_0` == 0 ~ 0,  # none
      `1408_0_0` == 1 ~ 1,  # little
      `1408_0_0` %in% c(2, 3) ~ 2,  # moderate
      `1408_0_0` %in% c(4, 5) ~ 3,  # high
      TRUE ~ `1408_0_0`),
    spread = case_when(
      `1428_0_0` %in% c(-3, -1) ~ NA_real_,  # recode missing values
      `1428_0_0` == 0 ~ 0,  # none
      `1428_0_0` == 1 ~ 1,  # little
      `1428_0_0` %in% c(2, 3) ~ 2,  # moderate
      `1428_0_0` %in% c(4, 5) ~ 3,  # high
      TRUE ~ `1428_0_0`), 
  ) 

# Step 2 - Get the saturated fat score by summing up all intake variables
data2 = data2 %>%
  dplyr::mutate(satfat_score = rowSums(select(., beef, lamb, pork, cheese, spread), na.rm = TRUE))

# Step 3 - Categorize the saturated fat score into quantiles
data2 = data2 %>%
  dplyr::mutate(satfat_3 = ntile(satfat_score, 3)) %>%
  # Recode missing values into the unknown category
  dplyr::mutate(satfat_3 = case_when(
    is.na(beef) & is.na(lamb) & is.na(pork) & is.na(cheese) & is.na(spread) ~ 9L,  # unknown
    TRUE ~ satfat_3
  ))


data2$satfat_3 %>% freq()
# Frequency  Percent
# 1        166230  33.1178
# 2        167312  33.3333
# 3        167312  33.3333
# 9          1082   0.2156
# Total    501936 100.0000

# Medical history -----
## Menopausal status  -------

# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=2724
# 2724: Had menopause

# Recode variables for menopausal status
data2 = data2 %>%
  dplyr::mutate(
    menopause = case_when(`2724_0_0` == 1 ~ 1,
                          `2724_0_0` == 0 ~ 0,
                          `2724_0_0`%in% c(2, 3, -3) ~ 0,
                          is.na(`2724_0_0`) ~ 9)) 

data2$menopause %>% freq()
#        Frequency Percent
# 0        107382   21.39
# 1        165301   32.92
# 9        229451   45.70
# Total    502134  100.00

# Family history ------

## A) CVD - includes Heart Disease, Stroke, high BP ------

# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20107
# 20107: Illnesses of father
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20110
# 20110: Illnesses of mother
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20111
# 20111: Illnesses of siblings


data2 = data2 %>% 
  dplyr::mutate(illness_father = `20107_0_0`, illness_mother = `20110_0_0`, 
                illness_siblings = `20111_0_0`) 

data2 = data2 %>% 
  dplyr::mutate(across(starts_with("illness_") , 
                       ~ str_split(.x, "\\|"))) %>% 
  unnest_wider(where(is.list), names_sep = "_") 

data2 = data2 %>% 
  dplyr::mutate(across(starts_with("illness_") , 
                       ~ as.numeric(.))) %>% 
  dplyr::mutate(across(starts_with("illness_") , 
                       ~ case_when(
                         . %in% c(1, 2, 8) ~ 1,
                         is.na(.) ~ NA_real_, 
                         TRUE ~ 0), .names = "cvd_{col}")) 

data2 = data2 %>% 
  dplyr::mutate(
    cvd_fam = case_when(
      cvd_illness_father_1 == 1 | cvd_illness_father_2 == 1 |  cvd_illness_father_3 == 1 | 
      cvd_illness_father_4 == 1 | cvd_illness_father_5 == 1 | cvd_illness_father_6 == 1 | 
      cvd_illness_father_7 == 1 | cvd_illness_father_8 == 1 | cvd_illness_father_9 == 1 | 
      cvd_illness_father_10 == 1 | 
      cvd_illness_mother_1 == 1 | cvd_illness_mother_2 == 1 |  cvd_illness_mother_3 == 1 | 
      cvd_illness_mother_4 == 1 | cvd_illness_mother_5 == 1 | cvd_illness_mother_6 == 1 | 
      cvd_illness_mother_7 == 1 | cvd_illness_mother_8 == 1 | cvd_illness_mother_9 == 1 | 
      cvd_illness_mother_10 == 1 | cvd_illness_mother_11 == 1 | 
      cvd_illness_siblings_1 == 1 | cvd_illness_siblings_2 == 1 |  cvd_illness_siblings_3 == 1 | 
      cvd_illness_siblings_4 == 1 | cvd_illness_siblings_5 == 1 | cvd_illness_siblings_6 == 1 | 
      cvd_illness_siblings_7 == 1 | cvd_illness_siblings_8 == 1 | cvd_illness_siblings_9 == 1 | 
      cvd_illness_siblings_10 == 1  ~ 1, 
      cvd_illness_father_1 == 0 | cvd_illness_father_2 == 0 |  cvd_illness_father_3 == 0 | 
      cvd_illness_father_4 == 0 | cvd_illness_father_5 == 0 | cvd_illness_father_6 == 0 | 
      cvd_illness_father_7 == 0 | cvd_illness_father_8 == 0 | cvd_illness_father_9 == 0 | 
      cvd_illness_father_10 == 0 | 
      cvd_illness_mother_1 == 0 | cvd_illness_mother_2 == 0 |  cvd_illness_mother_3 == 0 | 
      cvd_illness_mother_4 == 0 | cvd_illness_mother_5 == 0 | cvd_illness_mother_6 == 0 | 
      cvd_illness_mother_7 == 0 | cvd_illness_mother_8 == 0 | cvd_illness_mother_9 == 0 | 
      cvd_illness_mother_10 == 0 | cvd_illness_mother_11 == 0 | 
      cvd_illness_siblings_1 == 0 | cvd_illness_siblings_2 == 0 |  cvd_illness_siblings_3 == 0 | 
      cvd_illness_siblings_4 == 0 | cvd_illness_siblings_5 == 0 | cvd_illness_siblings_6 == 0 | 
      cvd_illness_siblings_7 == 0 | cvd_illness_siblings_8 == 0 | cvd_illness_siblings_9 == 0 | 
      cvd_illness_siblings_10 == 0 | cvd_illness_siblings_11 == 0 | cvd_illness_siblings_12 == 0  ~ 0, 
      TRUE ~ NA_real_)) 

data2$cvd_fam %>% freq()  
# Frequency Percent Valid Percent
# 0        125072  24.908         25.34
# 1        368496  73.386         74.66
# NA's       8566   1.706              
# Total    502134 100.000        100.00

## B) Diabetes ------

data2 = data2 %>% 
  dplyr::mutate(across(starts_with("illness_") , 
                       ~ case_when(
                         . == 9 ~ 1,
                         is.na(.) ~ NA_real_, 
                         TRUE ~ 0),
                       .names = "diab_{col}")) %>% 
  dplyr::mutate(
    diab_fam = case_when(
      diab_illness_father_1 == 1 | diab_illness_father_2 == 1 |  diab_illness_father_3 == 1 | 
        diab_illness_father_4 == 1 | diab_illness_father_5 == 1 | diab_illness_father_6 == 1 | 
        diab_illness_father_7 == 1 | diab_illness_father_8 == 1 | diab_illness_father_9 == 1 | 
        diab_illness_father_10 == 1 | 
        diab_illness_mother_1 == 1 | diab_illness_mother_2 == 1 |  diab_illness_mother_3 == 1 | 
        diab_illness_mother_4 == 1 | diab_illness_mother_5 == 1 | diab_illness_mother_6 == 1 | 
        diab_illness_mother_7 == 1 | diab_illness_mother_8 == 1 | diab_illness_mother_9 == 1 | 
        diab_illness_mother_10 == 1 | diab_illness_mother_11 == 1 | 
        diab_illness_siblings_1 == 1 | diab_illness_siblings_2 == 1 |  diab_illness_siblings_3 == 1 | 
        diab_illness_siblings_4 == 1 | diab_illness_siblings_5 == 1 | diab_illness_siblings_6 == 1 | 
        diab_illness_siblings_7 == 1 | diab_illness_siblings_8 == 1 | diab_illness_siblings_9 == 1 | 
        diab_illness_siblings_10 == 1  ~ 1, 
      diab_illness_father_1 == 0 | diab_illness_father_2 == 0 |  diab_illness_father_3 == 0 | 
        diab_illness_father_4 == 0 | diab_illness_father_5 == 0 | diab_illness_father_6 == 0 | 
        diab_illness_father_7 == 0 | diab_illness_father_8 == 0 | diab_illness_father_9 == 0 | 
        diab_illness_father_10 == 0 | 
        diab_illness_mother_1 == 0 | diab_illness_mother_2 == 0 |  diab_illness_mother_3 == 0 | 
        diab_illness_mother_4 == 0 | diab_illness_mother_5 == 0 | diab_illness_mother_6 == 0 | 
        diab_illness_mother_7 == 0 | diab_illness_mother_8 == 0 | diab_illness_mother_9 == 0 | 
        diab_illness_mother_10 == 0 | diab_illness_mother_11 == 0 | 
        diab_illness_siblings_1 == 0 | diab_illness_siblings_2 == 0 |  diab_illness_siblings_3 == 0 | 
        diab_illness_siblings_4 == 0 | diab_illness_siblings_5 == 0 | diab_illness_siblings_6 == 0 | 
        diab_illness_siblings_7 == 0 | diab_illness_siblings_8 == 0 | diab_illness_siblings_9 == 0 | 
        diab_illness_siblings_10 == 0 | diab_illness_siblings_11 == 0 | diab_illness_siblings_12 == 0  ~ 0, 
      TRUE ~ NA_real_))


data2$diab_fam %>% freq()
# Frequency Percent Valid Percent
# 0        384903  76.684         78.01
# 1        108469  21.610         21.99
# NA's       8564   1.706              
# Total    501936 100.000        100.00

## C) Hypertension only ------

data2 = data2 %>% 
  dplyr::mutate(across(starts_with("illness_") , 
                       ~ case_when(
                         . == 8 ~ 1,
                         is.na(.) ~ NA_real_, 
                         TRUE ~ 0),
                       .names = "htn_{col}")) %>% 
  dplyr::mutate(
    htn_fam = case_when(
      htn_illness_father_1 == 1 | htn_illness_father_2 == 1 |  htn_illness_father_3 == 1 | 
        htn_illness_father_4 == 1 | htn_illness_father_5 == 1 | htn_illness_father_6 == 1 | 
        htn_illness_father_7 == 1 | htn_illness_father_8 == 1 | htn_illness_father_9 == 1 | 
        htn_illness_father_10 == 1 | 
        htn_illness_mother_1 == 1 | htn_illness_mother_2 == 1 |  htn_illness_mother_3 == 1 | 
        htn_illness_mother_4 == 1 | htn_illness_mother_5 == 1 | htn_illness_mother_6 == 1 | 
        htn_illness_mother_7 == 1 | htn_illness_mother_8 == 1 | htn_illness_mother_9 == 1 | 
        htn_illness_mother_10 == 1 | htn_illness_mother_11 == 1 | 
        htn_illness_siblings_1 == 1 | htn_illness_siblings_2 == 1 |  htn_illness_siblings_3 == 1 | 
        htn_illness_siblings_4 == 1 | htn_illness_siblings_5 == 1 | htn_illness_siblings_6 == 1 | 
        htn_illness_siblings_7 == 1 | htn_illness_siblings_8 == 1 | htn_illness_siblings_9 == 1 | 
        htn_illness_siblings_10 == 1  ~ 1, 
      htn_illness_father_1 == 0 | htn_illness_father_2 == 0 |  htn_illness_father_3 == 0 | 
        htn_illness_father_4 == 0 | htn_illness_father_5 == 0 | htn_illness_father_6 == 0 | 
        htn_illness_father_7 == 0 | htn_illness_father_8 == 0 | htn_illness_father_9 == 0 | 
        htn_illness_father_10 == 0 | 
        htn_illness_mother_1 == 0 | htn_illness_mother_2 == 0 |  htn_illness_mother_3 == 0 | 
        htn_illness_mother_4 == 0 | htn_illness_mother_5 == 0 | htn_illness_mother_6 == 0 | 
        htn_illness_mother_7 == 0 | htn_illness_mother_8 == 0 | htn_illness_mother_9 == 0 | 
        htn_illness_mother_10 == 0 | htn_illness_mother_11 == 0 | 
        htn_illness_siblings_1 == 0 | htn_illness_siblings_2 == 0 |  htn_illness_siblings_3 == 0 | 
        htn_illness_siblings_4 == 0 | htn_illness_siblings_5 == 0 | htn_illness_siblings_6 == 0 | 
        htn_illness_siblings_7 == 0 | htn_illness_siblings_8 == 0 | htn_illness_siblings_9 == 0 | 
        htn_illness_siblings_10 == 0 | htn_illness_siblings_11 == 0 | htn_illness_siblings_12 == 0  ~ 0, 
      TRUE ~ NA_real_))


data2$htn_fam %>% freq()
# Frequency Percent Valid Percent
# 0        254447  50.673         51.55
# 1        239121  47.621         48.45
# NA's       8566   1.706              
# Total    502134 100.000        100.00

# saveRDS(data2, "/mnt/project/2025_04_13_covariates.RDS")

# New script 2025/07/21 to calculate medication variables 
# Medication ------
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=6177
# https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=100626
# 1	Cholesterol lowering medication
# 2	Blood pressure medication
# 3	Insulin
# 4	Hormone replacement therapy
# 5	Oral contraceptive pill or minipill
# -7	None of the above
# -1	Do not know
# -3	Prefer not to answer
## Medication for cholesterol -----
## Medication for blood pressure -----
## Medication for diabetes -----
## Medication for others reasons -----

data3 = data2 %>% 
  dplyr::rename(medication_male = `6177_0_0`, medication_female = `6153_0_0`) %>% 
  dplyr::mutate(across(starts_with("medication_") , 
                       ~ str_split(.x, "\\|"))) %>% 
  unnest_wider(where(is.list), names_sep = "_") 

data3 = data3 %>% 
  dplyr::mutate(across(starts_with("medication_") , 
                       ~ as.numeric(.x))) %>% 
  dplyr::mutate(
    medication_cholesterol = case_when(
     medication_male_1 == 1 | medication_male_2 == 1 | medication_male_3 == 1 | 
     medication_female_1 == 1 | medication_female_2 == 1 | medication_female_3 == 1 | 
     medication_female_4 == 1  ~ 1, TRUE ~ 0), 
    medication_bp = case_when(
      medication_male_1 == 2 | medication_male_2 == 2 | medication_male_3 == 2 | 
      medication_female_1 == 2 | medication_female_2 == 2 | medication_female_3 == 2 | 
      medication_female_4 == 2  ~ 2, TRUE ~ 0), 
    medication_diabetes = case_when(
      medication_male_1 == 3 | medication_male_2 == 3 | medication_male_3 == 3 | 
      medication_female_1 == 3 | medication_female_2 == 3 | medication_female_3 == 3 | 
      medication_female_4 == 3  ~ 1, TRUE ~ 0),
    medication_others = case_when(
        medication_female_1 %in% c(4, 5) | medication_female_2 %in% c(4, 5) | 
        medication_female_3 %in% c(4, 5) | medication_female_4 %in% c(4, 5)  ~ 1, TRUE ~ 0),
  )

data3$medication_cholesterol %>% freq()
# Frequency Percent
# 0        415280    82.7
# 1         86854    17.3
# Total    502134   100.0

data3$medication_bp %>% freq()
# Frequency Percent
# 0        398180    79.3
# 2        103954    20.7
# Total    502134   100.0
data3$medication_diabetes %>% freq()
# Frequency Percent
# 0        496525  98.883
# 1          5609   1.117
# Total    502134 100.000
data3$medication_others %>% freq()
# Frequency Percent
# 0        475681  94.732
# 1         26453   5.268
# Total    502134 100.000

# dir.create("upload")
# saveRDS(data2, "./upload/2025_07_09_covariates.RDS")

# Upload the files to RAP UK Bibank----
## Locate all the files generated in the upload folder 

# setwd("./upload")

# system("dx ls data/processed")

# system("pwd")

system("dx upload /home/rstudio-server/upload/*")
