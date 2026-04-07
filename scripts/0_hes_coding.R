

#NOTES
# I am using the variable 41270	"Diagnoses - ICD10"	Summary Diagnoses for ICD10
# and variable 41272	"Operative procedures - OPCS4"	Summary Operations for OPCS4 
# The variables 41270 and 41272 store all codes together separated by "|"
# I am using the variable 41280 "Date of first in-patient diagnosis - ICD10"	Summary Diagnoses for date of diagnosis in ICD10 codes 
# and variable 41282	"Date of first operative procedure - OPCS4"	Summary Operations for date of operative proc in OPCS4 codes. 
# The variables 41280 and 41282 have multiple subvariables (e.g., 41280-0.0, 41280-0.1, 41280-0.2, etc.) and I am assuming that the order of these 
#variables match the order of the codes in variables 41270 and 41272. For example, the second ICD10 code listed in variable 41270 corresponds to the 
#date recorded in subvariable 41280-0.1 
# For all variables 0 means absence of the disease or death, and 1 means presence

options(max.print=10000000)
library(readr)
library(dplyr)
library(descr)


# setwd("C:/Users/Carmen/ownCloud - CARMEN_MARÍA PIERNAS SÁNCHEZ@drive.ugr.es/UK Biobank UGR/Raw data (NEW)")  # change for the directory you are using


hes <- read_csv("/mnt/project/data/raw/Health Outcomes.csv")  # this is the path where the dataset is stored


######### CHD - Coronary Heart Disease #############
#function to extract codes of disease
extract_codes <- function(x) {
  if (is.na(x) || x == "") return(character(0))
  unlist(strsplit(x, "\\|"))
}

# target codes
target_chd_ic10 <- c("^I20", "^I21", "^I22", "^I23", "^I24", "^I25")
target_chd_opcs4 <- c("^K40", "^K41", "^K42", "^K43", "^K44", "^K45", "^K46", "^K49", "^K50", "^K75")

# function to find the targeted codes
matches_codes <- function(code, patterns) {
  any(sapply(patterns, function(p) grepl(p, code)))
}

# function to process every line in the dataset
process_row <- function(row) {
  codes_41270 <- extract_codes(row[["41270-0.0"]])
  codes_41272 <- extract_codes(row[["41272-0.0"]])
  
  # checking matched codes 
  matches_41270 <- if (length(codes_41270) > 0) sapply(codes_41270, function(x) matches_codes(x, target_chd_ic10)) else logical(0)
  matches_41272 <- if (length(codes_41272) > 0) sapply(codes_41272, function(x) matches_codes(x, target_chd_opcs4)) else logical(0)
  
  chd_diag <- 0
  earliest_date_chd <- NA
  all_dates <- c()  
  
  if (any(matches_41270) || any(matches_41272)) {
    chd_diag <- 1
    
    # identifying the position(index) of the targeted codes in the list ICD-10 
    if (length(matches_41270) > 0) {
      indices <- which(matches_41270)
      # Finding the correspondent date in variables 41280-0.X, according to the position of the code matched 
      for (index in indices) {
        date_col <- paste0("41280-0.", index - 1)  # because the variables with the dates starts in 0 not in 1
        if (date_col %in% names(row) && !is.na(row[[date_col]]) && row[[date_col]] != "") {
          if (grepl("^\\d{4}-\\d{2}-\\d{2}$", row[[date_col]])) {
            all_dates <- c(all_dates, row[[date_col]]) 
          }
        }
      }
    }
    
    # same process for OPCS4 codes
    if (length(matches_41272) > 0) {
      indices <- which(matches_41272)
      for (index in indices) {
        date_col <- paste0("41282-0.", index - 1)  
        if (date_col %in% names(row) && !is.na(row[[date_col]]) && row[[date_col]] != "") {
          if (grepl("^\\d{4}-\\d{2}-\\d{2}$", row[[date_col]])) {
            all_dates <- c(all_dates, row[[date_col]])
          }
        }
      }
    }
    
    # Selecting the earliest date founded
    if (length(all_dates) > 0) {
      earliest_date_chd <- min(as.Date(all_dates, format = "%Y-%m-%d"), na.rm = TRUE)
    }
  }
  
  return(c(chd = chd_diag, date_chd = as.character(earliest_date_chd)))
}

#applying the results to the dataset
resultado <- t(apply(hes, 1, process_row))
hes$chd <- as.numeric(resultado[, "chd"])
hes$date_chd <- as.Date(resultado[, "date_chd"])


########## CHF - Congestive Heart Failure ############# 
#I did not include "I43.0 Cardiomyopathy in infectious and parasitic diseases classified elsewhere" as in Carmen's do file

extract_codes <- function(x) {
  if (is.na(x) || x == "") return(character(0))
  unlist(strsplit(x, "\\|"))
}

target_chf <- c("I500", "I501", "I509", "I110", "I130", "I132", 
                "^I42", "I431", "I432", "I438")

matches_codes <- function(code, patterns) {
  any(sapply(patterns, function(p) grepl(p, code)))
}

process_row <- function(row) {
   codes_41270 <- extract_codes(row[["41270-0.0"]])
  
  matches_41270 <- sapply(codes_41270, function(x) matches_codes(x, target_chf))

  chf_diag <- 0
  earliest_date_chf <- NA

  if (any(matches_41270)) {
   chf_diag <- 1
    
    indices <- which(matches_41270)
  
    all_dates <- c()
    
    for (index in indices) {
      date_col <- paste0("41280-0.", index-1) 
      if (date_col %in% names(row) && !is.na(row[[date_col]]) && row[[date_col]] != "") {
        if (grepl("^\\d{4}-\\d{2}-\\d{2}$", row[[date_col]])) {
          all_dates <- c(all_dates, row[[date_col]])
        }
      }
    }
    
  if (length(all_dates) > 0) {
      earliest_date_chf <- min(as.Date(all_dates, format="%Y-%m-%d"), na.rm = TRUE)
    }
  }
  
  return(c(chf = chf_diag, date_chf = as.character(earliest_date_chf)))
}

resultado <- t(apply(hes, 1, process_row))
hes$chf <- as.numeric(resultado[, "chf"])
hes$date_chf <- as.Date(resultado[, "date_chf"])


############### STROKE #############
extract_codes <- function(x) {
  if (is.na(x) || x == "") return(character(0))
  unlist(strsplit(x, "\\|"))
}

target_stroke <- c("^I60", "^I61", "^I62", "^I63", "^I64")

matches_codes <- function(code, patterns) {
  any(sapply(patterns, function(p) grepl(p, code)))
}

process_row <- function(row) {
  codes_41270 <- extract_codes(row[["41270-0.0"]])
 
  matches_41270 <- sapply(codes_41270, function(x) matches_codes(x, target_stroke))
 
  stroke_diag <- 0
  earliest_date_stroke <- NA

  if (any(matches_41270)) {
    stroke_diag <- 1

    indices <- which(matches_41270)

    all_dates <- c()
    
    for (index in indices) {
      date_col <- paste0("41280-0.", index-1) 
      if (date_col %in% names(row) && !is.na(row[[date_col]]) && row[[date_col]] != "") {
        if (grepl("^\\d{4}-\\d{2}-\\d{2}$", row[[date_col]])) {
          all_dates <- c(all_dates, row[[date_col]])
        }
      }
    }
    
    if (length(all_dates) > 0) {
      earliest_date_stroke <- min(as.Date(all_dates, format="%Y-%m-%d"), na.rm = TRUE)
    }
  }
  
  return(c(stroke = stroke_diag, date_stroke = as.character(earliest_date_stroke)))
}

resultado <- t(apply(hes, 1, process_row))
hes$stroke <- as.numeric(resultado[, "stroke"])
hes$date_stroke <- as.Date(resultado[, "date_stroke"])




########## Diabetes mellitus - type 1 #############
extract_codes <- function(x) {
  if (is.na(x) || x == "") return(character(0))
  unlist(strsplit(x, "\\|"))
}

target_t1dm <- c("^E10")

matches_codes <- function(code, patterns) {
  any(sapply(patterns, function(p) grepl(p, code)))
}

process_row <- function(row) {
  codes_41270 <- extract_codes(row[["41270-0.0"]])
  matches_41270 <- sapply(codes_41270, function(x) matches_codes(x, target_t1dm))

  t1dm_diag <- 0
  earliest_date_t1dm <- NA

  if (any(matches_41270)) {
    t1dm_diag <- 1

    indices <- which(matches_41270)

    all_dates <- c()
    
    for (index in indices) {
      date_col <- paste0("41280-0.", index-1) 
      if (date_col %in% names(row) && !is.na(row[[date_col]]) && row[[date_col]] != "") {
        if (grepl("^\\d{4}-\\d{2}-\\d{2}$", row[[date_col]])) {
          all_dates <- c(all_dates, row[[date_col]])
        }
      }
    }

    if (length(all_dates) > 0) {
      earliest_date_t1dm <- min(as.Date(all_dates, format="%Y-%m-%d"), na.rm = TRUE)
    }
  }
  
  return(c(t1dm = t1dm_diag, date_t1dm = as.character(earliest_date_t1dm)))
}

resultado <- t(apply(hes, 1, process_row))
hes$t1dm <- as.numeric(resultado[, "t1dm"])
hes$date_t1dm <- as.Date(resultado[, "date_t1dm"])

########## Diabetes mellitus - type 2 #############
extract_codes <- function(x) {
  if (is.na(x) || x == "") return(character(0))
  unlist(strsplit(x, "\\|"))
}

target_t2dm <- c("^E11")

matches_codes <- function(code, patterns) {
  any(sapply(patterns, function(p) grepl(p, code)))
}

process_row <- function(row) {
  codes_41270 <- extract_codes(row[["41270-0.0"]])
  matches_41270 <- sapply(codes_41270, function(x) matches_codes(x, target_t2dm))

  t2dm_diag <- 0
  earliest_date_t2dm <- NA

  if (any(matches_41270)) {
    t2dm_diag <- 1

    indices <- which(matches_41270)

    all_dates <- c()
    
    for (index in indices) {
      date_col <- paste0("41280-0.", index-1) 
      if (date_col %in% names(row) && !is.na(row[[date_col]]) && row[[date_col]] != "") {
        if (grepl("^\\d{4}-\\d{2}-\\d{2}$", row[[date_col]])) {
          all_dates <- c(all_dates, row[[date_col]])
        }
      }
    }
    
    if (length(all_dates) > 0) {
      earliest_date_t2dm <- min(as.Date(all_dates, format="%Y-%m-%d"), na.rm = TRUE)
    }
  }
  
  return(c(t2dm = t2dm_diag, date_t2dm = as.character(earliest_date_t2dm)))
}

resultado <- t(apply(hes, 1, process_row))
hes$t2dm <- as.numeric(resultado[, "t2dm"])
hes$date_t2dm <- as.Date(resultado[, "date_t2dm"])

########## Diabetes mellitus - Malnutrition-related #############
extract_codes <- function(x) {
  if (is.na(x) || x == "") return(character(0))
  unlist(strsplit(x, "\\|"))
}

target_mdm <- c("^E12")

matches_codes <- function(code, patterns) {
  any(sapply(patterns, function(p) grepl(p, code)))
}

process_row <- function(row) {
  codes_41270 <- extract_codes(row[["41270-0.0"]])
  matches_41270 <- sapply(codes_41270, function(x) matches_codes(x, target_mdm))
  
  mdm_diag <- 0
  earliest_date_mdm <- NA

  if (any(matches_41270)) {
    mdm_diag <- 1

    indices <- which(matches_41270)

    all_dates <- c()
    
    for (index in indices) {
      date_col <- paste0("41280-0.", index-1) 
      if (date_col %in% names(row) && !is.na(row[[date_col]]) && row[[date_col]] != "") {
        if (grepl("^\\d{4}-\\d{2}-\\d{2}$", row[[date_col]])) {
          all_dates <- c(all_dates, row[[date_col]])
        }
      }
    }

    if (length(all_dates) > 0) {
      earliest_date_mdm <- min(as.Date(all_dates, format="%Y-%m-%d"), na.rm = TRUE)
    }
  }
  
  return(c(mdm = mdm_diag, date_mdm = as.character(earliest_date_mdm)))
}

resultado <- t(apply(hes, 1, process_row))
hes$mdm <- as.numeric(resultado[, "mdm"])
hes$date_mdm <- as.Date(resultado[, "date_mdm"])

########## Diabetes mellitus - Other specified #############
extract_codes <- function(x) {
  if (is.na(x) || x == "") return(character(0))
  unlist(strsplit(x, "\\|"))
}

target_osdm <- c("^E13")
matches_codes <- function(code, patterns) {
  any(sapply(patterns, function(p) grepl(p, code)))
}

process_row <- function(row) {
  codes_41270 <- extract_codes(row[["41270-0.0"]])
  matches_41270 <- sapply(codes_41270, function(x) matches_codes(x, target_osdm))

  osdm_diag <- 0
  earliest_date_osdm <- NA

  if (any(matches_41270)) {
    osdm_diag <- 1

    indices <- which(matches_41270)

    all_dates <- c()
    
    for (index in indices) {
      date_col <- paste0("41280-0.", index-1) 
      if (date_col %in% names(row) && !is.na(row[[date_col]]) && row[[date_col]] != "") {
        if (grepl("^\\d{4}-\\d{2}-\\d{2}$", row[[date_col]])) {
          all_dates <- c(all_dates, row[[date_col]])
        }
      }
    }

    if (length(all_dates) > 0) {
      earliest_date_osdm <- min(as.Date(all_dates, format="%Y-%m-%d"), na.rm = TRUE)
    }
  }
  
  return(c(osdm = osdm_diag, date_osdm = as.character(earliest_date_osdm)))
}

resultado <- t(apply(hes, 1, process_row))
hes$osdm <- as.numeric(resultado[, "osdm"])
hes$date_osdm <- as.Date(resultado[, "date_osdm"])

########## Diabetes mellitus - Unspecified #############
extract_codes <- function(x) {
  if (is.na(x) || x == "") return(character(0))
  unlist(strsplit(x, "\\|"))
}

target_usdm <- c("^E14")

matches_codes <- function(code, patterns) {
  any(sapply(patterns, function(p) grepl(p, code)))
}

process_row <- function(row) {
  codes_41270 <- extract_codes(row[["41270-0.0"]])
  matches_41270 <- sapply(codes_41270, function(x) matches_codes(x, target_usdm))
  
  usdm_diag <- 0
  earliest_date_usdm <- NA

  if (any(matches_41270)) {
    usdm_diag <- 1

    indices <- which(matches_41270)

    all_dates <- c()
    
    for (index in indices) {
      date_col <- paste0("41280-0.", index-1) 
      if (date_col %in% names(row) && !is.na(row[[date_col]]) && row[[date_col]] != "") {
        if (grepl("^\\d{4}-\\d{2}-\\d{2}$", row[[date_col]])) {
          all_dates <- c(all_dates, row[[date_col]])
        }
      }
    }

    if (length(all_dates) > 0) {
      earliest_date_usdm <- min(as.Date(all_dates, format="%Y-%m-%d"), na.rm = TRUE)
    }
  }
  
  return(c(usdm = usdm_diag, date_usdm = as.character(earliest_date_usdm)))
}

resultado <- t(apply(hes, 1, process_row))
hes$usdm <- as.numeric(resultado[, "usdm"])
hes$date_usdm <- as.Date(resultado[, "date_usdm"])

########## Overall CVD (chd, chf and stroke) #############

hes$cvd <- ifelse(hes$chd == 1 | 
                           hes$chf == 1 | 
                           hes$stroke == 1, 1, 0)

freq(hes$cvd)

hes$date_cvd <- as.Date(
  with(hes, 
       ifelse(cvd == 1,
              pmin(date_chd, date_chf, date_stroke, na.rm = TRUE),
              NA)
  ),
  origin = "1970-01-01"
)

freq(hes$date_cvd)

############ Overall DM (not included type 1 DM) ###############
hes$dm <- ifelse(hes$t2dm == 1 | 
                          hes$mdm == 1 | 
                          hes$osdm == 1 |
                          hes$usdm == 1, 1, 0)

freq(hes$dm)

hes$date_dm <- as.Date(
  with(hes, 
       ifelse(dm == 1,
              pmin(date_t2dm, date_mdm, date_osdm, date_usdm, na.rm = TRUE),
              NA)
  ),
  origin = "1970-01-01"
)

freq(hes$date_dm)

########### Death - CVD ###########

target_death_cvd <- c(
  "^I00", "^I01", "^I02", "^I05", "^I06", "^I07", "^I08",
  "^I09", "^I11", "^I10",
  "^I20", "^I21", "^I22", "^I23", "^I24", "^I25", 
  "^I27", "^I28", "^I30", "^I31", "^I32", "^I33", "^I34", "^I35",
  "^I36", "^I37", "^I38", "^I39", "^I40", "^I41", "^I42", "^I43",
  "^I44", "^I45", "^I46", "^I47", "^I48", "^I49",
  "^I50", "^I51", "^I52", 
  "^I12", "^I13", "^I14", "^I15",
  "^I60", "^I61", "^I62", "^I63", "^I64", "^I65", "^I66", "^I67",
  "^I68", "^I69", "^I70", "^I72", "^I73", "^I71", "^I74", "^I77",
  "^I78", "^I79", "^I80", "^I81", "^I82", "^I83", "^I84", "^I85",
  "^I86", "^I87","^I88", 
  "^I95", "^I96", "^I97", "^I98", "^I99"
  )


hes$death_cvd <- 0

#looking for the targeted codes in the variable 40001-0.0 (primary cause of death)
for (pattern in target_death_cvd) {
  matches_40001 <- grepl(pattern, hes$`40001-0.0`, fixed = FALSE)
  hes$death_cvd[!is.na(matches_40001) & matches_40001] <- 1
}

# looking for the targeted codes from variable 40002-0.1 to 40002-0.14 (secondary causes of death)
for (i in 1:14) {
  col_name <- paste0("40002-0.", i)
  
  if (col_name %in% names(hes)) {  
    for (pattern in target_death_cvd) {
      matches <- grepl(pattern, hes[[col_name]], fixed = FALSE)
      hes$death_cvd[!is.na(matches) & matches] <- 1
    }
  }
}


hes$date_death_cvd <- NA
is_death_cvd <- !is.na(hes$death_cvd) & hes$death_cvd == 1
hes$date_death_cvd[is_death_cvd] <- hes$`40000-0.0`[is_death_cvd]

hes$date_death_cvd <- as.Date(hes$date_death_cvd, origin = "1970-01-01")


########## ANY DEATH ###############
#Differently from Carmen's do file I built this variable if the person has a date of death without consider the cause,
#because 303 people have date of death but not the cause registered. 

hes$any_death  <- ifelse(is.na(hes$`40000-0.0`), 0, 1)

hes$any_death_date <- ifelse(hes$any_death==1, hes$`40000-0.0`, NA)
hes$any_death_date <- as.Date(hes$any_death_date, origin = "1970-01-01")
freq(hes$any_death_date)

#44500 deaths

######## LABELS ########

attr(hes$chd, "label") <- "Diagnosis of Coronary Heart Disease"
attr(hes$date_chd, "label") <- "Date - Diagnosis of Coronary Heart Disease"

attr(hes$chf, "label") <- "Diagnosis of Congestive Heart Failure"
attr(hes$date_chf, "label") <- "Date - Diagnosis of Congestive Heart Failure"

attr(hes$stroke, "label") <- "Diagnosis of Stroke"
attr(hes$date_stroke, "label") <- "Date - Diagnosis of Stroke"

attr(hes$cvd, "label") <- "Diagnosis of Cardiovascular Disease"
attr(hes$date_cvd, "label") <- "Date - Diagnosis of Cardiovascular Disease"

attr(hes$t1dm, "label") <- "Diagnosis of Type 1 Diabetes Mellitus"
attr(hes$date_t1dm, "label") <- "Date - Diagnosis of Type 1 Diabetes Mellitus"

attr(hes$t2dm, "label") <- "Diagnosis of Type 2 Diabetes Mellitus"
attr(hes$date_t2dm, "label") <- "Date - Diagnosis of Type 2 Diabetes Mellitus"

attr(hes$mdm, "label") <- "Diagnosis of Malnutrition-related Diabetes Mellitus"
attr(hes$date_mdm, "label") <- "Date - Diagnosis of Malnutrition-related Diabetes Mellitus"

attr(hes$osdm, "label") <- "Diagnosis of Other specified Diabetes Mellitus"
attr(hes$date_osdm, "label") <- "Date - Diagnosis of Other specified Diabetes Mellitus"

attr(hes$usdm, "label") <- "Diagnosis of Unspecified Diabetes Mellitus"
attr(hes$date_usdm, "label") <- "Date - Diagnosis of Unspecified Diabetes Mellitus"

attr(hes$dm, "label") <- "Diagnosis of Overall Diabetes Mellitus (except t1dm)"
attr(hes$date_dm, "label") <- "Date - Diagnosis of Overall Diabetes Mellitus (except t1dm)"

attr(hes$death_cvd, "label") <- "Death by CVD"
attr(hes$date_death_cvd, "label") <- "Date - Death by CVD"

attr(hes$any_death, "label") <- "Death by any cause"
attr(hes$any_death_date, "label") <- "Date - Death by any cause"


####### frequency tables for the variables ######
freq(hes$death_cvd)
freq(hes$any_death)
freq(hes$chd)
freq(hes$chf)
freq(hes$stroke)
freq(hes$cvd)
freq(hes$t1dm)
freq(hes$t2dm)
freq(hes$mdm)
freq(hes$osdm)
freq(hes$usdm)
freq(hes$dm)

freq(hes$any_death_date)
freq(hes$date_death_cvd)
freq(hes$date_chd)
freq(hes$date_chf)
freq(hes$date_stroke)
freq(hes$date_t1dm)
freq(hes$date_t2dm)
freq(hes$date_mdm)
freq(hes$date_osdm)
freq(hes$date_usdm)

# Save a file with the HES results -----
# saveRDS(hes, "/mnt/project/data/processed/HES.RDS")


# Before or after of the date of attending asessment centre ------
# Label:
# 1- it had diagnosis AFTER enrollment in the study
# 2- it had diagnosis BEFORE or AT THE SAME DAY of enrollment in the study
# NA - do not have a diagnosis

### to identify cases where the date of diagnosis is higher (not equal) to Date of attending assessment centre for the first time (variable 53-0.0)
### the variable 53-0.0 is in the dataset DEMO, so it will be necessary to merge DEMO and HES dataset. The codes for this are below:

hes = readRDS("/mnt/project/data/processed/HES.RDS")

#Creating a new dataset with only the variables of interest (to run quicker)
hes <- subset(hes, select= c("eid", "chd", "date_chd", "chf", "date_chf", "stroke", "date_stroke", "cvd", "date_cvd", "t1dm", "date_t1dm", 
                        "t2dm", "date_t2dm", "mdm", "date_mdm", "osdm", "date_osdm", "usdm", "date_usdm", "dm", "date_dm", "death_cvd", 
                        "date_death_cvd", "any_death", "any_death_date", "31-0.0")
)

# EDIT ACCORDING TO YOUR FILE with the variable 53-0.0 ------
demo <- read_csv("/mnt/project/data/raw/Demographics.csv")  # to read demo

#creating a subset only with the variable 53-0.0 to be quicker
demo <- subset(demo, select= c("eid", "53-0.0"))

#merging datasets
merged_hes <- merge(hes, demo, by="eid", all = TRUE)

#Label:
  #1- it had diagnosis AFTER enrollment in the study
  #2- it had diagnosis BEFORE or AT THE SAME DAY of enrollment in the study
  #NA - do not have a diagnosis
  
# for chd
merged_hes$chd_inc <- NA  
merged_hes$chd_inc <- ifelse(
  merged_hes$chd == 1 & merged_hes$date_chd > merged_hes$`53-0.0`, 1,
  ifelse(merged_hes$chd == 1 & merged_hes$date_chd <= merged_hes$`53-0.0`, 2, NA)
)
freq(merged_hes$chd_inc)

# for chf
merged_hes$chf_inc <- NA  
merged_hes$chf_inc <- ifelse(
  merged_hes$chf == 1 & merged_hes$date_chf > merged_hes$`53-0.0`, 1,
  ifelse(merged_hes$chf == 1 & merged_hes$date_chf <= merged_hes$`53-0.0`, 2, NA)
)
freq(merged_hes$chf_inc)

# for stroke
merged_hes$stroke_inc <- NA  
merged_hes$stroke_inc <- ifelse(
  merged_hes$stroke == 1 & merged_hes$date_stroke > merged_hes$`53-0.0`, 1,
  ifelse(merged_hes$stroke == 1 & merged_hes$date_stroke <= merged_hes$`53-0.0`, 2, NA)
)
freq(merged_hes$stroke_inc)

# for t1dm
merged_hes$t1dm_inc <- NA  
merged_hes$t1dm_inc <- ifelse(
  merged_hes$t1dm == 1 & merged_hes$date_t1dm > merged_hes$`53-0.0`, 1,
  ifelse(merged_hes$t1dm == 1 & merged_hes$date_t1dm <= merged_hes$`53-0.0`, 2, NA)
)
freq(merged_hes$t1dm_inc)

# for t2dm
merged_hes$t2dm_inc <- NA  
merged_hes$t2dm_inc <- ifelse(
  merged_hes$t2dm == 1 & merged_hes$date_t2dm > merged_hes$`53-0.0`, 1,
  ifelse(merged_hes$t2dm == 1 & merged_hes$date_t2dm <= merged_hes$`53-0.0`, 2, NA)
)
freq(merged_hes$t2dm_inc)

# for mdm
merged_hes$mdm_inc <- NA  
merged_hes$mdm_inc <- ifelse(
  merged_hes$mdm == 1 & merged_hes$date_mdm > merged_hes$`53-0.0`, 1,
  ifelse(merged_hes$mdm == 1 & merged_hes$date_mdm <= merged_hes$`53-0.0`, 2, NA)
)
freq(merged_hes$mdm_inc)

# for osdm
merged_hes$osdm_inc <- NA  
merged_hes$osdm_inc <- ifelse(
  merged_hes$osdm == 1 & merged_hes$date_osdm > merged_hes$`53-0.0`, 1,
  ifelse(merged_hes$osdm == 1 & merged_hes$date_osdm <= merged_hes$`53-0.0`, 2, NA)
)
freq(merged_hes$osdm_inc)

# for usdm
merged_hes$usdm_inc <- NA  
merged_hes$usdm_inc <- ifelse(
  merged_hes$usdm == 1 & merged_hes$date_usdm > merged_hes$`53-0.0`, 1,
  ifelse(merged_hes$usdm == 1 & merged_hes$date_usdm <= merged_hes$`53-0.0`, 2, NA)
)
freq(merged_hes$usdm_inc)

# saveRDS(merged_hes, "/mnt/project/data/processed/2025_demographics_with_hesoutcome.RDS")
# we renamed as HES.RDS
  
  
  
  
  
  
  
  
  
  
  
  
  
  
