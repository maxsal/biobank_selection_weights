# clean nhanes data
# requires: demographic, alcohol, and smoking questionnaire XPT files in
#           '/data/public/nhanes'
# outputs:  cleaned nhanes data
# author:   max salvatore
# date:     20221215

# libraries, functions, and options --------------------------------------------
library(data.table)
library(haven)
library(glue)
library(progress)
library(cli)

source("fn/cleaning-utils.R")

# specifications ---------------------------------------------------------------
nhanes_data_path   <- "data/public/nhanes/"
nhanes_data_prefix <- "P"
output_path        <- "data/public/nhanes/"
save               <- TRUE

# demographics -----------------------------------------------------------------
cli_alert_info("processing demographics...")
demo <- read_xpt(glue("{nhanes_data_path}{nhanes_data_prefix}_DEMO.XPT")) |>
  as.data.table()
demo_vars <- c("SEQN", "RIAGENDR", "RIDAGEYR", "RIDRETH1", "RIDRETH3",
               "DMDEDUC2", "DMDMARTZ", "WTINTPRP", "WTMECPRP",
               "SDMVPSU", "SDMVSTRA")
demo <- demo[, ..demo_vars]
setnames(demo,
         demo_vars,
         c("id", "gender", "age", "race_eth1", "race_eth3", "education",
           "marital_status", "interview_weight", "mec_weight", "psu", "strata"))
demo[, names(demo) := lapply(.SD, as.character)]
num_cols <- c("age", "interview_weight", "mec_weight")
demo[, (num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols]
demo <- demo[age >= 18, ]
demo[, gender := fcase(
  gender == "1", "Male",
  gender == "2", "Female"
  )]
demo[, race_eth1 := fcase(
  race_eth1 == "1", "Mexican American",
  race_eth1 == "2", "Other Hispanic",
  race_eth1 == "3", "Non-Hispanic White",
  race_eth1 == "4", "Non-Hispanic Black",
  race_eth1 == "5", "Other Race - Including Multi-Racial",
  default = "Unknown"
)]
demo[, race_eth3 := fcase(
  race_eth3 == "1", "Mexican American",
  race_eth3 == "2", "Other Hispanic",
  race_eth3 == "3", "Non-Hispanic White",
  race_eth3 == "4", "Non-Hispanic Black",
  race_eth3 == "6", "Non-Hipsanic Asian",
  race_eth3 == "7", "Other Race - Including Multi-Racial",
  default = "Unknown"
)]
demo[, nhw := as.numeric(race_eth1 == "Non-Hispanic White")]
demo[, education := fcase(
  education == "1", "Less than 9th grade",
  education == "2", "9-11th grade (includes 12th grade with no diploma)",
  education == "3", "High school graduate/GED or equivalent",
  education == "4", "Some college of AA degree",
  education == "5", "College graduate or above",
  default = "Unknown"
)]
demo[, marital_status := fcase(
  marital_status == "1", "Married/Living with Partner",
  marital_status == "2", "Widowed/Divorced/Separated",
  marital_status == "3", "Never married",
  default = "Unknown"
)]
demo[, female := as.numeric(gender == "Female")]

# alcohol ----------------------------------------------------------------------
cli_alert_info("processing alcohol...")
alc <- read_xpt(glue("{nhanes_data_path}{nhanes_data_prefix}_ALQ.XPT")) |>
  as.data.table()
alc_vars <- c("SEQN", "ALQ111", "ALQ121")
alc <- alc[, ..alc_vars]
setnames(alc,
         alc_vars,
         c("id", "drink_ever", "drink_last_year"))
alc[, names(alc) := lapply(.SD, as.character)]
alc[, drink_ever := fcase(
  drink_ever == "1", "Yes",
  drink_ever == "2", "No",
  default = "Unknown"
)]
alc[, drink_last_year := fcase(
  drink_last_year == "0", "No",
  drink_last_year %in% as.character(1:10), "At least once",
  default = "Unknown"
)]
alc[, drinker := fcase(
  drink_ever == "No", "Never",
  drink_ever == "Yes" & drink_last_year == "No", "Former",
  drink_last_year == "At least once", "Current",
  default = "Unknown"
)]

# smoking ----------------------------------------------------------------------
cli_alert_info("processing smoking...")
smk <- read_xpt(glue("{nhanes_data_path}{nhanes_data_prefix}_SMQ.XPT")) |>
  as.data.table()
smk_vars <- c("SEQN", "SMQ020", "SMQ040")
smk <- smk[, ..smk_vars]
setnames(smk,
         smk_vars,
         c("id", "smoke_100", "smoke_now"))
smk[, names(smk) := lapply(.SD, as.character)]
smk[, smoke_100 := fcase(
  smoke_100 == "1", "Yes",
  smoke_100 == "2", "No",
  default = "Unknown"
)]
smk[, smoke_now := fcase(
  smoke_now %in% c("1", "2"), "Yes",
  smoke_now == "3", "No",
  default = "Unknown"
)]
smk[, smoker := fcase(
  smoke_100 == "No", "Never",
  smoke_100 == "Yes" & smoke_now == "No", "Former",
  smoke_100 == "Yes" & smoke_now == "Yes", "Current",
  default = "Unknown"
)]

# medical conditions -----------------------------------------------------------
cli_alert_info("processing medical conditions...")
mcq <- read_xpt(glue("{nhanes_data_path}{nhanes_data_prefix}_MCQ.XPT")) |>
  as.data.table()
mcq_vars <- c("SEQN", "MCQ010", "MCQ080", "MCQ160E", "MCQ160F", "MCQ220",
              "MCQ230A")
mcq <- mcq[, ..mcq_vars]
setnames(mcq,
         mcq_vars,
         c("id", "asthma_ever", "overweight_ever", "heart_attack", "stroke",
           "cancer_ever", "cancer_first"))
yes_no_recode <- function(x) {
  fcase(
    x == 1, "Yes",
    x == 2, "No",
    default = "Unknown"
  )
}
mcq[, c("asthma_ever", "overweight_ever", "heart_attack", "stroke",
        "cancer_ever") := lapply(.SD, yes_no_recode),
    .SDcols = c("asthma_ever", "overweight_ever", "heart_attack", "stroke",
                "cancer_ever")]
mcq[, id := as.character(id)]
mcq[, cancer_first := fcase(
  cancer_first == 10, "Bladder",
  cancer_first == 11, "Blood",
  cancer_first == 12, "Bone",
  cancer_first == 13, "Brain",
  cancer_first == 14, "Breast",
  cancer_first == 15, "Cervix (cervical)",
  cancer_first == 16, "Colon",
  cancer_first == 17, "Esophagus (esophageal)",
  cancer_first == 18, "Gallbladder",
  cancer_first == 19, "Kidney",
  cancer_first == 20, "Larynx/ windpipe",
  cancer_first == 21, "Leukemia",
  cancer_first == 22, "Liver",
  cancer_first == 23, "Lung",
  cancer_first == 24, "Lymphoma/ Hodgkin's disease",
  cancer_first == 25, "Melanoma",
  cancer_first == 26, "Mouth/tongue/lip",
  cancer_first == 27, "Nervous system",
  cancer_first == 28, "Ovary (ovarian)",
  cancer_first == 29, "Pancreas (pancreatic)",
  cancer_first == 30, "Prostate",
  cancer_first == 31, "Rectum (recal)",
  cancer_first == 32, "Skin (non-melanoma)",
  cancer_first == 33, "Skin (don't know what kind)",
  cancer_first == 34, "Soft tissue (muscle or fat)",
  cancer_first == 35, "Stomach",
  cancer_first == 36, "Testis (testicular)",
  cancer_first == 37, "Thyroid",
  cancer_first == 38, "Uterus (uterine)",
  cancer_first %in% c(39, 66, 77, 99), "Other",
  default = "None reported"
)]


# merge ------------------------------------------------------------------------
cli_alert_info("merging...")
out <- Reduce(\(x, y) {
    merge.data.table(x = x, y = y, by = "id", all.x = TRUE)
  },
  list(demo, alc, smk, mcq))
replace_missing(data = out,
                cols = c("race_eth1", "race_eth3", "education",
                         "marital_status", "drink_ever", "drink_last_year",
                         "drinker", "smoke_100", "smoke_now", "smoker"),
                new_value = "Unknown")
if (save == TRUE) {
  fwrite(out,
         file = glue("{output_path}{nhanes_data_prefix}",
                     "_NHANES_CLEAN.txt"),
         sep  = "\t")
}

cli_alert_success(glue(
  "nhanes data cleaned and saved: {output_path}{nhanes_data_prefix}",
  "_NHANES_CLEAN.txt"))
