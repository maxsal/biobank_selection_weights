# clean nhanes data
# requires: demographic, alcohol, and smoking questionnaire XPT files in
#           '/data/public/nhanes'
# outputs:  cleaned nhanes data
# author:   max salvatore
# date:     20221207

# libraries, functions, and options --------------------------------------------
library(data.table)
library(haven)
library(glue)
library(progress)

# specifications ---------------------------------------------------------------
nhanes_data_path   <- "data/public/nhanes/"
nhanes_data_prefix <- "P"
output_path        <- "data/public/nhanes/"
save               <- TRUE
  
# demographics -----------------------------------------------------------------
demo <- haven::read_xpt(glue::glue("{nhanes_data_path}{nhanes_data_prefix}_DEMO.XPT")) |>
  data.table::as.data.table()
demo_vars <- c("SEQN", "RIAGENDR", "RIDAGEYR", "RIDRETH1", "RIDRETH3",
               "DMDEDUC2", "DMDMARTZ", "WTINTPRP", "WTMECPRP",
               "SDMVPSU", "SDMVSTRA")
demo <- demo[, ..demo_vars]
data.table::setnames(demo,
                     demo_vars,
                     c("id", "gender", "age", "race_eth1", "race_eth3",
                       "education", "marital_status", "interview_weight",
                       "mec_weight", "psu", "strata"))
demo[, names(demo) := lapply(.SD, as.character)]
num_cols <- c("age", "interview_weight", "mec_weight")
demo[, (num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols]
demo <- demo[age >= 18, ]
demo[, gender := data.table::fcase(gender == "1", "Male", gender == "2", "Female")]
demo[, race_eth1 := data.table::fcase(
  race_eth1 == "1", "Mexican American",
  race_eth1 == "2", "Other Hispanic",
  race_eth1 == "3", "Non-Hispanic White",
  race_eth1 == "4", "Non-Hispanic Black",
  race_eth1 == "5", "Other Race - Including Multi-Racial",
  default = "Unknown"
)]
demo[, race_eth3 := data.table::fcase(
  race_eth3 == "1", "Mexican American",
  race_eth3 == "2", "Other Hispanic",
  race_eth3 == "3", "Non-Hispanic White",
  race_eth3 == "4", "Non-Hispanic Black",
  race_eth3 == "6", "Non-Hipsanic Asian",
  race_eth3 == "7", "Other Race - Including Multi-Racial",
  default = "Unknown"
)]
demo[, education := data.table::fcase(
  education == "1", "Less than 9th grade",
  education == "2", "9-11th grade (includes 12th grade with no diploma)",
  education == "3", "High school graduate/GED or equivalent",
  education == "4", "Some college of AA degree",
  education == "5", "College graduate or above",
  default = "Unknown"
)]
demo[, marital_status := data.table::fcase(
  marital_status == "1", "Married/Living with Partner",
  marital_status == "2", "Widowed/Divorced/Separated",
  marital_status == "3", "Never married",
  default = "Unknown"
)]
demo[, female := as.numeric(gender == "Female")]

# alcohol ----------------------------------------------------------------------
alc <- haven::read_xpt(glue::glue("{nhanes_data_path}{nhanes_data_prefix}_ALQ.XPT")) |>
  data.table::as.data.table()
alc_vars <- c("SEQN", "ALQ111", "ALQ121")
alc <- alc[, ..alc_vars]
data.table::setnames(alc,
                     alc_vars,
                     c("id", "drink_ever", "drink_last_year"))
alc[, names(alc) := lapply(.SD, as.character)]
alc[, drink_ever := data.table::fcase(
  drink_ever == "1", "Yes",
  drink_ever == "2", "No",
  default = "Unknown"
)]
alc[, drink_last_year := data.table::fcase(
  drink_last_year == "0", "No",
  drink_last_year %in% as.character(1:10), "At least once",
  default = "Unknown"
)]
alc[, drinker := data.table::fcase(
  drink_ever == "No", "Never",
  drink_ever == "Yes" & drink_last_year == "No", "Former",
  drink_last_year == "At least once", "Current",
  default = "Unknown"
)]

# smoking ----------------------------------------------------------------------
smk <- haven::read_xpt(glue::glue("{nhanes_data_path}{nhanes_data_prefix}_SMQ.XPT")) |>
  data.table::as.data.table()
smk_vars <- c("SEQN", "SMQ020", "SMQ040")
smk <- smk[, ..smk_vars]
data.table::setnames(smk,
                     smk_vars,
                     c("id", "smoke_100", "smoke_now"))
smk[, names(smk) := lapply(.SD, as.character)]
smk[, smoke_100 := data.table::fcase(
  smoke_100 == "1", "Yes",
  smoke_100 == "2", "No",
  default = "Unknown"
)]
smk[, smoke_now := data.table::fcase(
  smoke_now %in% c("1", "2"), "Yes",
  smoke_now == "3", "No",
  default = "Unknown"
)]
smk[, smoker := data.table::fcase(
  smoke_100 == "No", "Never",
  smoke_100 == "Yes" & smoke_now == "No", "Former",
  smoke_100 == "Yes" & smoke_now == "Yes", "Current",
  default = "Unknown"
)]

# merge ------------------------------------------------------------------------
out <- Reduce(\(x, y) data.table::merge.data.table(x = x, y = y, by = "id", all.x = TRUE),
              list(demo, alc, smk))
replace_missing(data = out,
                cols = c("race_eth1", "race_eth3", "education",
                         "marital_status", "drink_ever", "drink_last_year",
                         "drinker", "smoke_100", "smoke_now", "smoker"),
                new_value = "Unknown")

if (save == TRUE) {
  data.table::fwrite(out,
                     file = glue::glue("{output_path}{nhanes_data_prefix}_NHANES_CLEAN.txt"),
                     sep  = "\t")
}
