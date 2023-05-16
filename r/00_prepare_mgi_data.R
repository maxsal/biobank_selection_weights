# Prepare MGI data including deriving variables for comorbidity status and
# descriptive variables derived from demographic data
# author:   max salvatore
# date:     20230418

# libraries, paths, and such ---------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(glue)
  library(qs)
  library(parallel)
  library(optparse)
})

# optparse list ----
option_list <- list(
  make_option("--mgi_version",
    type = "character", default = "20220822",
    help = "Cohort version in /net/junglebook/magic_data/EHRdata/ [default = %default]"
  )
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

message(glue("using cohort version {opt$mgi_version}; see /net/junglebook/magic_data/EHRdata/"))

out_path <- glue("/net/junglebook/home/mmsalva/projects/dissertation/aim_one/data/private/mgi/{opt$mgi_version}/")
if (!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)

source("fn/files-utils.R")
file_paths <- get_files(mgi_version = opt$mgi_version)

# load data --------------------------------------------------------------------
message("loading data...")
study <- fread("/net/junglebook/magic_data/Data_Pulls_from_Data_Office/MGI_Study_FirstEnrollment_20221102.txt")
MGIcohort <- fread(file_paths[["mgi"]][["cov_file"]])
load(file = file_paths[["mgi"]][["phe_overview_file"]])
load(file = file_paths[["mgi"]][["phecode_dsb_file"]])
cancer_phecodes <- fread("/net/junglebook/home/mmsalva/projects/dissertation/aim_one/data/public/cancer_phecodes.txt",
  colClasses = "character"
)[[1]]

### replace FirstDaySinceBirth and LastDaySinceBirth
# following a Slack conversation with Lars on 2/21/23, the original values for
# these variables include ICD codes that cannot be mapped to a phecode and
# non-ICD code records (e.g., OrderDate_DaysSinceBirth from LabResults)
MGIcohort <- MGIcohort[, !c("FirstDaySinceBirth", "LastDaySinceBirth")]

first_dsb <- unique(diagnoses_Phecodes[diagnoses_Phecodes[
  ,
  .I[which.min(DaysSinceBirth)],
  "IID"
][["V1"]]][
  ,
  .(
    DeID_PatientID         = IID,
    FirstDaySinceBirth     = DaysSinceBirth,
    age_at_first_diagnosis = round(DaysSinceBirth / 365.25, 1)
  )
])
last_dsb <- unique(diagnoses_Phecodes[diagnoses_Phecodes[, .I[which.max(DaysSinceBirth)], "IID"][["V1"]]][, .(
  DeID_PatientID        = IID,
  LastDaySinceBirth     = DaysSinceBirth,
  age_at_last_diagnosis = round(DaysSinceBirth / 365.25, 1)
)])

MGIcohort <- Reduce(
  \(x, y) merge.data.table(x, y, by = "DeID_PatientID"),
  list(MGIcohort, first_dsb, last_dsb)
)
###

MGIcohort <- merge.data.table(
  MGIcohort,
  study,
  by = "DeID_PatientID"
)[StudyName != "AOS", ]

comorbid <- list(
  "cad"                = list("phecodes" = c("411.4")),
  "diabetes"           = list("phecodes" = c("250")),
  "hypertension"       = list("phecodes" = c("272.12")),
  "mixed_hypertension" = list("phecodes" = c("272.13")),
  "vitamin_d"          = list("phecodes" = c("261.4")),
  "depression"         = list("phecodes" = c("296.2")),
  "anxiety"            = list("phecodes" = c("300")),
  "bipolar"            = list("phecodes" = c("296.1")),
  "cancer"             = list("phecodes" = cancer_phecodes)
)

# identify cases and create indicator variables --------------------------------
message("identifying cases and creating indicator variables...")
pb <- txtProgressBar(max = length(names(comorbid)), width = 50, style = 3)
for (i in names(comorbid)) {
  comorbid[[i]][["ids"]] <- diagnoses_Phecodes[phecode %in% comorbid[[i]][["phecodes"]], IID] |>
    unique()
  setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
}
close(pb)

pb <- txtProgressBar(max = length(names(comorbid)), width = 50, style = 3)
for (i in names(comorbid)) {
  set(MGIcohort, j = i, value = fifelse(MGIcohort[["DeID_PatientID"]] %in% comorbid[[i]][["ids"]], 1, 0))
  setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
}
close(pb)

MGIcohort[, triglycerides := fifelse(hypertension == 0 & mixed_hypertension == 0, 0, 1)]

MGIcohort[, nhanes_nhw := fifelse(Ethnicity != "Hispanic" & Race == "Caucasian", 1, 0)]

MGIcohort[, age_cat := between(AgeLastEntry, 0, 5.99) +
  2 * between(AgeLastEntry, 6, 11.99) +
  3 * between(AgeLastEntry, 12, 19.99) +
  4 * between(AgeLastEntry, 20, 39.99) +
  5 * between(AgeLastEntry, 40, 59.99) +
  6 * between(AgeLastEntry, 60, 150.99)]

setnames(MGIcohort, "BMI", "bmi")
MGIcohort[
  , bmi_cat := fcase(
    between(bmi, 0, 18.499), 1, # underweight
    between(bmi, 18.5, 24.999), 2, # "normal"
    between(bmi, 25.0, 29.999), 3, # overweight
    between(bmi, 30, 120), 4
  ) # obese
][, `:=`(
  bmi_under       = as.numeric(bmi_cat == 1),
  bmi_overweight  = as.numeric(bmi_cat == 3),
  bmi_obese       = as.numeric(bmi_cat == 4),
  smoking_current = as.numeric(SmokingStatus == "Current"),
  smoking_former  = as.numeric(SmokingStatus == "Past")
)]

MGIcohort[, female := as.numeric(Sex == "F")]

# saving files -----------------------------------------------------------------
message("saving processed files...")

save_qs(MGIcohort, file = glue("{out_path}data_{opt$mgi_version}_comb.qs"))
for (i in c("MGI", "MGI-MEND", "MIPACT", "MHB2")) {
  save_qs(MGIcohort[StudyName == i, ], file = glue("{out_path}data_{opt$mgi_version}_{tolower(i)}.qs"))
}

message(paste0("script success! see ", out_path, " for output files"))
