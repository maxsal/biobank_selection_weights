# libraries, paths, and such ---------------------------------------------------
library(data.table)
library(glue)
library(cli)
library(optparse)

# optparse list ---
option_list <- list(
  make_option("--cohort_version", type = "character", default = "20220822",
              help = "Cohort version in /net/junglebook/magic_data/EHRdata/ [default = '20220822']")
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args   <- parse_args(parser, positional_arguments = 0)
opt    <- args$options

cli_alert_info("using cohort version {opt$cohort_version}; see {.path /net/junglebook/magic_data/EHRdata/}")

data_path <- glue("/net/junglebook/magic_data/EHRdata/{opt$cohort_version}/")
out_path  <- glue("/net/junglebook/home/mmsalva/projects/dissertation/aim_one/data/private/mgi/{opt$cohort_version}/")
dir.create(out_path, recursive = TRUE)

# load data --------------------------------------------------------------------
cli_alert("loading data...")
study     <- fread("/net/junglebook/magic_data/Data_Pulls_from_Data_Office/MGI_Study_FirstEnrollment_20221102.txt")
load(file = glue("{data_path}MGI_20220822.Rsav"))
load(file = glue("{data_path}phenomes/UNFILTERED_20220822/UNFILTERED_20220822_Phenotype_Overview_All_Phecodes1plus.Rsav"))
load(file = glue("{data_path}phenomes/UNFILTERED_20220822/UNFILTERED_20220822_Phecodes_Birthyears.Rsav"))
cancer_phecodes <- fread("/net/junglebook/home/mmsalva/projects/dissertation/aim_one/data/public/cancer_phecodes.txt",
                         colClasses = "character")[[1]]

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
cli_alert("identifying cases and creating indicator variables...")
cli_progress_bar("comorbid ids", total = length(names(comorbid)))
for (i in names(comorbid)) {
  comorbid[[i]][["ids"]] <- diagnoses_Phecodes[phecode %in% comorbid[[i]][["phecodes"]], IID] |>
    unique()
  cli_progress_update()
}

cli_progress_bar("comorbid vars", total = length(names(comorbid)))
for (i in names(comorbid)) {
  MGIcohort[[i]] <- fifelse(MGIcohort[["DeID_PatientID"]] %in% comorbid[[i]][["ids"]], 1, 0)
  cli_progress_update()
}

MGIcohort[, Triglycerides := fifelse(Hypertension == 0 & `Mixed Hypertension` == 0, 0, 1)]

MGIcohort[, nhanes_nhw := fifelse(Ethnicity != "Hispanic" & Race == "Caucasian", 1, 0)]

# saving files -----------------------------------------------------------------
cli_alert("saving processed files...")
saveRDS(MGIcohort[StudyName == "MGI", ], file = glue("{out_path}data{opt$cohort_version}_bb.rds"))
saveRDS(MGIcohort[StudyName == "MHB2", ], file = glue("{out_path}data{opt$cohort_version}_mhb.rds"))
saveRDS(MGIcohort[StudyName == "MIPACT", ], file = glue("{out_path}data{opt$cohort_version}_mipact.rds"))
saveRDS(MGIcohort[StudyName == "MGI-MEND", ], file = glue("{out_path}data{opt$cohort_version}_mend.rds"))

cli_alert_success("script success! see {.path {out_path}} for output files")
