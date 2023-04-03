# quickly perform phecode-phecode phewas using MGI data for multiple
# time-thresholds for a selected outcome variable
# requires: time-threshold phecode indicator matrices must already exist
# outputs:  betas, sebetas, and p-values
# author:   max salvatore
# date:     20230220

# libraries, functions, and options --------------------------------------------
suppressPackageStartupMessages({
    library(data.table)
    library(MatchIt)
    library(logistf)
    library(glue)
    library(qs)
    library(cli)
    library(optparse)
    library(purrr)
})

set.seed(61787)

lapply(list.files("fn/", full.names = TRUE), source) |> # load functions
  invisible()

# optparse list ----------------------------------------------------------------
option_list <- list(
  make_option("--mgi_version", type = "character", default = "20220822",
              help = "Version of MGI data [default = %default]"),
  make_option("--mgi_cohort", type = "character", default = "comb",
              help = "Cohort of MGI used in weighting (comb, bb, mend, mhb) [default = %default]"),
  make_option("--ukb_version", type = "character", default = "20221117",
              help = "Version of UKB data [default = %default]")
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args   <- parse_args(parser, positional_arguments = 0)
opt    <- args$options
print(opt)

time_thresholds <- as.numeric(strsplit(opt$time_thresholds, ",")[[1]])

## extract file paths
file_paths <- get_files(mgi_version = opt$mgi_version,
                        ukb_version = opt$ukb_version)

# read data --------------------------------------------------------------------
## mgi
### demographics
mgi_cov <- read_qs(glue("data/private/mgi/{opt$mgi_version}/data_{opt$mgi_version}_comb.qs"))
setnames(mgi_cov,
         old = c("DeID_PatientID", "Age", "AliveYN", "Deceased_DaysSinceBirth",
                 "Ethnicity", "MaritalStatusCode", "Sex", "Race", "AgeFirstEntry",
                 "AgeLastEntry", "YearsInEHR", "FirstDaySinceBirth", "LastDaySinceBirth"),
         new = c("id", "age", "alive", "dead_dsb", "ethn", "marital",
                 "sex", "race", "age_at_first_diagnosis", "age_at_last_diagnosis",
                 "length_followup", "first_dsb", "last_dsb"))

### icd-phecode data
mgi_full_phe <- get(load(file_paths[["mgi"]][["phecode_dsb_file"]]))
if ("IID" %in% names(mgi_full_phe)) { setnames(mgi_full_phe, "IID", "id") }
if ("DaysSinceBirth" %in% names(mgi_full_phe)) { setnames(mgi_full_phe, "DaysSinceBirth", "dsb") }

## ukb
cli_alert("loading ukb data...")
### demographics
ukb_demo <- fread(file_paths[["ukb"]][["demo_file"]],
                              na.strings = c("", "NA", "."),
                              colClass = "character")
ukb_demo <- ukb_demo[, .(
  id   = as.character(id),
  dob  = as.Date(dob),
  age  = floor(as.numeric(as.Date("2022-11-17") - as.Date(dob)) / 365.25),
  ethn = ethnicity,
  sex)]
ukb_demo <- ukb_demo[complete.cases(ukb_demo), ]

### icd-phecode data
ukb_full_phe <- fread(file_paths[["ukb"]][["icd_phecode_file"]], colClasses = "character")
if ("IID" %in% names(ukb_full_phe)) { setnames(ukb_full_phe, "IID", "id") }
if ("DaysSinceBirth" %in% names(ukb_full_phe)) { setnames(ukb_full_phe, "DaysSinceBirth", "dsb") }
ukb_full_phe[, dsb := as.numeric(dsb)]

## other
cancer_phecodes <- fread("/net/junglebook/home/mmsalva/projects/dissertation/aim_one/data/public/cancer_phecodes.txt",
                         colClasses = "character")[[1]]

# prep data --------------------------------------------------------------------
## age categories
age_cats    <- seq(0, 80, 10)
tmp_age_cat <- age_grp_table(
    lower_ages   = age_cats,
    num_vec      = rep(NA, length(age_cats)),
    num_var_name = "counts"
)
## mgi
mgi_cov[, `:=` (
    race_eth      = factor(fifelse(ethn == "Hispanic", "Hispanic", fifelse(race == "", "Unknown", race))),
    SmokingStatus = relevel(factor(SmokingStatus, levels = c("Never", "Past", "Current", "Unknown")), ref = "Never"),
    bmi_verbose   = factor(fcase(
        bmi_cat == 1, "Underweight (<18.5)",
        bmi_cat == 2, "Healthy [18.5, 25)",
        bmi_cat == 3, "Overweight [25, 30)",
        bmi_cat == 4, "Obese [30+)"
    ), levels = c("Underweight (<18.5)", "Healthy [18.5, 25)", "Overweight [25, 30)", "Obese [30+)")),
    age_verbose = factor(cut(age_at_last_diagnosis,
                              breaks = c(0, tmp_age_cat[["upper"]]),
                              labels = tmp_age_cat[["group"]],
                              right  = FALSE),
                          levels = tmp_age_cat[["group"]])
)]
## ukb
ukb_demo[, `:=` (
    race_eth = fcase(
        ethn %in% c("Asian", "Chinese"), "Asian",
        ethn %in% c("Mixed"), "More than one population",
        ethn %in% c("Do not know", "Prefer not to answer"), "Unknown",
        ethn %in% c("Other ethnic group"), "Other",
        ethn == "White", "White",
        ethn == "Black", "Black"),
    age_verbose = factor(cut(age,
                            breaks = c(0, tmp_age_cat[["upper"]]),
                            labels = tmp_age_cat[["group"]],
                            right  = FALSE),
                        levels = tmp_age_cat[["group"]])
                          )
    ]


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

cli_alert("identifying cases and creating indicator variables...")
cli_progress_bar("comorbid ids", total = length(names(comorbid)))
for (i in names(comorbid)) {
  comorbid[[i]][["ids"]] <- ukb_full_phe[phecode %in% comorbid[[i]][["phecodes"]], id] |>
    unique()
  cli_progress_update()
}

cli_progress_bar("comorbid vars", total = length(names(comorbid)))
for (i in names(comorbid)) {
  set(ukb_demo, j = i, value = fifelse(ukb_demo[["id"]] %in% comorbid[[i]][["ids"]], 1, 0))
  cli_progress_update()
}

# summarize demo data ----------------------------------------------------------
(mgi_demo_summary <- mgi_cov[id %in% unique(mgi_full_phe[, id]), ] |>
    summarizer(col_names = c("age_at_last_diagnosis", "age_verbose", "sex", "race_eth", "bmi", "bmi_verbose",
                             "cancer", "diabetes", "cad", "anxiety", "depression", "SmokingStatus")))

(ukb_demo_summary <- ukb_demo[id %in% unique(ukb_full_phe[, id]), ] |>
    summarizer(col_names = c("age", "age_verbose", "sex", "race_eth",
                             "cancer", "diabetes", "cad", "anxiety", "depression")))

# summarize ehr data -----------------------------------------------------------
(mgi_ehr_summary <- phecode_dsb_summarizer(mgi_full_phe))

(ukb_ehr_summary <- phecode_dsb_summarizer(ukb_full_phe))

# stack summaries --------------------------------------------------------------
(mgi_stacked_summary <- rbindlist(list(
    mgi_demo_summary,
    mgi_ehr_summary
), use.names = TRUE, fill = TRUE))
(ukb_stacked_summary <- rbindlist(list(
    ukb_demo_summary,
    ukb_ehr_summary
), use.names = TRUE, fill = TRUE))

# save -------------------------------------------------------------------------
fwrite(
    x    = mgi_stacked_summary,
    file = glue("data/private/mgi/{opt$mgi_version}/mgi_demo_ehr_summary.txt"),
    sep  = "\t"
)
fwrite(
    x    = ukb_stacked_summary,
    file = glue("data/private/ukb/{opt$ukb_version}/ukb_demo_ehr_summary.txt"),
    sep  = "\t"
)
