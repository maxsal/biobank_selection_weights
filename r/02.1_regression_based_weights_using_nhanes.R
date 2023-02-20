# simplex regression-based weights
# requires: mgi data; nhanes wave info
# outputs:  nhanes and uncorrected cancer weights
# author:   max salvatore
# date:     20230203

# libraries --------------------------------------------------------------------
suppressPackageStartupMessages({
  library(haven)
  library(survey)
  library(dplyr)
  library(pracma)
  library(simplexreg)
  library(data.table)
  library(glue)
  library(fst)
  library(optparse)
})

# optparse list ---
option_list <- list(
  make_option("--cohort_version", type = "character", default = "20220822",
              help = "Cohort version in /net/junglebook/magic_data/EHRdata/ [default = '20220822']"),
  make_option("--mgi_cohort", type = "character", default = "comb",
              help = "Specific MGI cohort [default = %default]"),
  make_option("--nhanes_wave_letter", type = "character",
              default = "J",
              help = glue("NHANES data prefix corresponding to wave ",
                          "[default = %default]")),
  make_option("--nhanes_wave_years", type = "character",
              default = "2017-2018",
              help = glue("NHANES wave years corresponding to wave ",
                          "[default = %default]")),
  make_option("--nhanes_survey_names", type = "character",
              default = "DEMO,BMX,SMQ,DIQ,MCQ",
              help = glue("NHANES wave years corresponding to wave ",
                          "[default = %default]"))
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args   <- parse_args(parser, positional_arguments = 0)
opt    <- args$options

data_path <- glue("data/private/mgi/{opt$cohort_version}/")

lapply(list.files("fn", full.names = TRUE), source) |> # load functions
  invisible()

# load data --------------------------------------------------------------------
mgi <- read_fst(glue("{data_path}data_{opt$cohort_version}_{opt$mgi_cohort}.fst"),
               as.data.table = TRUE)
setnames(mgi, "DeID_PatientID", "id", skip_absent = TRUE)

nhanes_datasets <- unlist(strsplit(opt$nhanes_survey_names, ","))
nhanes_merged <- download_nhanes_data(
  wave_letter = opt$nhanes_wave_letter,
  wave_years  = opt$nhanes_wave_years,
  datasets    = nhanes_datasets
)

keep_vars <- c("SEQN", "RIAGENDR", "WTINT2YR", "RIDAGEYR", "RIDRETH1", "MCQ220",
               "BMXBMI", "SMQ040", "SMQ020", "DIQ010", "MCQ160C", "WTMEC2YR",
               "SDMVSTRA", "SDMVPSU")

if ("WTMECPRP" %in% names(nhanes_merged)) {
  setnames(nhanes_merged,
           "WTMECPRP",
           "WTMEC2YR")
}
if ("WTINTPRP" %in% names(nhanes_merged)) {
  setnames(nhanes_merged,
           "WTINTPRP",
           "WTINT2YR")
}

nhanes_merged <- nhanes_merged[, ..keep_vars]

prepped_nhanes <- prepare_nhanes_data(
  nhanes_data = nhanes_merged,
  mec_wt_var = "WTMEC2YR")

stacked <- rbindlist(list(
  prepped_nhanes,
  mgi[][, dataset := "MGI"]
), use.names = TRUE, fill = TRUE)

# estimate ipw and postratification weights ------------------------------------
cli_alert("estimating ipw weights...")
estimated_weights <- ipw(stacked_data = stacked)

cli_alert("estimating poststratification weights...")
post   <- poststratification(mgi_data = mgi, chop = TRUE)

merged <- Reduce(\(x, y) merge.data.table(x, y, by = "id"), list(mgi, estimated_weights, post))

# estimate cancer~female log(OR) -----------------------------------------------
cli_alert("estimating cancer~female log(OR) as sanity check...")
## ipw
m0 <- glm(cancer~as.numeric(Sex == "F"), family = quasibinomial(), data = merged)
m1 <- glm(cancer~as.numeric(Sex == "F"), family = quasibinomial(), data = merged, weights = no_cancer_ipw)
m2 <- glm(cancer~as.numeric(Sex == "F"), family = quasibinomial(), data = merged, weights = cancer_indirect_ipw)

m0_bb <- glm(cancer~as.numeric(Sex == "F"), family = quasibinomial(), data = merged[StudyName == "MGI", ])
m1_bb <- glm(cancer~as.numeric(Sex == "F"), family = quasibinomial(), data = merged[StudyName == "MGI", ], weights = no_cancer_ipw)
m2_bb <- glm(cancer~as.numeric(Sex == "F"), family = quasibinomial(), data = merged[StudyName == "MGI", ], weights = cancer_indirect_ipw)
##

## poststrat
m1_ps <- glm(cancer ~ as.numeric(Sex == "F"), family = quasibinomial(), data = merged, weights = no_cancer_postw)
m2_ps <- glm(cancer ~ as.numeric(Sex == "F"), family = quasibinomial(), data = merged, weights = cancer_postw)

m1_ps_bb <- glm(cancer ~ as.numeric(Sex == "F"), family = quasibinomial(), data = merged[StudyName == "MGI", ], weights = no_cancer_postw)
m2_ps_bb <- glm(cancer ~ as.numeric(Sex == "F"), family = quasibinomial(), data = merged[StudyName == "MGI", ], weights = cancer_postw)
##

extractr <- function(x, weight_name) {
  suppressMessages({y <- confint(x)})
  data.table(
    "weights" = weight_name,
    "est"     = coef(x)[[2]],
    "lower"   = y[2, 1],
    "upper"   = y[2, 2],
    "var"     = diag(summary(x)$cov.scaled)[[2]]
  )
}

log_or_est <- rbindlist(list(
  extractr(x = m0, weight_name = "None"),
  extractr(x = m1, weight_name = "No cancer"),
  extractr(x = m2, weight_name = "Cancer (uncorrected)"),
  extractr(x = m0_bb, weight_name = "None [BB]"),
  extractr(x = m1_bb, weight_name = "No cancer [BB]"),
  extractr(x = m2_bb, weight_name = "Cancer (uncorrected) [BB]"),
  extractr(x = m1_ps, weight_name = "Poststrat: without cancer"),
  extractr(x = m2_ps, weight_name = "Poststrat: with cancer"),
  extractr(x = m1_ps_bb, weight_name = "Poststrat: without cancer [BB]"),
  extractr(x = m2_ps_bb, weight_name = "Poststrat: with cancer [BB]")
))
log_or_est

fwrite(x = log_or_est, file = glue("{data_path}cancer_female_logor_est_{opt$cohort_version}_{opt$mgi_cohort}.csv"))

# save -------------------------------------------------------------------------
write_fst(
  x = merged[, .(id, no_cancer_ipw, cancer_indirect_ipw, no_cancer_postw, cancer_postw)],
  path = glue("{data_path}weights_{opt$cohort_version}_{opt$mgi_cohort}.fst")
  )

cli_alert_success("script success! see {.path {data_path}} and suffix {.emph {opt$mgi_cohort}}")
