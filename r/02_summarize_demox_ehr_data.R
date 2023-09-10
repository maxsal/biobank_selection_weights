# quickly perform phecode-phecode phewas using MGI data for multiple
# time-thresholds for a selected outcome variable
# author:   max salvatore
# date:     20230220

# libraries, functions, and options --------------------------------------------
ms::libri(
    data.table, MatchIt, logistf, glue, qs, optparse, purrr,
    ms, cli
)

set.seed(61787)

for (i in list.files("fn/", full.names = TRUE)) source(i)

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

## extract file paths
file_paths <- get_files(mgi_version = opt$mgi_version,
                        ukb_version = opt$ukb_version)

# read data --------------------------------------------------------------------
cli_alert("loading mgi data...")
## mgi
### demographics
mgi_cov <- read_qs(glue("data/private/mgi/{opt$mgi_version}/datax_{opt$mgi_version}_comb.qs"))
setnames(mgi_cov,
         old = c("DeID_PatientID", "Age", "AliveYN", "Deceased_DaysSinceBirth",
                 "Ethnicity", "MaritalStatusCode", "Sex", "Race",
                 "YearsInEHR", "FirstDaySinceBirth", "LastDaySinceBirth"),
         new = c("id", "age", "alive", "dead_dsb", "ethn", "marital",
                 "sex", "race", "length_followup", "first_dsb", "last_dsb"))

### icd-phecode data
mgi_full_phe <- qread(glue("data/private/mgi/{opt$mgi_version}/MGI_FULL_PHECODEX_DSB_{opt$mgi_version}.qs"))
# mgi_full_phe <- get(load(file_paths[["mgi"]][["phecode_dsb_file"]]))
if ("IID" %in% names(mgi_full_phe)) setnames(mgi_full_phe, "IID", "id")
if ("DaysSinceBirth" %in% names(mgi_full_phe)) setnames(mgi_full_phe, "DaysSinceBirth", "dsb")
mgi_full_phe <- mgi_full_phe[id %in% mgi_cov[, unique(id)], ]

### weights data
mgi_weights <- read_qs(glue("data/private/mgi/{opt$mgi_version}/weightsx_{opt$mgi_version}_comb.qs"))

mgi <- merge.data.table(
    mgi_cov,
    mgi_weights,
    by = "id",
    all.x = TRUE
)

## ukb
cli_alert("loading ukb data...")
### demographics
ukb_demo <- read_qs(glue("data/private/ukb/{opt$ukb_version}/datax_{opt$ukb_version}_comb.qs"))
ukb_demo <- ukb_demo[, .(
  id,
  dob  = birth_date,
  age  = age_at_consent,
  race_eth,
  sex,
  smoker = smk_status,
  bmi = bmi_med,
  bmi_cat = bmi_med_cat,
  drinker = as.numeric(alc_ev == "Ever"),
  cancer, diabetes, cad, anxiety, depression, in_phenome)]
# cc_vars <- c("id", "age", "ethn", "sex", "smoker", "bmi", "drinker")
# ukb_demo <- ukb_demo[complete.cases(ukb_demo[, ..cc_vars]), ]

### icd-phecode data
ukb_full_phe <- read_qs(glue("data/private/ukb/{opt$ukb_version}/UKB_FULL_PHECODEX_DSB_{opt$ukb_version}.qs"))
if ("IID" %in% names(ukb_full_phe)) setnames(ukb_full_phe, "IID", "id")
if ("DaysSinceBirth" %in% names(ukb_full_phe)) setnames(ukb_full_phe, "DaysSinceBirth", "dsb")
ukb_full_phe[, dsb := as.numeric(dsb)]

### weights data
ukb_weights <- fread("/net/junglebook/home/mmsalva/createUKBphenome/data/UKBSelectionWeights.tab",
    colClasses = "character")[, .(id = f.eid, ip_weight = as.numeric(LassoWeight))]

ukb <- merge.data.table(
    ukb_demo[id %in% unique(ukb_full_phe[, id]), ],
    ukb_weights,
    by = "id",
    all.x = TRUE
)

## other
cancer_phecodes <- fread("https://raw.githubusercontent.com/maxsal/public_data/main/phewas/cancer_phecodesx.csv",
                         colClasses = "character")[keep == 1, phecode]

# prep data --------------------------------------------------------------------
## age categories
age_cats    <- seq(0, 80, 10)
tmp_age_cat <- age_grp_table(
    lower_ages   = age_cats,
    num_vec      = rep(NA, length(age_cats)),
    num_var_name = "counts"
)
## mgi
mgi[, `:=` (
    race_eth = fcase(
        race == "Caucasian" & ethn == "Non-Hispanic", "NH White",
        race == "African American" & ethn == "Non-Hispanic", "NH Black",
        race == "Asian" & ethn == "Non-Hispanic", "NH Asian",
        ethn == "Hispanic", "Hispanic",
        default = "Other/Unknown"
    ),
    SmokingStatus = relevel(factor(SmokingStatus, levels = c("Never", "Past", "Current", "Unknown")), ref = "Never"),
    bmi_verbose   = factor(fcase(
        bmi_cat == 1, "Underweight (<18.5)",
        bmi_cat == 2, "Healthy [18.5, 25)",
        bmi_cat == 3, "Overweight [25, 30)",
        bmi_cat == 4, "Obese [30+)"
    ), levels = c("Underweight (<18.5)", "Healthy [18.5, 25)", "Overweight [25, 30)", "Obese [30+)")),
    age_verbose = factor(cut(age_at_last_diagnosisx,
                              breaks = c(0, tmp_age_cat[["upper"]]),
                              labels = tmp_age_cat[["group"]],
                              right  = FALSE),
                          levels = tmp_age_cat[["group"]])
)]
## ukb
ukb[, `:=` (
    age_verbose = factor(cut(age,
                            breaks = c(0, tmp_age_cat[["upper"]]),
                            labels = tmp_age_cat[["group"]],
                            right  = FALSE),
                        levels = tmp_age_cat[["group"]])
                          )
    ]

# unweighted demographics summary ----------------------------------------------
(mgi_demo_summary <- mgi[id %in% unique(mgi_full_phe[, id]), ] |>
    summarizer(col_names = c("age_at_last_diagnosisx", "age_verbose", "sex", "race_eth", "bmi", "bmi_verbose",
                             "cancerx", "diabetesx", "cadx", "anxietyx", "depressionx", "smoker")))
(ukb_demo_summary <- ukb[in_phenome == 1, ] |>
    summarizer(col_names = c("age", "age_verbose", "sex", "race_eth", "bmi", "bmi_cat",
                             "cancer", "diabetes", "cad", "anxiety", "depression", "smoker", "drinker")))

# weighted demographics summary ------------------------------------------------
(mgi_demo_summary_ip <- weighted_summary_wrapper(
    mgi,
    weight = "ip_selection",
    vars = c(
        "age_at_last_diagnosisx", "age_verbose", "sex", "race_eth", "bmi", "bmi_verbose",
        "cancerx", "diabetesx", "cadx", "anxietyx", "depressionx", "smoker"
    )
))
(mgi_demo_summary_ps <- weighted_summary_wrapper(
    mgi,
    weight = "ps_selection",
    vars = c(
        "age_at_last_diagnosisx", "age_verbose", "sex", "race_eth", "bmi", "bmi_verbose",
        "cancerx", "diabetesx", "cadx", "anxietyx", "depressionx", "smoker"
    )
))
(ukb_demo_summary_w <- weighted_summary_wrapper(
    data = ukb[in_phenome == 1, ],
    weight = "ip_weight",
    vars = c(
        "age", "age_verbose", "sex", "race_eth", "bmi", "bmi_cat",
        "cancer", "diabetes", "cad", "anxiety", "depression", "smoker", "drinker"
    )
))

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
    x    = mgi_demo_summary_ip,
    file = glue("data/private/mgi/{opt$mgi_version}/mgi_demo_summary_ip.txt"),
    sep  = "\t"
)
fwrite(
    x    = mgi_demo_summary_ps,
    file = glue("data/private/mgi/{opt$mgi_version}/mgi_demo_summary_ps.txt"),
    sep  = "\t"
)
fwrite(
    x    = ukb_stacked_summary,
    file = glue("data/private/ukb/{opt$ukb_version}/ukb_demo_ehr_summary.txt"),
    sep  = "\t"
)
fwrite(
    x    = ukb_demo_summary_w,
    file = glue("data/private/ukb/{opt$ukb_version}/ukb_demo_summary_w.txt"),
    sep  = "\t"
)

cli_alert_success("done! 🎉")
