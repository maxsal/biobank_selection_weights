# install and load packages -----------------------------------------------------
list_of_packages <- c(
    "data.table", "qs", "cli", "glue", "parallelly", "tidyverse",
    "bigrquery", "dplyr", "pracma", "simplexreg", "data.table",
    "lubridate", "colorblindr", "showtext", "optparse", "parallel", "parallelly",
    "tictoc"
)
new_packages <- list_of_packages[!(list_of_packages %in% installed.packages()[, "Package"])]
if (length(new_packages)) install.packages(new_packages)
for (i in list_of_packages) suppressPackageStartupMessages(library(i, character.only = TRUE))

font_add_google("Lato")
showtext_auto()
options(repr.plot.width = 10, repr.plot.height = 10, repr.plot.res = 100)

option_list <- list(
    make_option(c("-p", "--project_directory"),
        type = "character", default = "/home/jupyter/workspaces/describingehrlinkedbiobanks",
        help = "Path to project directory [default = %default]"
    )
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

setwd(opt$project_directory)

options(future.globals.maxSize=+Inf)

source("fn/aou_helpers.R")
source("fn/download_nhanes_data.R")
source("fn/prepare_nhanes_data.R")
source("fn/weight-utils.R")
source("fn/poststrat_nhanes.R")
source("fn/cleaning-utils.R")
source("fn/calculate_prevalences.R")
source("fn/calculate_weighted_prevalences.R")
source("fn/summarize-utils.R")
n_cores <- parallelly::availableCores()

# read data ---------------------------------------------------------------------
cli_alert("reading data...")
demo <- read_qs("data/processed/demo_cov_20230309.qs")
demo <- demo[!is.na(AgeLastEntry), ]

# estimate weights --------------------------------------------------------------
nhanes <- download_nhanes_data(
    datasets = c("DEMO", "BMX", "SMQ", "DIQ", "MCQ", "DPQ", "BPQ")
)

keep_vars <- c(
    "SEQN", "RIAGENDR", "WTINT2YR", "RIDAGEYR", "RIDRETH1", "MCQ220",
    "BMXBMI", "SMQ040", "SMQ020", "DIQ010", "MCQ160C", "WTMEC2YR",
    "SDMVSTRA", "SDMVPSU", paste0("DPQ0", 1:9, "0"), "BPQ020"
)

if ("WTMECPRP" %in% names(nhanes)) {
    setnames(
        nhanes_merged,
        "WTMECPRP",
        "WTMEC2YR"
    )
}
if ("WTINTPRP" %in% names(nhanes)) {
    setnames(
        nhanes_merged,
        "WTINTPRP",
        "WTINT2YR"
    )
}

nhanes_merged <- nhanes[, ..keep_vars]

prepped_nhanes <- prepare_nhanes_data(
    nhanes_data = nhanes,
    mec_wt_var = "WTMEC2YR"
)

stacked <- rbindlist(list(
    prepped_nhanes,
    demo[][, dataset := "AOU"]
), use.names = TRUE, fill = TRUE)

ipw_weights <- ipw(
    stacked_data = stacked,
    dataset_name = "AOU",
    covs = c(
        "as.numeric(age_cat == 5)", "as.numeric(age_cat == 6)",
        "cad", "diabetes", "smoking_current", "smoking_former",
        "bmi_under", "bmi_overweight", "bmi_obese", "nhanes_nhw",
        "cancer", "female"
    )
)

demo[, `:=`(
    smoker = fcase(
        SmokingStatus %in% c("Current", "Past"), 1,
        SmokingStatus == "Never", 0,
        default = NA_real_
    ),
    nhw = nhanes_nhw
)]

post <- poststrat_nhanes(
    int_data       = demo,
    chop           = TRUE,
    age_bin        = TRUE,
    age_bin_breaks = c(seq(0, 80, 10), 150),
    covs           = c("smoker", "cad", "diabetes", "nhw", "female")
)

ipw_weights <- na.omit(ipw_weights)
post <- na.omit(post)

weights <- merge.data.table(ipw_weights, post, by = "id", all = TRUE)

demow <- merge.data.table(
    demo,
    weights,
    by = "id"
)

# calculate prevalences ---------------------------------------------------------
cli_alert_info("calculating prevalences...")
# read pim
pim <- read_qs(
    file      = "data/processed/phenome/aou_20230309_phenome_pim0.qs",
    n_threads = n_cores
)
pim <- merge.data.table(
    pim,
    weights[, .(person_id = id, weights = ps_weight)],
    by = "person_id"
)

tic()
aou_prevs_w <- calculate_weighted_prevalences(
    pim_data    = pim[!is.na(weights), ],
    cov_data    = demo,
    pim_id_var   = "person_id",
    cov_id_var   = "id",
    sex_var = "sex",
    male_val     = "Male",
    female_val   = "Female",
    weight      = "weights",
    n_cores     = 2,
    parallelize = "future"
)
toc()

print(aou_prevs_w)

cli_alert_success("script success! 🎉")
