# quickly perform phecode-phecode phewas using MGI data for multiple
# time-thresholds for a selected outcome variable
# author:   max salvatore
# date:     20230220

# libraries, functions, and options --------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(MatchIt)
  library(logistf)
  library(glue)
  library(qs)
  library(optparse)
})

set.seed(61787)

for (i in list.files("fn/", full.names = TRUE)) source(i)

# optparse list ----------------------------------------------------------------
option_list <- list(
  make_option("--outcome",
    type = "character", default = "157",
    help = "Outcome phecode [default = %default]"
  ),
  make_option("--mgi_version",
    type = "character", default = "20220822",
    help = "Version of MGI data [default = %default]"
  ),
  make_option("--mgi_cohort",
    type = "character", default = "comb",
    help = "Cohort of MGI used in weighting (comb, bb, mend, mhb) [default = %default]"
  ),
  make_option("--ukb_version",
    type = "character", default = "20221117",
    help = "Version of UKB data [default = %default]"
  ),
  make_option("--time_thresholds",
    type = "character", default = "0,0.5,1,2,3,5",
    help = glue(
      "Time thresholds for the phenome data ",
      "[default = %default]"
    )
  ),
  make_option("--mod_type",
    type = "character", default = "glm",
    help = glue(
      "Type of model to use in cooccurrence analysis - ",
      "glm, logistf, or SPAtest [default = %default]"
    )
  )
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

time_thresholds <- as.numeric(strsplit(opt$time_thresholds, ",")[[1]])

## extract file paths
file_paths <- get_files(
  mgi_version = opt$mgi_version,
  ukb_version = opt$ukb_version
)

# read data --------------------------------------------------------------------
## mgi
mgi_tr_pims <- lapply(
  seq_along(time_thresholds),
  \(x) {
    glue(
      "data/private/mgi/{opt$mgi_version}/X", "{gsub('X', '', opt$outcome)}/",
      "time_restricted_phenomes/mgi_X{gsub('X', '', opt$outcome)}_t",
      "{time_thresholds[x]}_{opt$mgi_version}.qs"
    ) |>
      read_qs()
  }
)
names(mgi_tr_pims) <- glue("t{time_thresholds}_threshold")

mgi_covariates <- read_qs(glue(
  "data/private/mgi/{opt$mgi_version}/X{gsub('X', '', opt$outcome)}/",
  "matched_covariates.qs"
))

## ukb
ukb_tr_pims <- lapply(
  seq_along(time_thresholds),
  \(x) {
    glue(
      "data/private/ukb/{opt$ukb_version}/X", "{gsub('X', '', opt$outcome)}/",
      "time_restricted_phenomes/ukb_X{gsub('X', '', opt$outcome)}_t",
      "{time_thresholds[x]}_{opt$ukb_version}.qs"
    ) |>
      read_qs()
  }
)
names(ukb_tr_pims) <- glue("t{time_thresholds}_threshold")

ukb_covariates <- read_qs(
  glue(
    "data/private/ukb/{opt$ukb_version}/X{gsub('X', '', opt$outcome)}/",
    "matched_covariates.qs"
  )
)

## phenome
pheinfo <- fread("data/public/Phecode_Definitions_FullTable_Modified.txt",
  colClasses = "character"
)

# cooccurrence analysis --------------------------------------------------------
## mgi
mgi_results <- lapply(
  seq_along(time_thresholds),
  \(x) output_cooccurrence_results(
    pim_data     = mgi_tr_pims[[glue("t{time_thresholds[x]}_threshold")]],
    t_thresh     = time_thresholds[x],
    cov_data     = mgi_covariates,
    covariates   = c("age_at_threshold", "female", "length_followup"),
    all_phecodes = glue("X{pheinfo[, phecode]}"),
    model_type   = opt$mod_type,
    parallel     = TRUE
  )
)
names(mgi_results) <- glue("t{time_thresholds}")

## ukb
ukb_results <- lapply(
  seq_along(time_thresholds),
  \(i) output_cooccurrence_results(
    pim_data = ukb_tr_pims[[glue("t{time_thresholds[i]}_threshold")]],
    t_thresh = time_thresholds[i],
    cov_data = ukb_covariates,
    covariates = c("age_at_threshold", "female", "length_followup"),
    all_phecodes = glue("X{pheinfo[, phecode]}"),
    model_type = opt$mod_type,
    parallel = TRUE
  )
)
names(ukb_results) <- glue("t{time_thresholds}")

# save results -----------------------------------------------------------------
## mgi
for (i in seq_along(time_thresholds)) {
  save_qs(
    x = mgi_results[[i]],
    file = glue(
      "results/mgi/{opt$mgi_version}/X{gsub('X', '', opt$outcome)}/",
      "mgi_X{gsub('X', '', opt$outcome)}_t{time_thresholds[i]}_",
      "{opt$mgi_version}_results.qs"
    )
  )
}

## ukb
for (i in seq_along(time_thresholds)) {
  save_qs(
    x = ukb_results[[i]],
    file = glue(
      "results/ukb/{opt$ukb_version}/X{gsub('X', '', opt$outcome)}/",
      "ukb_X{gsub('X', '', opt$outcome)}_t{time_thresholds[i]}_",
      "{opt$ukb_version}_results.qs"
    )
  )
}
