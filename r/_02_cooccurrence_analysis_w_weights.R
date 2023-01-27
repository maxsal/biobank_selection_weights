# quickly perform phecode-phecode phewas using MGI data for multiple
# time-thresholds for a selected outcome variable
# requires: time-threshold phecode indicator matrices must already exist
# outputs:  betas, sebetas, and p-values
# author:   max salvatore
# date:     20230127

# libraries, functions, and options --------------------------------------------
library(data.table)
library(MatchIt)
library(logistf)
library(glue)
library(progress)
library(optparse)

set.seed(61787)

for (i in list.files("fn/")) source(paste0("fn/", i)) # load functions

# optparse list ----------------------------------------------------------------
option_list <- list(
  make_option("--outcome", type = "character", default = "",
              help = "Outcome phecode"),
  make_option("--mgi_version", type = "character", default = "20210318",
              help = "Version of MGI data [default = 20210318]"),
  make_option("--ukb_version", type = "character", default = "20221117",
              help = "Version of UKB data [default = 20221117]"),
  make_option("--time_thresholds", type = "character", default = "0,1,2,3,5",
              help = glue("Time thresholds for the phenome data ",
                          "[default = 0,1,2,3,5]")),
  make_option("--mod_type", type = "character", default = "logistf",
              help = glue("Type of model to use in cooccurrence analysis - ",
                          "logistf for SPAtest [default = logistf]")),
  make_option("--weighted", type = "logical", default = "",
              help = glue("Weighting variable to use for weighted analyses - ",
                          "weights_no_can [default = logistf]"))
)

parser <- OptionParser(usage="%prog [options]", option_list = option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

time_thresholds <- as.numeric(strsplit(opt$time_thresholds, ",")[[1]])

## extract file paths
file_paths <- get_files(mgi_version = opt$mgi_version,
                        ukb_version = opt$ukb_version)

# read data --------------------------------------------------------------------
## mgi
mgi_tr_pims <- list()
for (i in seq_along(time_thresholds)) {
  mgi_tr_pims[[i]] <- fread(
    glue("data/private/mgi/{opt$mgi_version}/X","{gsub('X', '', opt$outcome)}/",
         "time_restricted_phenomes/mgi_X{gsub('X', '', opt$outcome)}_t",
         "{time_thresholds[i]}_{opt$mgi_version}.txt")
    )
}
names(mgi_tr_pims) <- glue("t{time_thresholds}_threshold")

mgi_covariates <- fread(
  glue("data/private/mgi/{opt$mgi_version}/X{gsub('X', '', opt$outcome)}/",
             "matched_covariates.txt")
  )

## ukb
ukb_tr_pims <- list()
for (i in seq_along(time_thresholds)) {
  ukb_tr_pims[[i]] <- fread(
    glue("data/private/ukb/{opt$ukb_version}/X","{gsub('X', '', opt$outcome)}/",
               "time_restricted_phenomes/ukb_X{gsub('X', '', opt$outcome)}_t",
               "{time_thresholds[i]}_{opt$ukb_version}.txt")
  )
}
names(ukb_tr_pims) <- glue("t{time_thresholds}_threshold")

ukb_covariates <- fread(
  glue("data/private/ukb/{opt$ukb_version}/X{gsub('X', '', opt$outcome)}/",
             "matched_covariates.txt")
)

## phenome
pheinfo <- fread("data/public/Phecode_Definitions_FullTable_Modified.txt",
  colClasses = "character")

# cooccurrence analysis --------------------------------------------------------
## mgi
mgi_results <- list()
for (i in seq_along(time_thresholds)) {
  mgi_results[[i]] <- output_cooccurrence_results(
    pim_data = mgi_tr_pims[[glue("t{time_thresholds[i]}_threshold")]],
    t_thresh = time_thresholds[i],
    cov_data = mgi_covariates,
    covariates = c("age_at_threshold", "female", "length_followup"),
    all_phecodes = glue("X{pheinfo[, phecode]}"),
    model_type = opt$mod_type
  )
}
names(mgi_results) <- glue("t{time_thresholds}")

## ukb
ukb_results <- list()
for (i in seq_along(time_thresholds)) {
  ukb_results[[i]] <- output_cooccurrence_results(
    pim_data = ukb_tr_pims[[glue("t{time_thresholds[i]}_threshold")]],
    t_thresh = time_thresholds[i],
    cov_data = ukb_covariates,
    covariates = c("age_at_threshold", "female", "length_followup"),
    all_phecodes = glue("X{pheinfo[, phecode]}"),
    model_type = opt$mod_type
  )
}
names(ukb_results) <- glue("t{time_thresholds}")

# save results -----------------------------------------------------------------
## mgi
for (i in seq_along(time_thresholds)) {
  fwrite(
    x    = mgi_results[[i]],
    file = glue("results/mgi/{opt$mgi_version}/X{gsub('X', '', opt$outcome)}/",
                "mgi_X{gsub('X', '', opt$outcome)}_t{time_thresholds[i]}_",
                "{opt$mgi_version}_results.txt"),
    sep  = "\t"
  )
}

## ukb
for (i in seq_along(time_thresholds)) {
  fwrite(
    x    = ukb_results[[i]],
    file = glue("results/ukb/{opt$ukb_version}/X{gsub('X', '', opt$outcome)}/",
                "ukb_X{gsub('X', '', opt$outcome)}_t{time_thresholds[i]}_",
                "{opt$ukb_version}_results.txt"),
    sep  = "\t"
  )
}
