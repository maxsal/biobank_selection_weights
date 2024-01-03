# quickly perform phecode-phecode phewas using MGI data for multiple
# time-thresholds for a selected outcome variable
# requires: time-threshold phecode indicator matrices must already exist
# outputs:  betas, sebetas, and p-values
# author:   max salvatore
# date:     20221207

# libraries, functions, and options --------------------------------------------
library(data.table)
library(purrr)
library(logistf)
library(glue)
library(progress)

set.seed(61787)

for (i in list.files("fn/")) source(paste0("fn/", i)) # load functions

# specifications ---------------------------------------------------------------
outcome         <- "184.1"          # 155, liver cancer; 157, pancreatic cancer; 184.1, ovarian cancer
mgi_version     <- "20210318"       # mgi version
ukb_version     <- "20221117"       # ukb  version
time_thresholds <- c(0, 1, 2, 3, 5) # time thresholds
mod_type        <- "logistf"        # model type - use "SPAtest" or "logistf" for cooccur analyses

## extract file paths
file_paths <- get_files(mgi_version = mgi_version, ukb_version = ukb_version)

# read data --------------------------------------------------------------------
tr_pims <- list()
for (i in seq_along(time_thresholds)) {
  tr_pims[[i]] <- data.table::fread(
    glue::glue("data/{mgi_version}/processed/X","{gsub('X', '', outcome)}/",
               "time_restricted_phenomes/mgi_X{gsub('X', '', outcome)}_t",
               "{time_thresholds[i]}_{mgi_version}.txt")
    )
}
names(tr_pims) <- glue::glue("t{time_thresholds}_threshold")

covariates <- data.table::fread(
  glue::glue("data/{mgi_version}/processed/X{gsub('X', '', outcome)}/",
             "matched_covariates.txt")
  )
pheinfo <- fread(
  glue::glue("data/phecode_mapping/data/",
             "Phecode_Definitions_FullTable_Modified.txt"),
  colClasses = "character")

# cooccurrence analysis --------------------------------------------------------
results <- list()
for (i in seq_along(time_thresholds)) {
  results[[i]] <- output_mgi_cooccur_results(
    pim_data = tr_pims[[glue::glue("t{time_thresholds[i]}_threshold")]],
    t_thresh = time_thresholds[i],
    cov_data = covariates,
    covariates = c("age_at_threshold", "female", "length_followup"),
    all_phecodes = glue::glue("X{pheinfo[, phecode]}"),
    model_type = mod_type
  )
}
names(results) <- glue::glue("t{time_thresholds}")

# save results -----------------------------------------------------------------
for (i in seq_along(time_thresholds)) {
  data.table::fwrite(
    x    = results[[i]],
    file = glue::glue("results/{mgi_version}/X{gsub('X', '', outcome)}/mgi_X",
                      "{gsub('X', '', outcome)}_t{time_thresholds[i]}_",
                      "{mgi_version}_results.txt"),
    sep  = "\t"
  )
}

