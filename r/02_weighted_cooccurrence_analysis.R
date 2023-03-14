# quickly perform phecode-phecode phewas using MGI data for multiple
# time-thresholds with weights
# requires: time-threshold phecode indicator matrices, weights must already exist
# outputs:  betas, sebetas, and p-values
# author:   max salvatore
# date:     20230220

# libraries, functions, and options --------------------------------------------
library(data.table)
library(MatchIt)
library(logistf)
library(glue)
library(fst)
library(cli)
library(optparse)

set.seed(61787)

lapply(list.files("fn/", full.names = TRUE), source) |> # load functions
  invisible()

# optparse list ----------------------------------------------------------------
option_list <- list(
  make_option("--outcome", type = "character", default = "157",
              help = "Outcome phecode [default = %default]"),
  make_option("--mgi_version", type = "character", default = "20220822",
              help = "Version of MGI data [default = %default]"),
  make_option("--mgi_cohort", type = "character", default = "comb",
              help = "Cohort of MGI used in weighting (comb, bb, mend, mhb) [default = %default]"),
  make_option("--time_thresholds", type = "character", default = "0,0.5,1,2,3,5",
              help = glue("Time thresholds for the phenome data ",
                          "[default = %default]")),
  make_option("--mod_type", type = "character", default = "glm",
              help = glue("Type of model to use in cooccurrence analysis - ",
                          "glm, logistf or SPAtest [default = %default]")),
  make_option("--weights", type = "character", default = "all",
              help = glue("Weighting variable to use for weighted analyses - ",
                          "no_cancer_ipw, cancer_indirect_ipw, no_cancer_postw, cancer_postw or all [default = %default]"))
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args   <- parse_args(parser, positional_arguments = 0)
opt    <- args$options
print(opt)

time_thresholds <- as.numeric(strsplit(opt$time_thresholds, ",")[[1]])

## extract file paths
file_paths <- get_files(mgi_version = opt$mgi_version)

# read data --------------------------------------------------------------------
## mgi
mgi_tr_pims <- lapply(seq_along(time_thresholds),
                      \(x) {glue("data/private/mgi/{opt$mgi_version}/X","{gsub('X', '', opt$outcome)}/",
                                 "time_restricted_phenomes/mgi_X{gsub('X', '', opt$outcome)}_t",
                                 "{time_thresholds[x]}_{opt$mgi_version}.fst") |>
                          read_fst(as.data.table = TRUE)})
names(mgi_tr_pims) <- glue("t{time_thresholds}_threshold")

mgi_covariates <- read_fst(glue("data/private/mgi/{opt$mgi_version}/X{gsub('X', '', opt$outcome)}/",
                                "matched_covariates.fst"),
                           as.data.table = TRUE)

mgi_weights <- read_fst(glue("data/private/mgi/{opt$mgi_version}/weights_{opt$mgi_version}_{opt$mgi_cohort}.fst"),
                        as.data.table = TRUE)

## phenome
pheinfo <- fread("data/public/Phecode_Definitions_FullTable_Modified.txt",
                 colClasses = "character")

# weights ----------------------------------------------------------------------
if (opt$weights == "all") {
  weight_vars <- names(mgi_weights)[!names(mgi_weights) %in% c("id", "DeID_PatientID")]
} else {
  weight_vars <- unlist(strsplit(opt$weights, ","))
}

# cooccurrence analysis --------------------------------------------------------
out <- list()
for (w in seq_along(weight_vars)) {
  
  cli_alert("cooccurrence using {weight_vars[w]}...")
  
  out[[w]] <- lapply(
    seq_along(time_thresholds),
    \(i) output_cooccurrence_results(
      pim_data     = mgi_tr_pims[[glue("t{time_thresholds[i]}_threshold")]],
      t_thresh     = time_thresholds[i],
      cov_data     = mgi_covariates,
      covariates   = c("age_at_threshold", "female", "length_followup"),
      all_phecodes = glue("X{pheinfo[, phecode]}"),
      model_type   = opt$mod_type,
      w_data       = mgi_weights,
      w_var        = weight_vars[w],
      parallel     = TRUE
    )
  )
  
}
names(out) <- weight_vars

# save results -----------------------------------------------------------------
## mgi
for (w in seq_along(weight_vars)) {
  lapply(
    seq_along(time_thresholds),
    \(i) {
      write_fst(
        x = out[[w]][[i]],
        path = glue("results/mgi/{opt$mgi_version}/X{gsub('X', '', opt$outcome)}/",
                    "mgi_X{gsub('X', '', opt$outcome)}_t{time_thresholds[i]}_",
                    "{opt$mgi_version}_{weight_vars[w]}_results.fst")
      )
    }
  ) |> invisible()
}
