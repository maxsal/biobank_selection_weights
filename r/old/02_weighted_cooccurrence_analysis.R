# quickly perform phecode-phecode phewas using MGI data for multiple
# time-thresholds with weights
# author:   max salvatore
# date:     20230418

# libraries, functions, and options --------------------------------------------
ms::libri(
  ms, data.table, MatchIt, glue, qs, cli, optparse
)

set.seed(61787)

for (i in list.files("fn/", full.names = TRUE)) source(i)

# optparse list ----------------------------------------------------------------
option_list <- list(
  make_option("--outcome",
    type = "character", default = "153",
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
      "glm, logistf or SPAtest [default = %default]"
    )
  ),
  make_option("--weights",
    type = "character", default = "ip_selection_f,ps_nhw_f",
    help = glue(
      "Weighting variable to use for weighted analyses - ",
      "selection, all, or list of named weight variables [default = %default]"
    )
  )
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

time_thresholds <- as.numeric(strsplit(opt$time_thresholds, ",")[[1]])

## extract file paths
file_paths <- get_files(mgi_version = opt$mgi_version)

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

mgi_weights <- read_qs(glue("data/private/mgi/{opt$mgi_version}/weights_{opt$mgi_version}_{opt$mgi_cohort}.qs"))

## phenome
pheinfo <- fread("https://raw.githubusercontent.com/maxsal/public_data/main/phewas/Phecode_Definitions_FullTable_Modified.txt",
  colClasses = "character", showProgress = FALSE
)

# weights ----------------------------------------------------------------------
if (opt$weights == "all") {
  weight_vars <- names(mgi_weights)[!names(mgi_weights) %in% c("id", "DeID_PatientID")]
} else if (opt$weights == "selection") {
  weight_vars <- grep("selection_c", names(mgi_weights), value = TRUE)
} else {
  weight_vars <- unlist(strsplit(opt$weights, ","))
}

# merge data -------------------------------------------------------------------
weight_id <- c("id", weight_vars)

mgi_tr_merged <- lapply(
  names(mgi_tr_pims),
  \(x) {
    merge_list(list(mgi_tr_pims[[x]], mgi_covariates[, !c("case")], mgi_weights[, ..weight_id]), by_var = "id", join_fn = dplyr::left_join)[, `:=` (
      age_at_threshold = round(get(x) / 365.25, 1)
    )]
  }
)
names(mgi_tr_merged) <- glue("t{time_thresholds}_threshold")

# cooccurrence analysis --------------------------------------------------------
out <- list()
for (w in seq_along(weight_vars)) {
  cli_alert(glue("cooccurrence using {weight_vars[w]} [{w}/{length(weight_vars)}]..."))
  out[[w]] <- lapply(
    seq_along(time_thresholds),
    \(x) {
      print(x)
      cooccurrence_analysis(
        data       = mgi_tr_merged[[x]],
        covariates = c("age_at_threshold", "female", "length_followup"),
        weight_var = weight_vars[w],
        n_cores    = parallelly::availableCores() / 4,
        parallel   = TRUE
      )
    }
  )
}
names(out) <- weight_vars

# save results -----------------------------------------------------------------
for (w in seq_along(weight_vars)) {
  for (i in seq_along(time_thresholds)) {
    save_qs(
      x = out[[w]][[i]],
      file = glue(
        "results/mgi/{opt$mgi_version}/X{gsub('X', '', opt$outcome)}/",
        "mgi_X{gsub('X', '', opt$outcome)}_t{time_thresholds[i]}_",
        "{opt$mgi_version}_{weight_vars[w]}_results.qs"
      )
    )
  }
}

cli_alert_success("done! 🎉")
