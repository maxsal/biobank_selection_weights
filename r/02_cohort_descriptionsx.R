# Prepare MGI data including deriving variables for comorbidity status and
# descriptive variables derived from demographic data
# author:   max salvatore
# date:     20230809

options(stringsAsFactors = FALSE)

ms::libri(
  ms, data.table, qs, tidyverse, GGally, 
  scales, gridExtra, PheWAS/PheWAS, glue,
  survey, optparse, cli, parallelly,
  doParallel, foreach
)

option_list <- list(
  make_option("--mgi_version",
    type = "character", default = "20220822",
    help = "Cohort version in /net/junglebook/magic_data/EHRdata/ [default = %default]"
  ),
  make_option("--ukb_version",
    type = "character", default = "20221117",
    help = "Cohort version for UKB [default = %default]"
  ),
  make_option("--mgi_weight",
    type = "character", default = "ip_selection",
    help = "Weight variable to use for MGI [default = %default]"
  )
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args   <- parse_args(parser, positional_arguments = 0)
opt    <- args$options
print(opt)

for (i in c(
  "calculate_prevalences.R",
  "calculate_weighted_prevalences.R",
  "files-utils.R"
)) {
  source(paste0("fn/", i))
}

file_paths <- get_files(mgi_version = opt$mgi_version, ukb_version = opt$ukb_version)

# data
cli_alert("loading data...")
mgi_cov     <- qread(glue("data/private/mgi/{opt$mgi_version}/datax_{opt$mgi_version}_comb.qs"))
mgi_pim     <- qread(glue("data/private/mgi/{opt$mgi_version}/MGI_PIM0X_{opt$mgi_version}.qs"))[id %in% mgi_cov[, DeID_PatientID], ]
mgi_weights <- qread(glue("data/private/mgi/{opt$mgi_version}/weightsx_{opt$mgi_version}_comb.qs"))
mgi_pim     <- merge.data.table(mgi_pim, mgi_weights[, .(id, weights = get(opt$mgi_weight))], by.x = "id", by.y = "id", all.x = TRUE)

ukb_pim     <- read_qs(glue("/net/junglebook/home/mmsalva/projects/dissertation/aim_one/data/private/ukb/{opt$ukb_version}/UKB_PIM0X_{opt$ukb_version}.qs"))
ukb_cov     <- read_qs(glue("data/private/ukb/{opt$ukb_version}/datax_{opt$ukb_version}_comb.qs"))
ukb_weights <- fread("/net/junglebook/home/mmsalva/createUKBphenome/data/UKBSelectionWeights.tab",
                     colClasses = "character")[, .(id = f.eid, weights = as.numeric(LassoWeight))]
ukb_pim     <- merge.data.table(ukb_pim, ukb_weights, by = "id", all.x = TRUE)

pheinfo <- ms::pheinfox

# calculate prevalences --------------------------------------------------------
cli_alert("calculating unweighted prevalences...")
## unweighted
mgi_prevs <- calculate_prevalences(
  pim_data   = mgi_pim,
  cov_data   = mgi_cov,
  pim_id_var = "id",
  pheinfo    = ms::pheinfox
)
ukb_prevs <- calculate_prevalences(
  pim_data   = ukb_pim,
  cov_data   = ukb_cov,
  pheinfo    = ms::pheinfox,
  pim_id_var = "id",
  cov_id_var = "id",
  sex_var    = "sex",
  male_val   = "Male",
  female_val = "Female"
)

cli_alert("calculating weighted prevalences...")
## weighted
mgi_prevs_w <- calculate_weighted_prevalences(
  pim_data    = mgi_pim[!is.na(weights), ],
  cov_data    = mgi_cov,
  pim_id_var  = "id",
  pheinfo     = ms::pheinfox,
  weight      = "weights",
  parallelize = "doParallel",
  n_cores     = 8
)
ukb_prevs_w <- calculate_weighted_prevalences(
  pim_data    = ukb_pim[!is.na(weights), ],
  cov_data    = ukb_cov,
  pim_id_var  = "id",
  cov_id_var  = "id",
  sex_var     = "sex",
  male_val    = "Male",
  female_val  = "Female",
  pheinfo     = ms::pheinfox,
  weight      = "weights",
  parallelize = "doParallel",
  n_cores     = 8
)

# merge weighted and unweighted prevalences ------------------------------------
mgi_merged_prevs <- merge.data.table(
  mgi_prevs,
  mgi_prevs_w,
  by = "phecode",
  suffixes = c("_unweighted", "_weighted")
)
ukb_merged_prevs <- merge.data.table(
  ukb_prevs,
  ukb_prevs_w,
  by = "phecode",
  suffixes = c("_unweighted", "_weighted")
)

# save prevalences
save_qs(
  x    = mgi_merged_prevs,
  file = glue("results/mgi/{opt$mgi_version}/mgi_prevsx.qs")
)
save_qs(
  x    = ukb_merged_prevs,
  file = glue("results/ukb/{opt$ukb_version}/ukb_prevsx.qs")
)

cli_alert_success("done! 🎉")
