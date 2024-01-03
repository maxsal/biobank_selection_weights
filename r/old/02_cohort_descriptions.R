# Prepare MGI data including deriving variables for comorbidity status and
# descriptive variables derived from demographic data
# author:   max salvatore
# date:     20230809

options(stringsAsFactors = FALSE)

ms::libri(
  ms, data.table, qs, tidyverse, GGally,
  ggnetwork, network, scales, gridExtra,
  qgraph, igraph, circlize, ComplexHeatmap,
  PheWAS/PheWAS, glue, survey, optparse,
  cli
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
    type = "character", default = "ip_selection_f",
    help = "Weight variable to use for MGI [default = %default]"
  )
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args   <- parse_args(parser, positional_arguments = 0)
opt    <- args$options
print(opt)

for (i in c(
  "phecode_partial_correlation_chord_diagram.R",
  "phenome_partial_correlation_network.R",
  "phenome_partial_correlation_neoqgraph.R",
  "calculate_prevalences.R",
  "calculate_weighted_prevalences.R",
  "files-utils.R"
)) {
  source(paste0("fn/", i))
}

file_paths <- get_files(mgi_version = opt$mgi_version, ukb_version = opt$ukb_version)

# data
cli_alert("loading data...")
mgi <- read_qs(glue(
  "data/private/mgi/{opt$mgi_version}/",
  "mgi_phenome_partial_correlations_{opt$mgi_version}.qs"
))
mgi_cov     <- qread(file_paths[["mgi"]][["cov_processed_file"]])
mgi_pim     <- qread(file_paths[["mgi"]][["pim0_file"]])[IID %in% mgi_cov[, DeID_PatientID], ]
mgi_weights <- read_qs(paste0("data/private/mgi/", opt$mgi_version, "/weights_", opt$mgi_version, "_comb.qs"))
mgi_pim     <- merge.data.table(mgi_pim, mgi_weights[, .(id, weights = get(opt$mgi_weight))], by.x = "IID", by.y = "id", all.x = TRUE)

ukb <- read_qs(glue(
  "data/private/ukb/{opt$ukb_version}/",
  "ukb_phenome_partial_correlations_{opt$ukb_version}.qs"
))

ukb_pim     <- read_qs(glue("/net/junglebook/home/mmsalva/projects/dissertation/aim_one/data/private/ukb/{opt$ukb_version}/UKB_PHENOME_PIM0_{opt$ukb_version}.qs"))
ukb_cov     <- read_qs(file_paths[["ukb"]][["demo_file"]])
ukb_weights <- fread("/net/junglebook/home/mmsalva/createUKBphenome/data/UKBSelectionWeights.tab",
                     colClasses = "character")[, .(id = f.eid, weights = as.numeric(LassoWeight))]
ukb_pim     <- merge.data.table(ukb_pim, ukb_weights, by = "id", all.x = TRUE)

pheinfo <- ms::pheinfo

# calculate prevalences --------------------------------------------------------
cli_alert("calculating unweighted prevalences...")
## unweighted
mgi_prevs <- calculate_prevalences(
  pim_data = mgi_pim,
  cov_data = mgi_cov
)
ukb_prevs <- calculate_prevalences(
  pim_data   = ukb_pim,
  cov_data   = ukb_cov,
  cov_id_var = "id",
  sex_var    = "sex",
  male_val   = "Male",
  female_val = "Female"
)

cli_alert("calculating weighted prevalences...")
## weighted
mgi_prevs_w <- calculate_weighted_prevalences(
  pim_data = mgi_pim[!is.na(weights), ],
  cov_data = mgi_cov,
  weight   = "weights"
)
ukb_prevs_w <- calculate_weighted_prevalences(
  pim_data = ukb_pim[!is.na(weights), ],
  cov_data = ukb_cov,
  weight   = "weights",
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
  file = glue("results/mgi/{opt$mgi_version}/mgi_prevs.qs")
)
save_qs(
  x    = ukb_merged_prevs,
  file = glue("results/ukb/{opt$ukb_version}/ukb_prevs.qs")
)

# network plot
cli_alert("plotting networks...")
phenome_partial_correlation_network(
  x          = mgi,
  savefile   = glue("results/mgi/{opt$mgi_version}/MGI_network.pdf"),
  prevs      = mgi_prevs,
  prev_var   = "prev_unweighted",
  plot_title = "MGI correlations"
)
phenome_partial_correlation_network(
  x          = ukb,
  savefile   = glue("results/ukb/{opt$ukb_version}/UKB_network.pdf"),
  prevs      = ukb_prevs,
  prev_var   = "prev_unweighted",
  plot_title = "UKB correlations"
)

# chord diagram
cli_alert("plotting chord diagrams...")
phenome_partial_correlation_chord_diagram(
  x = mgi,
  savefile = glue("results/mgi/{opt$mgi_version}/MGI_partialchord.pdf")
)
phenome_partial_correlation_chord_diagram(
  x = ukb,
  savefile = glue("results/ukb/{opt$ukb_version}/UKB_partialchord.pdf")
)

# neoplasm qgraph
cli_alert("plotting qgraphs...")
phenome_partial_correlation_neoqgraph(
  x        = mgi,
  savefile = glue("results/mgi/{opt$mgi_version}/MGI_qgraph.pdf")
)
phenome_partial_correlation_neoqgraph(
  x        = ukb,
  savefile = glue("results/ukb/{opt$ukb_version}/ukb_qgraph.pdf")
)

cli_alert_success("done! 🎉")
