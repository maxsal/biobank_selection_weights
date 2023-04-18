# Prepare MGI data including deriving variables for comorbidity status and
# descriptive variables derived from demographic data
# author:   max salvatore
# date:     20230418

suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(tidyverse)
  library(GGally)
  library(ggnetwork)
  library(network)
  library(scales)
  library(gridExtra)
  library(qgraph)
  library(igraph)
  library(circlize)
  library(ComplexHeatmap)
  library(PheWAS)
  library(glue)
  library(optparse)
})


option_list <- list(
  make_option("--mgi_version",
    type = "character", default = "20220822",
    help = "Cohort version in /net/junglebook/magic_data/EHRdata/ [default = %default]"
  ),
  make_option("--ukb_version",
    type = "character", default = "20221117",
    help = "Cohort version for UKB [default = %default]"
  )
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

for (i in c(
  "phecode_partial_correlation_chord_diagram.R",
  "phenome_partial_correlation_network.R",
  "phenome_partial_correlation_neoqgraph.R",
  "calculate_prevalences.R",
  "files-utils.R"
)) {
  source(paste0("fn/", i))
}

file_paths <- get_files(mgi_version = opt$mgi_version, ukb_version = opt$ukb_version)

# data
mgi <- read_qs(glue(
  "data/private/mgi/{opt$mgi_version}/",
  "mgi_phenome_partial_correlations_w_geno_pcs_{opt$mgi_version}.qs"
))
mgi_pim <- fread(gsub("_0", "", file_paths[["mgi"]][["pim0_file"]]))
mgi_cov <- fread(file_paths[["mgi"]][["cov_file"]])

ukb <- read_qs(glue(
  "data/private/ukb/{opt$ukb_version}/",
  "ukb_phenome_partial_correlations_{opt$ukb_version}.qs"
))
ukb_pim <- fread(file_paths[["ukb"]][["pim0_file"]])
ukb_cov <- fread(file_paths[["ukb"]][["demo_file"]])

pheinfo <- fread("https://gitlab.com/maxsal/public_data/-/raw/main/phewas/Phecode_Definitions_FullTable_Modified.txt",
  colClasses = "character", showProgress = FALSE
)

# calculate prevalences
mgi_prevs <- calculate_prevalences(pim_data = mgi_pim, cov_data = mgi_cov)
ukb_prevs <- calculate_prevalences(
  pim_data = ukb_pim,
  cov_data = ukb_cov,
  cov_id_var = "id",
  sex_var = "sex",
  male_val = "Male",
  female_val = "Female"
)

# save prevalences
save_qs(
  x    = mgi_prevs,
  file = glue("results/mgi/{opt$mgi_version}/mgi_prevs.qs")
)
save_qs(
  x    = ukb_prevs,
  file = blue("results/ukb/{opt$ukb_version}/ukb_prevs.qs")
)

# network plot
phenome_partial_correlation_network(
  x          = mgi,
  savefile   = glue("results/mgi/{opt$mgi_version}/MGI_network.pdf"),
  prevs      = mgi_prevs,
  plot_title = "MGI correlations"
)
phenome_partial_correlation_network(
  x          = ukb,
  savefile   = glue("results/ukb/{opt$ukb_version}/UKB_network.pdf"),
  prevs      = ukb_prevs,
  plot_title = "UKB correlations"
)

# chord diagram
phenome_partial_correlation_chord_diagram(
  x = mgi,
  savefile = glue("results/mgi/{opt$mgi_version}/MGI_partialchord.pdf")
)
phenome_partial_correlation_chord_diagram(
  x = ukb,
  savefile = glue("results/ukb/{opt$ukb_version}/UKB_partialchord.pdf")
)

# neoplasm qgraph
phenome_partial_correlation_neoqgraph(
  x        = mgi,
  savefile = glue("results/mgi/{opt_mgi_version}/MGI_qgraph.pdf")
)
phenome_partial_correlation_neoqgraph(
  x        = ukb,
  savefile = glue("results/ukb/{opt$ukb_version}/ukb_qgraph.pdf")
)
