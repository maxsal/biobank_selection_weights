# Calculate partial correlations in UKB
# outputs: table of pairwise partial correlations
# author:  max salvatore
# date:    20230418

# 1. libraries, functions, and options -----------------------------------------
options(stringsAsFactors = FALSE)

ms::libri(
  ms, data.table, qs, glue, cli, doMC, optparse, parallelly, foreach
)
set.seed(61787)

source("fn/files-utils.R") # load partial correlation function

# optparse list ----------------------------------------------------------------
option_list <- list(
  make_option("--use_geno", type = "logical", default = FALSE,
              help = "Adjust for genotype PCs [default = %default]"),
  make_option("--ukb_version", type = "character", default = "20221117",
              help = "Version of MGI data [default = %default]"),
  make_option("--core_prop",
    type = "numeric", default = "0.125",
    help = "Proportion of available cores to use for parallel work [default = %default]"
  )
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args   <- parse_args(parser, positional_arguments = 0)
opt    <- args$options
print(opt)

file_paths <- get_files(ukb_version = opt$ukb_version)
n_cores    <- availableCores() * opt$core_prop

# 3. read data -----------------------------------------------------------------
cli_alert("reading data...")
## phecode indicator matrix (PEDMASTER_0)
pim0 <- read_qs(file_paths[["ukb"]][["pim0_file"]])

## demographics data
demo <- read_qs(file_paths[["ukb"]][["demo_file"]])[
  , .(id, age = age_at_consent, female = as.numeric(sex == "Female"))
] |> na.omit()

## pc data
if (opt$use_geno == TRUE) {
  stop("script not ready to incorporate UKB genotype PCs data")
  pcs <- fread(glue("/net/junglebook/magic_data/MGI_GenotypeData_Freeze4/",
                    "MGI_Freeze4_Sample_Info_60215samples.txt"),
                    colClasses = "character")
  setnames(pcs, old = "Deid_ID", new = "id")
  # merge pcs in demo data
  merged <- merge.data.table(demo, pcs)
} else {
  merged <- demo
}

use_these_ids <- intersect(pim0[, id], merged[, id])

sub_pim <- pim0[id %in% use_these_ids, ]    # subset pim
merged  <- merged[id %in% use_these_ids, ] # subset merged

# sex-specific
male_ids   <- demo[female == 0, id]
female_ids <- demo[female == 1, id]

male_vars   <- intersect(paste0("X", ms::pheinfo[sex == "Male", phecode]), names(sub_pim))
female_vars <- intersect(paste0("X", ms::pheinfo[sex == "Female", phecode]), names(sub_pim))

sub_pim[id %in% male_ids, (female_vars) := lapply(.SD, \(x) NA), .SDcols = female_vars]
sub_pim[id %in% female_ids, (male_vars) := lapply(.SD, \(x) NA), .SDcols = male_vars]

ukb_weights <- fread(file_paths[["ukb"]][["weight_file"]], colClasses = "character")[
  , .(id = f.eid, weight = as.numeric(LassoWeight))
]
ukb_weights <- merge.data.table(
  sub_pim[, .(id)],
  ukb_weights,
  by = "id"
)
windex <- which(sub_pim[, id] %in% ukb_weights[, id])

# covariates -------------------------------------------------------------------
if (opt$use_geno == TRUE) {
  x1 <- merged[, .(
    age, female,
    PC1, PC2, PC3, PC4)]
  x2 <- merged[, .(age, PC1, PC2, PC3, PC4)]
} else {
  x1 <- merged[, .(age, female)]
  x2 <- merged[, .(age)]
}

# MGI Partial Correlations -----------------------------------------------------
sub_pim <- sub_pim[, 2:ncol(sub_pim)]

cli_alert("calculating unweighted partial correlations....")
res_list <- ms::fcor(
  x        = sub_pim,
  covs     = x1,
  covs_alt = x2,
  n_cores  = n_cores,
  verbose  = TRUE
)

cli_alert_success("unweighted correlation complete! creating output table...")
if (is.data.table(res_list)) {
  res_table <- res_list
} else {
  res_table <- rbindlist(res_list, fill = TRUE)
}

cli_alert("calculating weighted partial correlations....")
wres_list <- ms::fcor(
  x        = sub_pim[windex, ],
  covs     = x1[windex, ],
  covs_alt = x2[windex, ],
  n_cores  = n_cores,
  w        = ukb_weights[["weight"]],
  verbose  = TRUE
)

cli_alert_success("weighted correlation complete! creating output table...")
if (is.data.table(res_list)) {
  wres_table <- wres_list
} else {
  wres_table <- rbindlist(wres_list, fill = TRUE)
}


# save results -----------------------------------------------------------------
output_file <- glue("data/private/ukb/{opt$ukb_version}/",
                    "ukb_phenome_partial_correlations_",
                    "{ifelse(opt$use_geno == TRUE, 'w_geno_pcs_', '')}",
                    "{opt$ukb_version}.qs")
message("saving results to: {.path {output_file}}")

save_qs(
  x    = res_table,
  file = output_file
)

woutput_file <- glue(
  "data/private/ukb/{opt$ukb_version}/",
  "ukb_phenome_weighted_partial_correlations_",
  "{ifelse(opt$use_geno == TRUE, 'w_geno_pcs_', '')}",
  "{opt$ukb_version}.qs"
)
message("saving results to: {.path {output_file}}")

save_qs(
  x    = wres_table,
  file = woutput_file
)

cli_alert_success("script success! see {.path {dirname(output_file)}}")
