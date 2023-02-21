# calculate partial correlations
# outputs: table of pairwise partial correlations
# author:  max salvatore
# date:    20230220

# 1. libraries, functions, and options -----------------------------------------
options(stringsAsFactors = FALSE)

library(data.table)
library(fst)
library(ppcor)
library(glue)
library(cli)
library(doMC)
library(optparse)
library(progressr)

set.seed(61787)

source("fn/files-utils.R") # load partial correlation function
source("fn/partial_corr_veloce.R") # load partial correlation function

# optparse list ----------------------------------------------------------------
option_list <- list(
  make_option("--use_geno", type = "logical", default = TRUE,
              help = "Adjust for genotype PCs [default = %default]"),
  make_option("--ukb_version", type = "character", default = "20221117",
              help = "Version of MGI data [default = %default]")
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args   <- parse_args(parser, positional_arguments = 0)
opt    <- args$options
print(opt)

file_paths <- get_files(ukb_version = opt$ukb_version)

handlers(global = TRUE)

# 3. read data -----------------------------------------------------------------
cli_alert("reading data...")
## phecode indicator matrix (PEDMASTER_0)
pim0 <- fread(file_paths[["ukb"]][["pim0_file"]])
setnames(pim0, old = "IID", new = "id")

## demographics data
icd_phecode <- fread(file_paths[["ukb"]][["icd_phecode_file"]])
icd_phecode <- icd_phecode[ icd_phecode[, .I[which.max(dsb)], by = "id"][["V1"]] ][, .(id, dsb)][, age_at_last_first := round(dsb / 365.25, 3)][]

demo <- fread(file_paths[["ukb"]][["demo_file"]])[, `:=` (
  age_today = as.numeric(round((as.Date("2022-11-17") - as.Date(dob))/365.25, 3))
  )][id %in% icd_phecode[, id]][]

merge.data.table(
  unique(demo[, .(id, age_today)][!is.na(age_today), ]),
  icd_phecode[, .(id, age_at_last_first)],
  by = "id"
)

## pc data
if (opt$use_geno == TRUE) {
  pcs <- fread(glue("/net/junglebook/magic_data/MGI_GenotypeData_Freeze4/",
                    "MGI_Freeze4_Sample_Info_60215samples.txt"))
  setnames(pcs, old = "Deid_ID", new = "id")
  # merge pcs in demo data
  merged <- merge.data.table(demo, pcs)
} else {
  merged <- demo
}

sub_pim <- pim0[id %in% merged[, id]]    # subset pim
merged  <- merged[id %in% sub_pim[, id]] # subset merged

# covariates -------------------------------------------------------------------
if (opt$use_geno == TRUE) {
  x1_mgi <- merged[, .(
    Age = age,
    FEMALE = as.numeric(sex == "F"),
    PC1, PC2, PC3, PC4)]
  x2_mgi <- merged[, .(Age = age, PC1, PC2, PC3, PC4)]
} else {
  x1_mgi <- merged[, .(Age = age, FEMALE = as.numeric(sex == "F"))]
  x2_mgi <- merged[, .(Age = age)]
}

# MGI Partial Correlations -----------------------------------------------------
sub_pim <- sub_pim[, 2:ncol(sub_pim)]
cli_alert("identifying pairwise combinations...")
combos  <- combn(names(sub_pim), 2, simplify = FALSE)

cli_alert("calculating pairwise partial correlations....")
res_list <- partial_corr_veloce(
  pim   = sub_pim,
  ncore = detectCores()/2,
  covs1 = x1_mgi,
  covs2 = x2_mgi
)

cli_alert_success("calculation complete! creating output table...")
res_table <- rbindlist(res_list, fill = TRUE)

# save results -----------------------------------------------------------------
output_file <- glue("data/private/mgi/{opt$mgi_version}/",
                    "mgi_phenome_partial_correlations_",
                    "{ifelse(opt$use_geno == TRUE, 'w_geno_pcs_', '')}",
                    "{opt$mgi_version}.fst")
cli_alert_info("saving results to: {.path {output_file}}")

write_fst(x    = res_table,
          path = output_file)

cli_alert_success("script success! see {.path {output_file}}")
