# calculate partial correlations
# outputs: table of pairwise partial correlations
# author:  max salvatore
# date:    20221212

# 1. libraries, functions, and options -----------------------------------------
options(stringsAsFactors=F)

library(data.table)
library(ppcor)
library(parallel)
library(pbmcapply)
library(glue)
library(cli)

set.seed(61787)

source("fn/quick_mgi_partial_correlation.R") # load partial correlation function

# 2. specifications ------------------------------------------------------------
use_geno    <- TRUE
mgi_version <- "20210318" # mgi data version
n_cores     <- parallel::detectCores()*0.5 # use 50% of available cores

file_paths <- get_files(mgi_version = mgi_version)

# 3. read data -----------------------------------------------------------------
cli::cli_alert_info("reading data...")
## phecode indicator matrix (PEDMASTER_0)
pim0 <- data.table::fread(file_paths[["mgi"]]$pim0_file)
data.table::setnames(pim0, old = "IID", new = "id")

## demographics data
demo <- data.table::fread(file_paths[["mgi"]]$demo_file)[, .(id = Deid_ID, age = Age, alive = as.numeric(AliveYN == "Y"), dead_dsb = Deceased_DaysSinceBirth, ethn = EthnicityName, marital = MaritalStatusCode, sex = Sex, race = RaceName)]

## pc data 
if (use_geno == TRUE) {
  pcs <- data.table::fread("/net/junglebook/magic_data/MGI_GenotypeData_Freeze4/MGI_Freeze4_Sample_Info_60215samples.txt")
  data.table::setnames(pcs, old = "Deid_ID", new = "id")
  # merge pcs in demo data
  merged <- merge(demo, pcs)
} else {
  merged <- demo
}

sub_pim <- pim0[id %in% merged[, id]]    # subset pim
merged  <- merged[id %in% sub_pim[, id]] # subset merged

# covariates -------------------------------------------------------------------
if (use_geno == TRUE) {
  X1_MGI = merged[, .(
    Age = age,
    FEMALE = as.numeric(sex == "F"),
    PC1, PC2, PC3, PC4)]
  X2_MGI = merged[, .(Age = age, PC1, PC2, PC3, PC4)]
} else {
  X1_MGI = merged[, .(Age = age, FEMALE = as.numeric(sex == "F"))]
  X2_MGI = merged[, .(Age = age)]
}

# MGI Partial Correlations -----------------------------------------------------
sub_pim <- sub_pim[, 2:ncol(sub_pim)]
cli::cli_alert_info("identifying pairwise combinations...")
combos  <- combn(names(sub_pim), 2, simplify = FALSE)

cli::cli_alert_info("calculating pairwise partial correlations....")
res_list <- pbmclapply(combos,
                       quick_mgi_partial_correlation,
                       mc.cores       = n_cores,
                       mc.preschedule = FALSE)

cli::cli_alert_success("calculation complete! creating output table...")
res_table <- rbindlist(res_list)

# save results -----------------------------------------------------------------
output_file <- glue::glue("data/private/mgi/{mgi_version}/mgi_phenome_partial_correlations_{ifelse(use_geno == TRUE, 'w_geno_pcs_', '')}{mgi_version}.txt")
cli::cli_alert_info("saving results to: {output_file}...")

data.table::fwrite(x    = res_table,
                   file = output_file,
                   sep  = "\t")
cli::cli_alert_success("file saved and script complete!")