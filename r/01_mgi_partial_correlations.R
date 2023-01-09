# calculate partial correlations
# outputs: table of pairwise partial correlations
# author:  max salvatore
# date:    20230109

# 1. libraries, functions, and options -----------------------------------------
options(stringsAsFactors = FALSE)

library(data.table)
library(ppcor)
library(glue)
library(cli)
library(doMC)

set.seed(61787)

source("fn/files-utils.R") # load partial correlation function
source("fn/partial_corr_veloce.R") # load partial correlation function

# 2. specifications ------------------------------------------------------------
use_geno    <- TRUE
mgi_version <- "20210318" # mgi data version
parallelize <- FALSE

file_paths <- get_files(mgi_version = mgi_version)

# 3. read data -----------------------------------------------------------------
cli_alert_info("reading data...")
## phecode indicator matrix (PEDMASTER_0)
pim0 <- fread(file_paths[["mgi"]]$pim0_file)
setnames(pim0, old = "IID", new = "id")

## demographics data
demo <- fread(file_paths[["mgi"]]$demo_file)[, .(
  id       = Deid_ID,
  age      = Age,
  alive    = as.numeric(AliveYN == "Y"),
  dead_dsb = Deceased_DaysSinceBirth,
  ethn     = EthnicityName,
  marital  = MaritalStatusCode,
  sex      = Sex,
  race     = RaceName
  )]

## pc data
if (use_geno == TRUE) {
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
if (use_geno == TRUE) {
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
cli_alert_info("identifying pairwise combinations...")
combos  <- combn(names(sub_pim), 2, simplify = FALSE)

cli_alert_info("calculating pairwise partial correlations....")
res_list <- partial_corr_veloce(
  pim   = sub_pim,
  ncore = detectCores()/2,
  covs1 = x1_mgi,
  covs2 = x2_mgi
)

cli_alert_success("calculation complete! creating output table...")
res_table <- rbindlist(res_list, fill = TRUE)

# save results -----------------------------------------------------------------
output_file <- glue("data/private/mgi/{mgi_version}/",
                    "mgi_phenome_partial_correlations_",
                    "{ifelse(use_geno == TRUE, 'w_geno_pcs_', '')}",
                    "{mgi_version}.txt")
cli_alert_info("saving results to: {output_file}...")

fwrite(x    = res_table,
       file = output_file,
       sep  = "\t")
cli_alert_success("file saved and script complete!")
