# Calculate partial correlations in MGI
# author:  max salvatore
# date:    20230418

# 1. libraries, functions, and options -----------------------------------------
options(stringsAsFactors = FALSE)

ms::libri(
  ms, data.table, qs, glue, cli, doMC, optparse, parallelly, foreach, tictoc
)

set.seed(61787)

source("fn/files-utils.R")

# optparse list ----------------------------------------------------------------
option_list <- list(
  make_option("--use_geno",
    type = "logical", default = FALSE,
    help = "Adjust for genotype PCs [default = %default]"
  ),
  make_option("--mgi_version",
    type = "character", default = "20220822",
    help = "Version of MGI data [default = %default]"
  ),
  make_option("--weights",
    type = "character", default = "ip_selection",
    help = glue(
      "Weighting variable to use for weighted analyses - ",
      "selection, all, or list of named weight variables [default = %default]"
    )
  ),
  make_option("--core_prop",
    type = "numeric", default = "0.125",
    help = "Proportion of available cores to use for parallel work [default = %default]"
  )
)

parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

file_paths <- get_files(mgi_version = opt$mgi_version)
n_cores    <- availableCores() * opt$core_prop

# 3. read data -----------------------------------------------------------------
cli_alert("reading data...")
## phecode indicator matrix (PEDMASTER_0)
pim0 <- qread(glue("data/private/mgi/{opt$mgi_version}/MGI_PIM0X_{opt$mgi_version}.qs"))
if ("IID" %in% names(pim0)) setnames(pim0, old = "IID", new = "id")

# restrict to phecodes with at least 20 occurrences in all cohorts
common_codes <- fread("data/public/phecodex_20plus.csv")[plus20 == 1, phecode]
restrict_vars <- c("id", common_codes)
pim0 <- pim0[, ..restrict_vars]

## demographics data
demo <- qread(glue("data/private/mgi/{opt$mgi_version}/datax_{opt$mgi_version}_comb.qs"))[, .(
  id       = DeID_PatientID,
  age      = Age,
  alive    = as.numeric(AliveYN == "Y"),
  dead_dsb = Deceased_DaysSinceBirth,
  ethn     = Ethnicity,
  marital  = MaritalStatusCode,
  sex      = Sex,
  race     = Race
)]

## pc data
if (opt$use_geno == TRUE) {
  pcs <- fread(glue(
    "/net/junglebook/magic_data/MGI_GenotypeData_Freeze4/",
    "MGI_Freeze4_Sample_Info_60215samples.txt"
  ))
  setnames(pcs, old = "Deid_ID", new = "id")
  # merge pcs in demo data
  merged <- merge.data.table(demo, pcs)
} else {
  merged <- demo
}

sub_pim <- pim0[id %in% merged[, id]] # subset pim
merged  <- merged[id %in% sub_pim[, id]] # subset merged

# sex-specific
male_ids   <- demo[sex == "M", id]
female_ids <- demo[sex == "F", id]

male_vars   <- intersect(ms::pheinfox[sex == "Male", phecode], names(sub_pim))
female_vars <- intersect(ms::pheinfox[sex == "Female", phecode], names(sub_pim))

sub_pim[id %in% male_ids, (female_vars) := lapply(.SD, \(x) NA), .SDcols = female_vars]
sub_pim[id %in% female_ids, (male_vars) := lapply(.SD, \(x) NA), .SDcols = male_vars]

mgi_weights <- read_qs(glue("data/private/mgi/{opt$mgi_version}/weightsx_{opt$mgi_version}_comb.qs"))
weight_vars <- c("id", opt$weights)
mgi_weights <- mgi_weights[!is.na(get(opt$weights)), ..weight_vars]
mgi_weights <- merge.data.table(
  sub_pim[, .(id)],
  mgi_weights[, ..weight_vars],
  by = "id"
)
windex <- which(sub_pim[, id] %in% mgi_weights[, id])

# covariates -------------------------------------------------------------------
if (opt$use_geno == TRUE) {
  x1_mgi <- merged[, .(
    Age    = age,
    FEMALE = as.numeric(sex == "F"),
    PC1, PC2, PC3, PC4
  )]
  x2_mgi <- merged[, .(Age = age, PC1, PC2, PC3, PC4)]
} else {
  x1_mgi <- merged[, .(Age = age, FEMALE = as.numeric(sex == "F"))]
  x2_mgi <- merged[, .(Age = age)]
}

# MGI Partial Correlations -----------------------------------------------------
sub_pim <- sub_pim[, 2:ncol(sub_pim)]

cli_alert("calculating unweighted partial correlations....")
tic()
res_list <- ms::fcor(
  x        = sub_pim,
  covs     = x1_mgi,
  covs_alt = x2_mgi,
  n_cores  = n_cores,
  verbose  = TRUE
)
toc()
cli_alert_success("unweighted correlation complete! creating output table...")
if (is.data.table(res_list)) {
  res_table <- res_list
} else {
  res_table <- rbindlist(res_list, fill = TRUE)
}

output_file <- glue(
  "data/private/mgi/{opt$mgi_version}/",
  "mgi_phenomex_partial_correlations_",
  "{ifelse(opt$use_geno == TRUE, 'w_geno_pcs_', '')}",
  "{opt$mgi_version}.qs"
)

save_qs(
  x    = res_table,
  file = output_file
)


cli_alert("calculating weighted partial correlations....")
wres_list <- ms::fcor(
  x        = sub_pim[windex, ],
  covs     = x1_mgi[windex, ],
  covs_alt = x2_mgi[windex, ],
  n_cores  = n_cores,
  w        = mgi_weights[[opt$weights]],
  verbose  = TRUE
)

cli_alert_success("weighted correlation complete! creating output table...")
if (is.data.table(res_list)) {
  wres_table <- wres_list
} else {
  wres_table <- rbindlist(wres_list, fill = TRUE)
}

woutput_file <- glue(
  "data/private/mgi/{opt$mgi_version}/",
  "mgi_phenomex_weighted_partial_correlations_",
  "{ifelse(opt$use_geno == TRUE, 'w_geno_pcs_', '')}",
  "{opt$mgi_version}.qs"
)

save_qs(
  x    = wres_table,
  file = woutput_file
)

cli_alert_success("script success! see {.path {dirname(output_file)}}")
