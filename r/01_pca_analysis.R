# Conduct principal components analysis of the phecode indicator matrix in MGI
# author:  max salvatore
# date:    20220418

# 1. libraries, functions, and options -----------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(qs)
  library(glue)
  library(optparse)
})
options(datatable.print.class = TRUE)

set.seed(61787)

for (i in list.files("fn/", full.names = TRUE)) source(i)

# optparse list ----
option_list <- list(
  make_option("--mgi_version",
    type = "character", default = "20220822",
    help = "MGI cohort version in /net/junglebook/magic_data/EHRdata/ [default = %default]"
  ),
  make_option("--ukb_version",
    type = "character", default = "20221117",
    help = "UKB cohort version in /net/junglebook/magic_data/EHRdata/ [default = %default]"
  )
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

# 2. specifications ------------------------------------------------------------

## extract file paths
file_paths <- get_files(
  mgi_version = opt$mgi_version,
  ukb_version = opt$ukb_version
)

# 3. read data -----------------------------------------------------------------
## mgi
### demographics
mgi_demo <- read_qs(file_paths[["mgi"]][["cov_processed_file"]])[, .(
  id  = DeID_PatientID,
  age = Age,
  sex = Sex
)]
mgi_demo[complete.cases(mgi_demo), ]

### phecode indicator matrix (PEDMASTER_0)
mgi_pim0 <- fread(file_paths[["mgi"]][["pim0_file"]])
mgi_pim0 <- merge(
  mgi_pim0,
  data.frame(IID = mgi_demo[, id]),
  by = "IID",
  all.x = FALSE
)
short_mgi_pim <- mgi_pim0[, !c("IID")]
short_mgi_pim[is.na(short_mgi_pim)] <- 0

## ukb
### demographics
ukb_demo <- read_qs(file_paths[["ukb"]][["demo_file"]])[in_phenome == 1, .(
  id   = as.character(id),
  age  = age_at_consent,
  ethn = race_eth,
  sex 
)]
ukb_demo <- ukb_demo[complete.cases(ukb_demo), ]

### phecode indicator matrix (PEDMASTER_0)
ukb_pim0 <- read_qs(file_paths[["ukb"]][["pim0_file"]])
ukb_pim0 <- merge(
  ukb_pim0,
  ukb_demo[, .(id)],
  by = "id",
  all.x = FALSE
)
short_ukb_pim <- ukb_pim0[, !c("id")]
short_ukb_pim[is.na(short_ukb_pim)] <- 0

# 4. calculates pcs ------------------------------------------------------------
## mgi
### calculate
mgi_pca     <- prcomp(short_mgi_pim, center = FALSE, scale. = FALSE)
mgi_pca_sum <- summary(mgi_pca)$importance
### save data
save_qs(
  x = as.data.table(mgi_pca_sum),
  file = glue("results/mgi/{opt$mgi_version}/mgi_pca_importance_{opt$mgi_version}.qs")
)
### save plot
mgi_pca_plot <- stacked_pca_plot(x = mgi_pca_sum, cohort = "mgi")
ggsave(
  plot = mgi_pca_plot,
  filename = glue(
    "results/mgi/{opt$mgi_version}/",
    "stacked_pca_plot_{opt$mgi_version}.pdf"
  ),
  device = cairo_pdf,
  width = 6,
  height = 6
)

## ukb
### calculate
ukb_pca <- prcomp(short_ukb_pim, center = FALSE, scale. = FALSE)
ukb_pca_sum <- summary(ukb_pca)$importance
### save data
save_qs(
  x = as.data.table(ukb_pca_sum),
  file = glue("results/ukb/{opt$ukb_version}/ukb_pca_importance_{opt$ukb_version}.qs")
)
### save plot
ukb_pca_plot <- stacked_pca_plot(x = ukb_pca_sum, cohort = "ukb")
ggsave(
  plot = ukb_pca_plot,
  filename = glue(
    "results/ukb/{opt$ukb_version}/",
    "stacked_pca_plot_{opt$ukb_version}.pdf"
  ),
  device = cairo_pdf,
  width = 6,
  height = 6
)
