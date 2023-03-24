# calculate PCs of phenome data
# outputs: pcs
# author:  max salvatore
# date:    20221208

# 1. libraries, functions, and options -----------------------------------------
library(data.table)
library(ggplot2)
library(patchwork)
library(glue)
library(qs)

set.seed(61787)

for (i in list.files("fn/")) source(paste0("fn/", i)) # load functions

# 2. specifications ------------------------------------------------------------
mgi_version     <- "20210318"       # mgi version
ukb_version     <- "20221117"       # ukb  version

## extract file paths
file_paths <- get_files(mgi_version = mgi_version, ukb_version = ukb_version)

# 3. read data -----------------------------------------------------------------
## mgi
### demographics
mgi_demo <- fread(file_paths[["mgi"]]$demo_file)[, .(
  id       = Deid_ID,
  age      = Age,
  sex      = Sex,
  race     = RaceName)][, .(id, age, sex)]
mgi_demo <- mgi_demo[complete.cases(mgi_demo), ]

### phecode indicator matrix (PEDMASTER_0)
mgi_pim0 <- fread(file_paths[["mgi"]]$pim0_file)
mgi_pim0 <- merge(
  mgi_pim0,
  data.frame(IID = mgi_demo[, id]),
  by = "IID",
  all.x = FALSE)
short_mgi_pim <- mgi_pim0[, !c("IID")]
short_mgi_pim[is.na(short_mgi_pim)] <- 0

## ukb
### demographics
ukb_demo <- fread(file_paths[["ukb"]]$demo_file,
                              na.strings = c("", "NA", "."),
                              colClass = "character")
ukb_demo <- ukb_demo[, .(
  id   = as.character(id),
  dob  = as.Date(dob),
  age  = floor(as.numeric(as.Date("2022-11-17") - as.Date(dob)) / 365.25),
  ethn = ethnicity,
  sex)][, .(id, age, sex)]
ukb_demo <- ukb_demo[complete.cases(ukb_demo), ]

### phecode indicator matrix (PEDMASTER_0)
ukb_pim0 <- fread(file_paths[["ukb"]]$pim0_file)
ukb_pim0[, IID := as.character(IID)]
ukb_pim0 <- merge(
  ukb_pim0,
  data.frame(IID = ukb_demo[, id]),
  by = "IID",
  all.x = FALSE)
short_ukb_pim <- ukb_pim0[, !c("IID")]
short_ukb_pim[is.na(short_ukb_pim)] <- 0

# 4. calculates pcs ------------------------------------------------------------
## mgi
### calculate
mgi_pca      <- prcomp(short_mgi_pim, center = FALSE, scale. = FALSE)
mgi_pca_sum  <- summary(mgi_pca)$importance
### save data
save_qs(
  x = as.data.table(mgi_pca_sum),
  file = glue("results/mgi/{mgi_version}/mgi_pca_importance_{mgi_version}.qs"))
### save plot
mgi_pca_plot <- stacked_pca_plot(x = mgi_pca_sum, cohort = "mgi")
ggsave(
  plot = mgi_pca_plot,
  filename = glue("results/mgi/{mgi_version}/",
                        "stacked_pca_plot_{mgi_version}.pdf"),
  device = cairo_pdf,
  width = 6,
  height = 6
)

## ukb
### calculate
ukb_pca      <- prcomp(short_ukb_pim, center = FALSE, scale. = FALSE)
ukb_pca_sum  <- summary(ukb_pca)$importance
### save data
save_qs(
  x = as.data.table(ukb_pca_sum),
  file = glue("results/ukb/{ukb_version}/ukb_pca_importance_{ukb_version}.qs"))
### save plot
ukb_pca_plot <- stacked_pca_plot(x = ukb_pca_sum, cohort = "ukb")
ggsave(
  plot = ukb_pca_plot,
  filename = glue("results/ukb/{ukb_version}/",
                        "stacked_pca_plot_{ukb_version}.pdf"),
  device = cairo_pdf,
  width = 6,
  height = 6
)
