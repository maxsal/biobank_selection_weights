# Conduct principal components analysis of the phecode indicator matrix in MGI
# author:  max salvatore
# date:    20220809

# 1. libraries, functions, and options -----------------------------------------
ms::libri(
  maxsal/ms, data.table, ggplot2, patchwork, qs, glue, optparse,
  survey, cli
)

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
  ),
  make_option("--mgi_weights",
    type = "character", default = "ip_selection",
    help = "Name of weight variable to use for MGI [default = %default]"
  )
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args   <- parse_args(parser, positional_arguments = 0)
opt    <- args$options
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
mgi_demo <- read_qs(glue("data/private/mgi/{opt$mgi_version}/datax_{opt$mgi_version}_comb.qs"))[, .(
  id  = DeID_PatientID,
  age = Age,
  sex = Sex
)]
mgi_demo <- mgi_demo[complete.cases(mgi_demo), ]

### phecode indicator matrix (PEDMASTER_0)
common_codes <- fread("data/public/phecodex_20plus.csv")[plus20 == 1, phecode]

mgi_pim0 <- qread(glue("data/private/mgi/{opt$mgi_version}/MGI_PIM0X_{opt$mgi_version}.qs"))
pim_vars <- c("id", common_codes)
mgi_pim0 <- mgi_pim0[, ..pim_vars]
mgi_pim0 <- merge(
  mgi_pim0,
  data.frame(id = mgi_demo[, id]),
  by = "id",
  all.x = FALSE
)

mgi_weights <- read_qs(glue("data/private/mgi/{opt$mgi_version}/weights_{opt$mgi_version}_comb.qs"))
setnames(mgi_weights, old = opt$mgi_weights, new = "weight")
mgi_weights <- mgi_weights[!is.na(weight), ]
weight_vars <- c("id", "weight")
mgi_weights <- merge.data.table(
  mgi_pim0[, .(id)],
  mgi_weights[, ..weight_vars],
  by = "id"
)

short_mgi_pim <- mgi_pim0[, !c("id")]
short_mgi_pim[is.na(short_mgi_pim)] <- 0

wshort_mgi_pim <- merge.data.table(
  mgi_pim0[id %in% unique(mgi_weights$id), ],
  mgi_weights,
  by.x = "id", by.y = "id"
)[, !c("id")]
wshort_mgi_pim[is.na(wshort_mgi_pim)] <- 0

mgi_design <- svydesign(
  id      = ~1,
  weights = ~weight,
  data    = wshort_mgi_pim
)

## ukb
### demographics
ukb_demo <- read_qs(glue("data/private/ukb/{opt$ukb_version}/datax_{opt$ukb_version}_comb.qs"))[in_phenome == 1, .(
  id   = as.character(id),
  age  = age_at_consent,
  ethn = race_eth,
  sex 
)]
ukb_demo <- ukb_demo[complete.cases(ukb_demo), ]

### phecode indicator matrix (PEDMASTER_0)
ukb_pim0 <- read_qs(glue("/net/junglebook/home/mmsalva/projects/dissertation/aim_one/data/private/ukb/{opt$ukb_version}/UKB_PIM0X_{opt$ukb_version}.qs"))
ukb_pim0 <- merge(
  ukb_pim0,
  ukb_demo[, .(id)],
  by = "id",
  all.x = FALSE
)
ukb_pim0 <- ukb_pim0[, ..pim_vars]

ukb_weights <- fread(file_paths[["ukb"]][["weight_file"]], colClasses = "character")[
  , .(id = f.eid, weight = as.numeric(LassoWeight))
]
ukb_weights <- ukb_weights[!is.na(weight), ]

short_ukb_pim <- ukb_pim0[, !c("id")]
short_ukb_pim[is.na(short_ukb_pim)] <- 0

wshort_ukb_pim <- merge.data.table(
  ukb_pim0[id %in% unique(ukb_weights$id), ],
  ukb_weights,
  by = "id"
)[, !c("id")]
wshort_ukb_pim[is.na(wshort_ukb_pim)] <- 0

ukb_design <- svydesign(
  id      = ~1,
  weights = ~weight,
  data    = wshort_ukb_pim
)

# 4. calculates pcs ------------------------------------------------------------
## mgi
### calculate
mgi_pca     <- prcomp(short_mgi_pim, center = FALSE, scale. = FALSE)
mgi_pca_sum <- summary(mgi_pca)$importance

mgi_vars <- setdiff(names(wshort_mgi_pim), "weight")
mgi_wpca <- svyprcomp(
  as.formula(paste0("~", paste0(mgi_vars, collapse = " + "))),
  design = mgi_design, scale. = FALSE, center = FALSE
)
mgi_wpca_sum <- summary(mgi_wpca)$importance

### save data
save_qs(
  x = as.data.table(mgi_pca_sum),
  file = glue("results/mgi/{opt$mgi_version}/mgi_pcax_importance_{opt$mgi_version}.qs")
)
save_qs(
  x = as.data.table(mgi_wpca_sum),
  file = glue("results/mgi/{opt$mgi_version}/mgi_wpcax_importance_{opt$mgi_version}.qs")
)
### save plot
mgi_pca_plot <- stacked_pca_plot(x = mgi_pca_sum, cohort = "mgi")
ggsave(
  plot = mgi_pca_plot,
  filename = glue(
    "results/mgi/{opt$mgi_version}/",
    "stacked_pcax_plot_{opt$mgi_version}.pdf"
  ),
  device = cairo_pdf,
  width = 6,
  height = 6
)
mgi_wpca_plot <- stacked_pca_plot(x = mgi_wpca_sum, cohort = "mgi")
ggsave(
  plot = mgi_wpca_plot,
  filename = glue(
    "results/mgi/{opt$mgi_version}/",
    "stacked_wpcax_plot_{opt$mgi_version}.pdf"
  ),
  device = cairo_pdf,
  width = 6,
  height = 6
)

## ukb
### calculate
ukb_pca     <- prcomp(short_ukb_pim, center = FALSE, scale. = FALSE)
ukb_pca_sum <- summary(ukb_pca)$importance

ukb_vars <- setdiff(names(wshort_ukb_pim), "weight")
ukb_wpca <- svyprcomp(
  as.formula(paste0("~", paste0(ukb_vars, collapse = " + "))),
  design = ukb_design, scale. = FALSE, center = FALSE
)
ukb_wpca_sum <- summary(ukb_wpca)$importance

### save data
save_qs(
  x = as.data.table(ukb_pca_sum),
  file = glue("results/ukb/{opt$ukb_version}/ukb_pcax_importance_{opt$ukb_version}.qs")
)
save_qs(
  x = as.data.table(ukb_wpca_sum),
  file = glue("results/ukb/{opt$ukb_version}/ukb_wpcax_importance_{opt$ukb_version}.qs")
)
### save plot
ukb_pca_plot <- stacked_pca_plot(x = ukb_pca_sum, cohort = "ukb")
ggsave(
  plot = ukb_pca_plot,
  filename = glue(
    "results/ukb/{opt$ukb_version}/",
    "stacked_pcax_plot_{opt$ukb_version}.pdf"
  ),
  device = cairo_pdf,
  width = 6,
  height = 6
)

ukb_wpca_plot <- stacked_pca_plot(x = ukb_wpca_sum, cohort = "ukb")
ggsave(
  plot = ukb_wpca_plot,
  filename = glue(
    "results/ukb/{opt$ukb_version}/",
    "stacked_wpcax_plot_{opt$ukb_version}.pdf"
  ),
  device = cairo_pdf,
  width = 6,
  height = 6
)

cli_alert_success("script complete! 🎉")
