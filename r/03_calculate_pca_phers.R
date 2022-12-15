# phers construction code
# an adaptaion from Lauren Beesley's D1_Construct PheRS PanCan.R script
# author: max salvatore
# date:   20221214

# libraries, functions, and options --------------------------------------------
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ppcor)
  library(ggplot2)
  library(gridExtra)
  library(MatchIt)
  library(glue)
  library(SPAtest)
  library(ggrepel)
  library(cli)
  library(optparse)
})

set.seed(61787)

source("fn/files-utils.R") # load partial correlation function
source(glue("https://raw.githubusercontent.com/umich-cphds/",
            "createUKBphenome/master/scripts/function.expandPhecodes.r"))

mytheme <- gridExtra::ttheme_default(
  core = list(fg_params = list(cex = 0.5)),
  colhead = list(fg_params = list(cex = 0.5)),
  rowhead = list(fg_params = list(cex = 0.5)))

# optparse list ----------------------------------------------------------------
option_list <- list(
  make_option("--outcome", type = "character", default = "",
              help = "Outcome phecode"),
  make_option("--mgi_version", type = "character", default = "20210318",
              help = "Version of MGI data [default = 20210318]"),
  make_option("--ukb_version", type = "character", default = "20221117",
              help = "Version of UKB data [default = 20221117]"),
  make_option("--time_threshold", type = "numeric", default = "0",
              help = glue("Time threshold for the phenome data ",
                          "[default = 0]")),
  make_option("--pc_var_explain", type = "numeric", default = "0.95",
              help = glue("Cumulative variation explained by PCs threshold ",
                          "[default = 0.95]"))
)

parser <- OptionParser(usage="%prog [options]", option_list = option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

# specifications ---------------------------------------------------------------
mgi_version       <- opt$mgi_version
ukb_version        <- opt$ukb_version
outcome           <- opt$outcome
pc_var_explain    <- opt$pc_var_explain
time_threshold    <- opt$time_threshold
derivation_cohort <- "mgi"

mgi_data_path <- glue("data/private/mgi/{mgi_version}/",
                      "X{gsub('X', '', outcome)}")
mgi_results_path <- glue("results/mgi/{mgi_version}/X{gsub('X', '', outcome)}/",
                         "pca")
if (!dir.exists(mgi_results_path)) {
  dir.create(mgi_results_path, recursive = TRUE)
}

## pull file paths corresponding to the data version specified
file_paths <- get_files(mgi_version = mgi_version, ukb_version = ukb_version)

# read in data -----------------------------------------------------------------
cli_alert_info("reading data...")

## mgi covariate data
mgi_cov <- fread(glue("{mgi_data_path}/matched_covariates.txt"))[
  get(glue("t{time_threshold}_indicator")) == 1, .(
    id, female,
    age_at_threshold = round(get(glue("t{time_threshold}_threshold")) /
                               365.25, 1),
    length_followup, case)]
mgi_cov <- mgi_cov[complete.cases(mgi_cov), ]

# mgi time-restricted pim
mgi_pim <- fread(glue("{mgi_data_path}/time_restricted_phenomes/",
                      "mgi_X{gsub('X', '', outcome)}_t{time_threshold}",
                      "_{mgi_version}.txt"))
mgi_pim <- merge(mgi_pim, data.table(id = mgi_cov$id),
                 by = "id", all.x = FALSE)
short_mgi_pim <- mgi_pim[, !c("id")]
short_mgi_pim[is.na(short_mgi_pim)] <- 0

## mgi master PIM - read 1 row for variable names only
mgi_pim0 <- fread(file_paths[["mgi"]]$pim0_file, nrow = 1)

## ukb master PIM - read 1 row for variable names only
ukb_pim0 <- fread(file_paths[["ukb"]]$pim0_file, nrow = 1)

## phecode info
pheinfo <- fread("data/public/Phecode_Definitions_FullTable_Modified.txt",
                 colClasses = "character")

# exclusions -------------------------------------------------------------------
outcome_vector <- short_mgi_pim[["case"]]

## exclusion range phecodes (from PhewasCatalog)
exclusionRange <- pheinfo[phecode == outcome, phecode_exclude_range]
exclusions1    <- pheinfo[phecode %in% unlist(unname(sapply(
  strsplit(exclusionRange,', {0,1}')[[1]],
  expandPhecodes))), phecode]

## phecodes not defined in both cohorts
exclusions2 <- pheinfo[!(phecode %in% gsub("X", "",
                                           intersect(names(mgi_pim0),
                                                     names(ukb_pim0)))),
                       phecode]

exclusionsX <- c("case", glue("X{c(gsub('X', '', outcome), ",
                              "union(exclusions1, exclusions2))}"))

# subset phenotype data --------------------------------------------------------
included <- names(short_mgi_pim)[!(names(short_mgi_pim) %in% exclusionsX)]
exclude_mgi_pim <- short_mgi_pim[, ..included]

# obtain principal components --------------------------------------------------
cli_alert_info("calculating principal components")
mgi_pca <- prcomp(exclude_mgi_pim, center = FALSE, scale. = FALSE)

pcs <- mgi_pca$x[
  , 1:which.min(abs(summary(mgi_pca)$importance[3, ] - pc_var_explain))]

maxs <- apply(abs(pcs), 2, max)
pcs_mod <- sweep(pcs, MARGIN = 2, maxs, "/")
pcs_mod <- pcs_mod + 1

rotations <- data.table(mgi_pca$rotation)
rotations <- rotations[
  , 1:which.min(abs(summary(mgi_pca)$importance[3, ] - pc_var_explain))]

# run SPAtest ------------------------------------------------------------------
cli_alert_info("running SPAtest")
## method 2: pca
test2 <- ScoreTest_SPA(
  genos = t(pcs_mod),
  pheno = outcome_vector,
  cov   = mgi_cov[id %in% mgi_pim[, id], !c("id")],
  method = "fastSPA",
  beta.out = TRUE, beta.Cutoff = 1)

# organize results -------------------------------------------------------------
cli_alert_info("organizing results and generating plots...")
results_pc <- data.frame(
  phecode      = glue("PC{c(1:length(pcs_mod[1, ]))}"),
  pvals_m2     = test2$p.value,
  alphas_m2    = test2$beta,
  alphasmod_m2 = sweep(as.matrix(test2$beta), MARGIN = 1, maxs, "/")
) |> as.data.table()

results_pc[, alphasmod_m2sig := fifelse(
  pvals_m2 > 0.05/length(results_pc[, 1]) | is.na(pvals_m2),
  rep(0, length(results_pc[, 1])),
  alphasmod_m2)]

results <- data.table(
  phecode = names(exclude_mgi_pim)
)

results$betas_m2 <- as.matrix(rotations) %*%
  as.matrix(results_pc[, alphasmod_m2])
results$betas_m2sig <- as.matrix(rotations) %*%
  as.matrix(results_pc[, alphasmod_m2sig])

### save rotations, results_pc, and results

# manhattan plot ---------------------------------------------------------------
results_short <- data.table(
  log10pvals = -log10(results_pc[, pvals_m2]),
  betas = results_pc[, alphasmod_m2]
)
results_short[, pch := fifelse(
  log10pvals > -log10(0.05/length(log10pvals)),
  fifelse(
    betas > 0,
    24,
    25
  ),
  21
)]

p <- results_short |>
  ggplot() +
  geom_hline(yintercept = -log10(0.05/length(results_short[, log10pvals])),
             linetype = 1, color = "red") +
  geom_point(aes(x = 1:nrow(results_short),
                 y = log10pvals),
             shape = results_short[, pch]) +
  labs(
    x = "Principal components",
    y = "-log10(p-value)",
    title = glue("Associations with phecode {outcome}: ",
                 "{pheinfo[phecode == gsub('X', '', outcome), description]}"),
    caption = glue("Cumulative variation explained by PCs: {pc_var_explain}; ",
                   "n PCs = {format(length(pcs_mod[1, ]), big.mark = ',')}; ",
                   "t = {time_threshold}; ",
                   "derivation cohort: {derivation_cohort}; ",
                   "application cohort: mgi")
  ) +
  theme_minimal() +
  theme(
    legend.position  = "",
    legend.text      = element_text(size = 10),
    legend.title     = element_blank(),
    axis.text.x      = element_text(angle = 60, hjust = 1, vjust = 1),
    plot.caption     = element_text(hjust = 0),
    text             = element_text(size = 12)
    )

ggsave(
  filename = glue("{mgi_results_path}/mgi_X{gsub('X', '', outcome)}_",
                  "t{time_threshold}_pve{pc_var_explain}_",
                  "pc_mahanttan.pdf"),
  plot     = p,
  width = 8, height = 6,
  device = cairo_pdf)

# beta plot --------------------------------------------------------------------
results <- merge.data.table(
  results,
  pheinfo[, .(
    phecode = glue_data(.SD, "X{phecode}"),
    color,
    desc = description
    )],
  by = "phecode"
)

p2 <- results |>
  ggplot() +
  geom_point(aes(
    x = 1:nrow(results),
    y = betas_m2sig,
    col = color
  ), shape = 16) +
  scale_colour_manual(values = unique(as.character(results[, color]))) +
  geom_label_repel(aes(x = c(1:nrow(results)), y = betas_m2sig, label = desc),
                   label.size = 0.1, force = 2, size = 2, 
                   label.padding = 0.1, point.padding = unit(0.2, "lines"),
                   segment.alpha  = 0.3) +
  labs(
    x = "Phenotypes",
    y = "Beta",
    title = glue("Associations with phecode {outcome}: ",
                 "{pheinfo[phecode == gsub('X', '', outcome), description]}"),
    caption = glue("Cumulative variation explained by PCs: {pc_var_explain}; ",
                   "n PCs = {format(length(pcs_mod[1, ]), big.mark = ',')}; ",
                   "t = {time_threshold}; ",
                   "derivation cohort: {derivation_cohort}; ",
                   "application cohort: mgi")
  ) +
  theme_minimal() +
  theme(
    plot.caption     = element_text(hjust = 0),
    legend.position  = "",
    legend.text      = element_text(size = 10),
    legend.title     = element_blank(),
    axis.text.x      = element_text(angle = 60, hjust = 1, vjust = 1),
    text = element_text(size = 12))

ggsave(
  filename = glue("{mgi_results_path}/mgi_X{gsub('X', '', outcome)}_",
                  "t{time_threshold}_pve{pc_var_explain}_",
                  "pc_betas.pdf"),
  plot     = p2,
  width = 8, height = 6,
  device = cairo_pdf)

# calculate phers --------------------------------------------------------------
cli_alert_info("calculating phers...")
mgi_phers_vec     <- as.matrix(exclude_mgi_pim) %*% as.matrix(results[, betas_m2sig])
mgi_phers_vec_std <- (mgi_phers - mean(mgi_phers)) / sd(mgi_phers)

mgi_phers <- data.frame(
  id        = mgi_pim[, id],
  phers     = mgi_phers_vec,
  phers_std = mgi_phers_ved_std
) |> as.data.table()

fwrite(
  x = mgi_phers,
  file = glue("{mgi_results_path}/mgi_X{gsub('X', '', outcome)}_",
              "t{time_threshold}_pve{pc_var_explain}_",
              "pc_phers.txt"),
  sep = "\t"
)

cli_alert_success("PCA-PheRS for {outcome} a success!")
