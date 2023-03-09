# construct phers
# author: max salvatore
# date:   20230118

suppressPackageStartupMessages({
  library(cli)
})

# 1. libraries, functions, and options (outcome agnostic) ----------------------
cli_alert("loading packages and initializing...")
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
  library(glue)
  library(logistf)
  library(fst)
  library(optparse)
  library(colorblindr)
})

set.seed(61787)

lapply(list.files("fn/", full.names = TRUE), source) |> # load functions
  invisible()

# optparse list ----------------------------------------------------------------
option_list <- list(
  make_option("--outcome", type = "character", default = "157",
              help = "Outcome phecode [default = %default]"),
  make_option("--mgi_version", type = "character", default = "20220822",
              help = "Version of MGI data [default = %default]"),
  make_option("--ukb_version", type = "character", default = "20221117",
              help = "Version of UKB data [default = %default]"),
  make_option("--time_threshold", type = "numeric", default = "0",
              help = glue("Time threshold for the phenome data ",
                          "[default = %default]")),
  make_option("--discovery_cohort", type = "character", default = "mgi",
              help = glue("Use co-occurrence results from discovery cohort in phers ",
                          "[default = %default]")),
  make_option("--method", type = "character", default = "pwide_sig",
              help = glue("Method for determining phecodes for PheRS ('pwide_sig' or 'tophits') ",
                          "[default = %default]")),
  make_option("--tophits_n", type = "numeric", default = "50",
              help = glue("Number of top hits to use in top hits PheRS ",
                          "[default = %default]")),
  make_option("--weights", type = "character", default = NULL,
              help = glue("Name of weight suffix for weighted cooccurrence ",
                          "[default = %default]"))
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args   <- parse_args(parser, positional_arguments = 0)
opt    <- args$options
print(opt)

# check output folder exists ---------------------------------------------------
out_path <- glue("results/{coh}/{coh_version}/X{outc}/naive/",
                 coh = opt$discovery_cohort,
                 coh_version = ifelse(opt$discovery_cohort == "mgi",
                                      opt$mgi_version,
                                      ifelse(opt$discovery_cohort == "ukb",
                                             opt$ukb_version, NA)),
                 outc = gsub("X", "", opt$outcome))
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

# 2. specifications (specifies outcome) ----------------------------------------
external_cohort <- ifelse(opt$discovery_cohort == "mgi", "ukb", "mgi")
w               <- ifelse(is.null(opt$weights), "naive", opt$weights)

mgi_out_prefix <- glue("{opt$discovery_cohort}d_mgi_X{gsub('X', '', opt$outcome)}_t{opt$time_threshold}_{opt$method}_{w}_")
ukb_out_prefix <- glue("{opt$discovery_cohort}d_ukb_X{gsub('X', '', opt$outcome)}_t{opt$time_threshold}_{opt$method}_{w}_")
comb_out_prefix <- glue("{opt$discovery_cohort}d_{external_cohort}e_X{gsub('X', '', opt$outcome)}_t{opt$time_threshold}_{opt$method}_{w}_")

## extract file paths
file_paths <- get_files(mgi_version = opt$mgi_version,
                        ukb_version = opt$ukb_version)

# 3. read data -----------------------------------------------------------------
cli_alert_info("loading data...")
## mgi
mgi_pim0    <- fread(file_paths[["mgi"]]$pim0_file)
mgi_pim     <- read_fst(glue("data/private/mgi/{opt$mgi_version}/",
                          "X{gsub('X', '', opt$outcome)}/time_restricted_phenomes/",
                          "mgi_X{gsub('X', '', opt$outcome)}_",
                          "t{opt$time_threshold}_{opt$mgi_version}.fst"),
                     as.data.table = TRUE)
mgi_pim[is.na(mgi_pim)] <- 0
mgi_covariates <- read_fst(glue("data/private/mgi/{opt$mgi_version}/",
                             "X{gsub('X', '', opt$outcome)}/",
                             "matched_covariates.fst"),
                           as.data.table = TRUE)

## ukb
ukb_pim0    <- fread(file_paths[["ukb"]]$pim0_file)
ukb_pim     <- read_fst(glue("data/private/ukb/{opt$ukb_version}/",
                          "X{gsub('X', '', opt$outcome)}/",
                          "time_restricted_phenomes/",
                          "ukb_X{gsub('X', '', opt$outcome)}_",
                          "t{opt$time_threshold}_{opt$ukb_version}.fst"),
                        as.data.table = TRUE)
ukb_covariates <- read_fst(glue("data/private/ukb/{opt$ukb_version}/",
                             "X{gsub('X', '', opt$outcome)}/",
                             "matched_covariates.fst"),
                           as.data.table = TRUE)

## cooccur
if (opt$discovery_cohort == "mgi") {
  cooccur <- read_fst(glue("results/mgi/{opt$mgi_version}/",
                           "X{gsub('X', '', opt$outcome)}/",
                           "mgi_X{gsub('X', '', opt$outcome)}_",
                           "t{opt$time_threshold}_{opt$mgi_version}_",
                           "{ifelse(is.null(opt$weights), '', paste0(opt$weights, '_'))}results.fst"),
                      as.data.table = TRUE)
} else {
  cooccur <- read_fst(glue("results/ukb/{opt$ukb_version}/",
                           "X{gsub('X', '', opt$outcome)}/",
                           "ukb_X{gsub('X', '', opt$outcome)}_",
                           "t{opt$time_threshold}_{opt$ukb_version}",
                           "{ifelse(is.null(opt$weights), '', paste0(opt$weights, '_'))}_results.fst"),
                      as.data.table = TRUE)
}

## phecode info
pheinfo <- fread("data/public/Phecode_Definitions_FullTable_Modified.txt",
                 colClasses = "character")

# 4. phecode exclusions ------------------------------------------------------
## exclusion range from PhewasCatalog
exclusionRange <- pheinfo[phecode == opt$outcome, phecode_exclude_range]
exclusions1    <- pheinfo[phecode %in% unlist(unname(sapply(
  strsplit(exclusionRange,', {0,1}')[[1]],
  expandPhecodes))), phecode]

## phecodes not defined in both cohorts
exclusions2 <- pheinfo[!(phecode %in% gsub("X", "",
                                           intersect(names(mgi_pim0),
                                                     names(ukb_pim0)))),
                       phecode]

exclusionsX <- glue("X{union(exclusions1, exclusions2)}")

## helper functions
pretty_round <- function(x, r) {
  format(round(x, r), big.mark = ",", nsmall = r)
}
pretty_print <- function(x, r = 3) {
  paste0( pretty_round(x[2], r), " (", pretty_round(x[1], r), ", ",
          pretty_round(x[3], r), ")")
}
extractr_or <- function(x, r = 2) {
  suppressMessages(y <- confint(x))
  data.table(
    or_est = exp(coef(x)[["phers"]]),
    or_lo = exp(y["phers", 1]),
    or_hi = exp(y["phers", 2])
  )[, print := paste0(format(round(or_est, r), big.mark = ",", nsmall = r), " (",
                      format(round(or_lo, r), big.mark = ",", nsmall = r), ", ",
                      format(round(or_hi, r), big.mark = ",", nsmall = r), ")")][]
}

# 5. calculate naive phers -----------------------------------------------------
cli_alert_info("calculating naive phers for phecode {gsub('X', '', opt$outcome)}...")
## naive
### using MGI data
#### mgi
mgi_phers <- calculate_phers(
  pim        = mgi_pim,
  res        = cooccur[!(phecode %in% exclusionsX), ],
  method     = opt$method,
  tophits_n  = opt$tophits_n
)


#### ukb
ukb_phers <- calculate_phers(
  pim        = ukb_pim,
  res        = cooccur[!(phecode %in% exclusionsX), ],
  method     = opt$method,
  tophits_n  = opt$tophits_n
)

suppressMessages({
  mgi_roc <- pROC::roc(mgi_phers[["data"]][, case], mgi_phers[["data"]][, phers])
  mgi_auc <- pROC::ci.auc(mgi_phers[["data"]][, case], mgi_phers[["data"]][, phers])
  
  ukb_roc <- pROC::roc(ukb_phers[["data"]][, case], ukb_phers[["data"]][, phers])
  ukb_auc <- pROC::ci.auc(ukb_phers[["data"]][, case], ukb_phers[["data"]][, phers])
})


mgi_stuff <- data.table(
  sensitivity = mgi_roc$sensitivities,
  specificity = mgi_roc$specificities,
  test_data   = glue("MGI ({ifelse(opt$discovery_cohort == 'mgi', 'discovery', 'external')})")
)
ukb_stuff <- data.table(
  sensitivity = ukb_roc$sensitivities,
  specificity = ukb_roc$specificities,
  test_data   = glue("UKB ({ifelse(opt$discovery_cohort == 'mgi', 'external', 'discovery')})")
)

auc_sum <- data.table(
  test_data = c(glue("MGI ({ifelse(opt$discovery_cohort == 'mgi', 'discovery', 'external')})"), glue("UKB ({ifelse(opt$discovery_cohort == 'mgi', 'external', 'discovery')})")),
  cohort    = c("mgi", "ukb"),
  auc_est   = c(mgi_auc[2], ukb_auc[2]),
  auc_lo    = c(mgi_auc[1], ukb_auc[1]),
  auc_hi    = c(mgi_auc[3], ukb_auc[3]),
  auc_print = c(pretty_print(mgi_auc), pretty_print(ukb_auc))
)

auc_plot <- rbindlist(list(mgi_stuff, ukb_stuff))  |>
  ggplot(aes(x = 1 - specificity, y = sensitivity, color = test_data)) +
  geom_abline(lty = 3) +
  geom_path(linewidth = 2, alpha = 0.2) +
  geom_smooth(method = "loess", formula = "y ~ x",
              span = 0.5, se = FALSE, linewidth = 1) +
  scale_color_OkabeIto() +
  annotate(
    geom = "label", x = rep(0.75, 2), y = c(0.275, 0.225), label.size = NA,
    label = auc_sum[, auc_print], color = rev(palette_OkabeIto[1:nrow(auc_sum)])
  ) +
  labs(
    title    = glue("AUC for X{gsub('X', '', opt$outcome)} at t{opt$time_threshold}"),
    caption = str_wrap(glue("Discovery cohort: {toupper(opt$discovery_cohort)}; external cohort: {toupper(external_cohort)}; N_phecodes: {length(mgi_phers$phecodes$phecode)}; method: {opt$method}"), 60)
  ) +
  coord_equal() +
  cowplot::theme_minimal_grid() +
  theme(
    plot.caption    = element_text(hjust = 0),
    legend.position = "top",
    legend.title    = element_blank()
  )
ggsave(plot = auc_plot,
       filename = glue("{out_path}{comb_out_prefix}_auc.pdf"),
       width = 6, height = 6, device = cairo_pdf)

mgi_mod <- logistf(case ~ phers, data = mgi_phers$data)
ukb_mod <- logistf(case ~ phers, data = ukb_phers$data)

or_sum <- rbindlist(list(
  cbind(data.table(cohort = "mgi"), extractr_or(mgi_mod)),
  cbind(data.table(cohort = "ukb"), extractr_or(ukb_mod))
))

total_sum <- merge.data.table(
  auc_sum,
  or_sum,
  by = "cohort"
)

total_sum

mgi_phers_dist_plot <- top_or_plotr(phers_data = mgi_phers[["data"]],
                                     .title = glue("X{gsub('X', '', opt$outcome)}",
                                                   " PheRS distribution by case status at t{opt$time_threshold}"),
                                     .subtitle = glue("In MGI discovery cohort"),
                                     .caption = str_wrap(glue("Discovery cohort = {toupper(opt$discovery_cohort)}"),
                                                         width = 100))
ggsave(plot = mgi_phers_dist_plot,
       filename = glue("{out_path}{mgi_out_prefix}_phers_dist.pdf"),
       width = 8, height = 6, device = cairo_pdf)


ukb_phers_dist_plot <- top_or_plotr(phers_data = ukb_phers[["data"]],
                                    .title = glue("X{gsub('X', '', opt$outcome)}",
                                                  " PheRS distribution by case status at t{opt$time_threshold}"),
                                    .subtitle = glue("In UKB external cohort"),
                                    .caption = str_wrap(glue("Discovery cohort = {toupper(opt$discovery_cohort)}"),
                                                        width = 100))
ggsave(plot = ukb_phers_dist_plot,
       filename = glue("{out_path}{ukb_out_prefix}_phers_dist.pdf"),
       width = 8, height = 6, device = cairo_pdf)

## SAVE OUTPUT

saveRDS(
  mgi_phers,
  file = glue("{out_path}{mgi_out_prefix}_summary.rds")
)

saveRDS(
  ukb_phers,
  file = glue("{out_path}{ukb_out_prefix}_summary.rds")
)