# perform random forest analysis for given threshold and discovery cohort
# requires: time-threshold phecode indicator matrices must already exist
# outputs:  auc plot, vip plot, summary list object, optimal random forest object
# author:   max salvatore
# date:     20230223

# libraries --------------------------------------------------------------------
library(fst)
library(vip)
library(data.table)
library(snakecase)
library(stringr)
library(cowplot)
library(cli)
library(optparse)
library(glue)
library(colorblindr)
library(rsample)      # data splitting
library(ranger)       # a faster implementation of randomForest
library(caret)        # an aggregator package for performing many 

# optparse list ---
option_list <- list(
  make_option("--outcome", type = "character", default = "157",
              help = "Outcome phecode [default = %default]"),
  make_option("--mgi_version", type = "character", default = "20220822",
              help = "Version of MGI data [default = %default]"),
  make_option("--ukb_version", type = "character", default = "20221117",
              help = "Version of UKB data [default = %default]"),
  make_option("--time_threshold", type = "numeric", default = "0",
              help = glue("Time threshold for the phenome data ",
                          "[default = 0]")),
  make_option("--split_prop", type = "numeric", default = "0.7",
              help = glue("Proportion of data in training set ",
                          "[default = %default]")),
  make_option("--strata", type = "character", default = "case",
              help = glue("Strata (outcome variable) for distributing between train/test sets ",
                          "[default = %default]")),
  make_option("--n_trees", type = "numeric", default = "500",
              help = glue("Number of trees to grow in random forest ",
                          "[default = %default]")),
  make_option("--n_vip", type = "numeric", default = "20",
              help = glue("Number of variables to include in variable importance plot ",
                          "[default = %default]")),
  make_option("--seed", type = "numeric", default = "1234",
              help = glue("Set seed for reproducibility ",
                          "[default = %default]")),
  make_option("--n_core_prop", type = "numeric", default = "0.5",
              help = glue("Proportion of available cores to use ",
                          "[default = %default]")),
  make_option("--discovery_cohort", type = "character", default = "mgi",
              help = glue("Cohort to use as discovery cohort (mgi / ukb) ",
                          "[default = %default]"))
)
parser <- OptionParser(usage="%prog [options]", option_list = option_list)
args   <- parse_args(parser, positional_arguments = 0)
opt    <- args$options
print(opt)

if (parallel::detectCores() == 1) {
  cli_alert_warning("Only 1 core detected - this could be a while!")
  n_cores <- 1
} else {
  n_cores <- parallel::detectCores() * opt$n_core_prop
}
source("fn/expandPhecodes.R")
source("fn/files-utils.R")

## extract file paths - UNNECESSARY?!
file_paths <- get_files(mgi_version = opt$mgi_version,
                        ukb_version = opt$ukb_version)

# check output folder exists ---------------------------------------------------
out_path <- glue("results/{coh}/{coh_version}/X{outc}/random_forest/",
     coh = opt$discovery_cohort,
     coh_version = ifelse(opt$discovery_cohort == "mgi", opt$mgi_version, ifelse(opt$discovery_cohort == "ukb", opt$ukb_version, NA)),
     outc = gsub("X", "", opt$outcome))
if (!dir.exists(out_path)) {
  dir.create(out_path, recursive = TRUE)
}

# data -------------------------------------------------------------------------
cli_alert("reading data...")
d <- read_fst(glue("data/private/mgi/{opt$mgi_version}/X{gsub('X', '', opt$outcome)}/time_restricted_phenomes/mgi_X{gsub('X', '', opt$outcome)}_t{opt$time_threshold}_{opt$mgi_version}.fst"),
              as.data.table = TRUE)
d_c <- read_fst(glue("data/private/mgi/{opt$mgi_version}/X{gsub('X', '', opt$outcome)}/matched_covariates.fst"),
                as.data.table = TRUE)
d <- merge.data.table(
  d,
  d_c[, .(id, age_at_threshold = round(get(glue("t{opt$time_threshold}_threshold")) / 365.25, 3), female)],
  by = "id",
  all.x = TRUE
)
d <- d[, !c("id")]

u <- read_fst(glue("data/private/ukb/{opt$ukb_version}/X{gsub('X', '', opt$outcome)}/time_restricted_phenomes/ukb_X{gsub('X', '', opt$outcome)}_t{opt$time_threshold}_{opt$ukb_version}.fst"),
              as.data.table = TRUE)
u_c <- read_fst(glue("data/private/ukb/{opt$ukb_version}/X{gsub('X', '', opt$outcome)}/matched_covariates.fst"),
                as.data.table = TRUE)
u <- merge.data.table(
  u,
  u_c[, .(id, age_at_threshold = round(get(glue("t{opt$time_threshold}_threshold")) / 365.25, 3), female)],
  by = "id",
  all.x = TRUE
)
u <- u[, !c("id")]

p <- fread("data/public/Phecode_Definitions_FullTable_Modified.txt",
           colClasses = "character")

## exclusion range from PhewasCatalog
exclusionRange <- p[phecode == gsub("X", "", opt$outcome), phecode_exclude_range]
exclusions1    <- p[phecode %in% unlist(unname(sapply(
  strsplit(exclusionRange,', {0,1}')[[1]],
  expandPhecodes))), phecode]
exclusions2 <- p[leaf == 0, phecode]

p[, `:=` (
  phecode = paste0("X", phecode),
  group = snakecase::to_sentence_case(group)
)]

keep_these_vars <- intersect(names(d), names(u))
keep_these_vars <- keep_these_vars[!(keep_these_vars %in% paste0("X", union(exclusions1, exclusions2)))]

cli_alert_info("using {format(length(keep_these_vars) - 1, big.mark = ',')} predictors (only 'leaf' phecodes) available in both cohorts...")
d <- d[, ..keep_these_vars]
u <- u[, ..keep_these_vars]

if (opt$discovery_cohort == "mgi") {
  data     <- d
  external <- u
} else if (opt$discover_cohort == "ukb") {
  data     <- u
  external <- d
} else {
  stop("opt$discovery_cohort must be one of 'mgi' or 'ukb'")
}

# below is function or script
cli_alert("splitting data...")
if (!is.null(opt$strata)) {
  data_split <- initial_split(data,
                              prop = opt$split_prop,
                              strata = opt$strata)
} else {
  data_split <- initial_split(data,
                              prop = opt$split_prop)
}

data_train <- training(data_split)
data_test  <- testing(data_split)

# full grid search
hyper_grid <- expand.grid(
  mtry        = round(seq(5, round(ncol(data_train) / 2), length.out = 5)),
  node_size   = round(seq(3, round(nrow(data_train) / 100), length.out = 5)),
  sample_size = seq(0.5, 0.8, length.out = 4),
  oob_rmse    = 0
) |> as.data.table()


cli_progress_bar(total = nrow(hyper_grid), "Optimizing hyperparameters")
for (i in 1:nrow(hyper_grid)) {
  model <- ranger(
    formula         = case ~ .,
    data            = data_train,
    num.trees       = opt$n_trees,
    mtry            = hyper_grid$mtry[i],
    min.node.size   = hyper_grid$node_size[i],
    sample.fraction = hyper_grid$sample_size[i],
    num.threads     = n_cores
  )
  hyper_grid$oob_rmse[i] <- sqrt(model$prediction.error)
  cli_progress_update()
}

hyper_grid[order(oob_rmse),][1, ]

oob_rmse <- vector(mode = "numeric", length = 100)

cli_progress_bar(total = length(oob_rmse), "run using optimal settings...")
for (i in seq_along(oob_rmse)) {
  optimal_ranger <- ranger(
    formula         = case ~ .,
    data            = data_train,
    num.trees       = opt$n_trees,
    mtry            = hyper_grid[order(oob_rmse), ][1, mtry],
    min.node.size   = hyper_grid[order(oob_rmse), ][1, node_size],
    sample.fraction = hyper_grid[order(oob_rmse), ][1, sample_size],
    importance      = "impurity",
    num.threads     = n_cores
  )
  oob_rmse[i] <- sqrt(optimal_ranger$prediction.error)
  cli_progress_update()
}

# predicting -------------------------------------------------------------------
suppressMessages({
  pred_opt     <- predict(optimal_ranger, data_test, type = "response")
  pred_opt_roc <- pROC::roc(data_test[, case], pred_opt[["predictions"]])
  pred_opt_auc <- pROC::ci.auc(data_test[, case], pred_opt[["predictions"]])
  
  pred_oth     <- predict(optimal_ranger, external, type = "response")
  pred_oth_roc <- pROC::roc(external[, case], pred_oth[["predictions"]])
  pred_oth_auc <- pROC::ci.auc(external[, case], pred_oth[["predictions"]])
})

pretty_print <- function(x, r = 3) {
  paste0( round(x[2], r), " (", round(x[1], r), ", ", round(x[3], r), ")")
}

in_stuff <- data.table(
  sensitivity = pred_opt_roc$sensitivities,
  specificity = pred_opt_roc$specificities,
  test_data   = "Hold-out test data"
)
out_stuff <- data.table(
  sensitivity = pred_oth_roc$sensitivities,
  specificity = pred_oth_roc$specificities,
  test_data   = "External data"
)

auc_sum <- data.table(
  test_data = c("Hold-out test data", "External data"),
  cohort    = c(opt$discovery_cohort, ifelse(opt$discovery_cohort == "mgi", "ukb", "mgi")),
  auc_est   = c(pred_opt_auc[2], pred_oth_auc[2]),
  auc_lo    = c(pred_opt_auc[1], pred_oth_auc[1]),
  auc_hi    = c(pred_opt_auc[3], pred_oth_auc[3]),
  auc_print = c(pretty_print(pred_opt_auc), pretty_print(pred_oth_auc))
)

auc_plot <- rbindlist(list(in_stuff, out_stuff))  |>
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
    subtitle = glue("discovery = {opt$discovery_cohort}, external = {ifelse(opt$discovery_cohort == 'mgi', 'ukb', 'mgi')}"),
    caption  = str_wrap(glue("N_trees = {opt$n_trees}, mtry = {hyper_grid[order(oob_rmse), ][1, mtry]}, node size = {hyper_grid[order(oob_rmse), ][1, node_size]}, sample fraction = {hyper_grid[order(oob_rmse), ][1, sample_size]}"), width = 100)
  ) +
  coord_equal() +
  cowplot::theme_minimal_grid() +
  theme(
    legend.position = "top",
    legend.title    = element_blank()
  )
ggsave(plot = auc_plot,
       filename = glue("{out_path}mgid_ukbe_X{gsub('X', '', opt$outcome)}_t{opt$time_threshold}_rf_auc.pdf"),
       width = 6, height = 6, device = cairo_pdf)

## VIP plot
vip_data <- data.table(
  variable = names(optimal_ranger$variable.importance),
  importance = optimal_ranger$variable.importance
)[order(-importance), ][1:opt$n_vip, ] |>
  merge.data.table(p[, .(phecode, description, group, color)], by.x = "variable", by.y = "phecode", all.x = TRUE)
vip_data[, full_name := ifelse(is.na(description), variable, paste0(variable, ": ", description))]

cols <- unique(vip_data[, .(group, color)])
colors <- cols[, color]
names(colors) <- cols[, group]

vip_plot <- vip_data |>
  ggplot(aes(x = reorder(full_name, importance), y = importance, fill = group)) +
  geom_bar(stat = "identity") +
  labs(
    title    = glue("VIP for X{gsub('X', '', opt$outcome)} at t{opt$time_threshold}"),
    subtitle = glue("discovery = {opt$discovery_cohort}, external = {ifelse(opt$discovery_cohort == 'mgi', 'ukb', 'mgi')}"),
    caption  = str_wrap(glue("N_trees = {opt$n_trees}, mtry = {hyper_grid[order(oob_rmse), ][1, mtry]}, node size = {hyper_grid[order(oob_rmse), ][1, node_size]}, sample fraction = {hyper_grid[order(oob_rmse), ][1, sample_size]}"), width = 100),
    x = "",
    y = "Importance"
  ) +
  scale_fill_manual(values = colors) +
  scale_x_discrete(labels = \(x) stringr::str_wrap(x, width = 50)) +
  guides(fill = guide_legend(ncol = 3, byrow = TRUE)) +
  coord_flip() +
  cowplot::theme_minimal_grid() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  )
ggsave(plot = vip_plot,
       filename = glue("{out_path}mgid_ukbe_X{gsub('X', '', opt$outcome)}_t{opt$time_threshold}_vip.pdf"),
       width = 10, height = 8, device = cairo_pdf)

list(
  "discovery_cohort"      = opt$discovery_cohort,
  "external_cohort"       = ifelse(opt$discovery_cohort == "mgi", "ukb", "mgi"),
  "auc_summary"           = auc_sum,
  "discovery_sense_spec"  = in_stuff,
  "external_sense_spec"   = out_stuff,
  "vip_table"             = data.table(
    variable = names(optimal_ranger$variable.importance),
    importance = optimal_ranger$variable.importance
  )[order(-importance), ],
  "optimal_forest_object" = optimal_ranger,
  "hyper_grid"            = hyper_grid,
  "auc_plot"              = auc_plot,
  "vip_plot"              = vip_plot,
  "out_path"              = out_path,
  "optparse_list"         = opt
) |> saveRDS(file = glue("{out_path}mgid_ukbe_X{gsub('X', '', opt$outcome)}_t{opt$time_threshold}_summary.rds"))
  
saveRDS(optimal_ranger, file = glue("{out_path}mgid_ukbe_X{gsub('X', '', opt$outcome)}_t{opt$time_threshold}_optimal_rf.rds"))

cli_alert_success("script success! see output in {.path {out_path}}")