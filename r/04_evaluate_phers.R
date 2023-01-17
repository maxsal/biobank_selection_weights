# evaluate phers
# author: max salvatore
# date:   20230117

# 1. libraries, functions, and options -----------------------------------------
library(data.table)
library(caret)
library(purrr)
library(progress)
library(pROC)
library(glue)
library(logistf)
library(optparse)

set.seed(61787)

for (i in list.files("fn/", full.names = TRUE)) source(i) # load functions
source(glue("https://raw.githubusercontent.com/umich-cphds/",
            "createUKBphenome/master/scripts/function.expandPhecodes.r"))

# optparse list ----------------------------------------------------------------
option_list <- list(
  make_option("--outcome", type = "character", default = "",
              help = "Outcome phecode"),
  make_option("--mgi_version", type = "character", default = "20210318",
              help = "Version of MGI data [default = 20210318]"),
  make_option("--ukb_version", type = "character", default = "20221117",
              help = "Version of UKB data [default = 20221117]"),
  make_option("--time_thresholds", type = "character", default = "0,1,2,3,5",
              help = glue("Time thresholds for the phenome data ",
                          "[default = 0,1,2,3,5]")),
  make_option("--mod_type", type = "character", default = "logistf",
              help = glue("Type of model to use in cooccurrence analysis - ",
                          "logistf for SPAtest [default = logistf]"))
)

parser <- OptionParser(usage="%prog [options]", option_list = option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

# 2. specifications ------------------------------------------------------------
outcome         <- opt$outcome
time_thresholds <- as.numeric(strsplit(opt$time_thresholds, ",")[[1]])

# model type - use "SPAtest" or "logistf" for cooccur analyses
mod_type <- "logistf"

# method
method <- "tophits"
# method <- "pwide_sig"
tophits_n <- 50

### files
file_paths <- get_files(mgi_version = opt$mgi_version,
                        ukb_version = opt$ukb_version)

# read data --------------------------------------------------------------------
## mgi
mgi_tr_pims <- list()
for (i in seq_along(time_thresholds)) {
  mgi_tr_pims[[i]] <- fread(
    glue("data/private/mgi/{opt$mgi_version}/X","{gsub('X', '', opt$outcome)}/",
         "time_restricted_phenomes/mgi_X{gsub('X', '', opt$outcome)}_t",
         "{time_thresholds[i]}_{opt$mgi_version}.txt")
  )
}
names(mgi_tr_pims) <- glue("t{time_thresholds}_threshold")

mgi_covariates <- fread(
  glue("data/private/mgi/{opt$mgi_version}/X{gsub('X', '', opt$outcome)}/",
       "matched_covariates.txt")
)

## ukb
ukb_pim0    <- fread(file_paths[["ukb"]]$pim0_file)
ukb_pim     <- fread(glue("data/private/ukb/{ukb_version}/",
                          "X{gsub('X', '', outcome)}/",
                          "time_restricted_phenomes/",
                          "ukb_X{gsub('X', '', outcome)}_",
                          "t{time_threshold}_{ukb_version}.txt"))
ukb_cooccur <- fread(glue("results/ukb/{ukb_version}/",
                          "X{gsub('X', '', outcome)}/",
                          "ukb_X{gsub('X', '', outcome)}_",
                          "t{time_threshold}_{ukb_version}_results.txt"))
ukb_covariates <- fread(glue("data/private/ukb/{ukb_version}/",
                             "X{gsub('X', '', outcome)}/",
                             "matched_covariates.txt"))

## phenome
pheinfo    <- fread("./data/phecode_mapping/data/Phecode_Definitions_FullTable_Modified.txt", colClasses = "character")

########################
### calculate phers ####
########################
phers_data <- purrr::map(time_thresholds,
                         ~calculate_phers(
                           pim       = tr_pims[[paste0("t", .x, "_threshold")]],
                           res       = cooccur_res[[paste0("t", .x, "_threshold")]],
                           method    = method,
                           tophits_n = tophits_n
                         ))
names(phers_data) <- paste0("t", time_thresholds)

for (i in time_thresholds) {
  phers_data[[paste0("t", i)]][["data"]] <- merge.data.table(
    phers_data[[paste0("t", i)]][["data"]],
    covariates,
    by = "id", all.x = TRUE
  )[, age_at_threshold := round(get(paste0("t", i, "_threshold"))/365.22, 1)]
}

################
### evaluate ###
################
eval_data <- purrr::map(time_thresholds,
                        ~evaluate_phers(
                          eval_pim   = phers_data[[paste0("t", .x)]][["data"]],
                          predictor  = "phers",
                          phers_phes = phers_data[[paste0("t", .x)]][["phecodes"]][, phecode],
                          covariates = c("age_at_threshold", "female"),
                          out_phe    = outcome,
                          pctile_or  = TRUE
                        ))
names(eval_data) <- paste0("t", time_thresholds)

for (i in time_thresholds) {
  phers_data[[paste0("t", i)]][["evaluation"]] <- eval_data[[paste0("t", i)]]
}

### save results!!
# save results as named nested list object
saveRDS(phers_data,
        file = paste0("./results/", version, "/X", gsub("X", "", outcome), "/mgi_X", gsub("X", "", outcome), "_", version, "_", ifelse(method == "tophits", paste0(method, tophits_n), method), "_eval.rds"))


# save evaluations as txt files
output_tsv <- function(data, names, ver = NULL, outc = NULL){ 
  if (is.null(ver)) {stop("Missing version (`ver`) argument. Check is `version` is specified in environment.")}
  if (is.null(outc)) {stop("Missing outcome (`outc`) arguement. Check if `outcome` is specified in environmnet.")}
  fwrite(data, paste0("./results/", ver, "/X", gsub("X", "", outc), "/mgi_X", gsub("X", "", outc), "_", names, "_", ver, "_", ifelse(method == "tophits", paste0(method, tophits_n), method), "_eval.txt"), sep = "\t")
}

for (i in time_thresholds) {
  output_tsv(data = phers_data[[paste0("t", i)]][["evaluation"]], names = paste0("t", i), ver = version, outc = outcome)
}



