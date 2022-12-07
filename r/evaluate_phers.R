#################
### libraries ###
#################
library(data.table)
library(caret)
library(purrr)
library(progress)
library(pROC)
library(glue)
library(logistf)

set.seed(61787)

#############
### specs ###
#############
# data version
version <- "20210318"

# ukb phenome version
ukb_phenome_version <- "20221020"

# outcome phecode
# outcome <- "157" # pancreatic cancer
# outcome <- "155" # liver cancer
outcome <- "184.1" # ovarian cancer

# time thresholds
time_thresholds <- c(0, 1, 2, 3, 5)

# model type - use "SPAtest" or "logistf" for cooccur analyses
mod_type <- "logistf"

# method
method <- "tophits"
# method <- "pwide_sig"
tophits_n <- 50

######################
### load functions ###
######################
purrr::walk(list.files("fn/"), ~source(paste0("fn/", .x)))

### files
file_paths <- get_files(data_version = version)

#################
### read data ###
#################
tr_pims <- purrr::map(time_thresholds,
                      ~fread(paste0("./data/", version, "/processed/X", gsub("X", "", outcome), "/time_restricted_phenomes/mgi_X", gsub("X", "", outcome), "_t", .x, "_", version, ".txt")))
names(tr_pims) <- paste0("t", time_thresholds, "_threshold")

cooccur_res <- purrr::map(time_thresholds,
                          ~fread(paste0("./results/", version, "/X", gsub("X", "", outcome), "/mgi_X", gsub("X", "", outcome), "_t", .x, "_", version, "_results.txt")))
names(cooccur_res) <- paste0("t", time_thresholds, "_threshold")

covariates <- fread(paste0("./data/", version, "/processed/X", gsub("X", "", outcome), "/matched_covariates.txt"))
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



