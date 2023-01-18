# construct phers
# author: max salvatore
# date:   20230118

# 1. libraries, functions, and options (outcome agnostic) ----------------------
options(stringsAsFactors = FALSE)

library(data.table)
library(caret)
library(purrr)
library(progress)
library(pROC)
library(glue)
library(logistf)
library(cli)
library(optparse)

set.seed(61787)

for (i in list.files("fn/", full.names = TRUE)) source(i) # load functions

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
  make_option("--tophits_n", type = "numeric", default = "50",
              help = glue("Number of top hits to use in top hits PheRS ",
                          "[default = 50]"))
)

parser <- OptionParser(usage="%prog [options]", option_list = option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

# 2. specifications (specifies outcome) ----------------------------------------
mgi_version    <- opt$mgi_version  # mgi phenome version
ukb_version    <- opt$ukb_version  # ukb phenome version
outcome        <- opt$outcome      # outcome phecode
time_threshold <- opt$time_threshold

## extract file paths
file_paths <- get_files(mgi_version = opt$mgi_version,
                        ukb_version = opt$ukb_version)

## make sure output directories exist
if (!dir.exists(glue("results/mgi/{opt$mgi_version}/",
                "X{gsub('X', '', opt$outcome)}/phers/"))) {
  dir.create(glue("results/mgi/{opt$mgi_version}/",
                  "X{gsub('X', '', opt$outcome)}/phers/"),
             recursive = TRUE)
}
if (!dir.exists(glue("results/ukb/{opt$ukb_version}/",
                     "X{gsub('X', '', opt$outcome)}/phers/"))) {
  dir.create(glue("results/ukb/{opt$ukb_version}/",
                  "X{gsub('X', '', opt$outcome)}/phers/"),
             recursive = TRUE)
}

# 3. read data -----------------------------------------------------------------
cli_alert_info("loading data...")
## mgi
mgi_pim0    <- fread(file_paths[["mgi"]]$pim0_file)
mgi_pim     <- fread(glue("data/private/mgi/{mgi_version}/",
                          "X{gsub('X', '', outcome)}/time_restricted_phenomes/",
                          "mgi_X{gsub('X', '', outcome)}_",
                          "t{time_threshold}_{mgi_version}.txt"))
mgi_pim[is.na(mgi_pim)] <- 0
mgi_cooccur <- fread(glue("results/mgi/{mgi_version}/",
                          "X{gsub('X', '', outcome)}/",
                          "mgi_X{gsub('X', '', outcome)}_",
                          "t{time_threshold}_{mgi_version}_results.txt"))
mgi_covariates <- fread(glue("data/private/mgi/{mgi_version}/",
                             "X{gsub('X', '', outcome)}/",
                             "matched_covariates.txt"))

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

## phecode info
pheinfo <- fread("data/public/Phecode_Definitions_FullTable_Modified.txt",
                 colClasses = "character")

# 4. phecode exclusions ------------------------------------------------------
## exclusion range from PhewasCatalog
exclusionRange <- pheinfo[phecode == outcome, phecode_exclude_range]
exclusions1    <- pheinfo[phecode %in% unlist(unname(sapply(
  strsplit(exclusionRange,', {0,1}')[[1]],
  expandPhecodes))), phecode]

## phecodes not defined in both cohorts
exclusions2 <- pheinfo[!(phecode %in% gsub("X", "",
                                           intersect(names(mgi_pim0),
                                                     names(ukb_pim0)))),
                       phecode]

exclusionsX <- glue("X{union(exclusions1, exclusions2)}")

# 5. calculate naive phers -----------------------------------------------------
cli_alert_info("calculating naive phers for phecode {gsub('X', '', opt$outcome)}...")
## naive - pwide significant
### using MGI data
#### mgi
mgi_phers_dm_bn <- calculate_phers(
  pim        = mgi_pim,
  res        = mgi_cooccur[!(phecode %in% exclusionsX), ],
  method     = "pwide_sig",
  phers_name = glue("phers0_t{time_threshold}_dm_bn")
)
saveRDS(
  mgi_phers_dm_bn,
  file = glue("results/mgi/{opt$mgi_version}/",
              "X{gsub('X', '', opt$outcome)}/phers/",
              "mgi_phers_t{opt$time_threshold}_dm_bn.rds")
)

#### ukb
ukb_phers_dm_bn <- calculate_phers(
  pim        = ukb_pim,
  res        = mgi_cooccur[!(phecode %in% exclusionsX), ],
  method     = "pwide_sig",
  phers_name = glue("phers0_t{time_threshold}_dm_bn")
)
saveRDS(
  ukb_phers_dm_bn,
  file = glue("results/ukb/{opt$ukb_version}/",
              "X{gsub('X', '', opt$outcome)}/phers/",
              "ukb_phers_t{opt$time_threshold}_dm_bn.rds")
)

### using UKB data
#### mgi
mgi_phers_du_bn <- calculate_phers(
  pim        = mgi_pim,
  res        = ukb_cooccur[!(phecode %in% exclusionsX), ],
  method     = "pwide_sig",
  phers_name = glue("phers0_t{time_threshold}_du_bn")
)
saveRDS(
  mgi_phers_du_bn,
  file = glue("results/mgi/{opt$mgi_version}/",
              "X{gsub('X', '', opt$outcome)}/phers/",
              "mgi_phers_t{opt$time_threshold}_du_bn.rds")
)

#### ukb
ukb_phers_du_bn <- calculate_phers(
  pim        = ukb_pim,
  res        = ukb_cooccur[!(phecode %in% exclusionsX), ],
  method     = "pwide_sig",
  phers_name = glue("phers0_t{time_threshold}_du_bn")
)
saveRDS(
  ukb_phers_du_bn,
  file = glue("results/ukb/{opt$ukb_version}/",
              "X{gsub('X', '', opt$outcome)}/phers/",
              "ukb_phers_t{opt$time_threshold}_du_bn.rds")
)

## naive - # independent tests (i.e., # PCA explaining 99% of variance)
### using mgi data
mgi_pca <- fread(glue("results/mgi/{mgi_version}/",
                      "mgi_pca_importance_{mgi_version}.txt"))
mgi_pca[, stat := c("sd", "prop", "cum_prop")]
mgi_pca <- melt(mgi_pca, id.vars = "stat")[, pc := as.numeric(
  gsub("PC", "", variable)
)]
mgi_m_0.99 <- mgi_pca[stat == "cum_prop"][order(pc), ][value > 0.99, ][1, pc]

#### mgi
mgi_phers_dm_bm <- calculate_phers(
  pim        = mgi_pim,
  res        = mgi_cooccur[!(phecode %in% exclusionsX), ],
  method     = "pwide_sig",
  bonf_tests = mgi_m_0.99,
  phers_name = glue("phers0_t{time_threshold}_dm_bm")
)
saveRDS(
  mgi_phers_dm_bm,
  file = glue("results/mgi/{opt$mgi_version}/",
              "X{gsub('X', '', opt$outcome)}/phers/",
              "mgi_phers_t{opt$time_threshold}_dm_bm.rds")
)

#### ukb
ukb_phers_dm_bm <- calculate_phers(
  pim        = ukb_pim,
  res        = mgi_cooccur[!(phecode %in% exclusionsX), ],
  method     = "pwide_sig",
  bonf_tests = mgi_m_0.99,
  phers_name = glue("phers0_t{time_threshold}_dm_bm")
)
saveRDS(
  ukb_phers_dm_bm,
  file = glue("results/ukb/{opt$ukb_version}/",
              "X{gsub('X', '', opt$outcome)}/phers/",
              "ukb_phers_t{opt$time_threshold}_dm_bm.rds")
)

### using ukb data
ukb_pca <- fread(glue("results/ukb/{ukb_version}/",
                      "ukb_pca_importance_{ukb_version}.txt"))
ukb_pca[, stat := c("sd", "prop", "cum_prop")]
ukb_pca <- melt(ukb_pca, id.vars = "stat")[, pc := as.numeric(
  gsub("PC", "", variable)
)]
ukb_m_0.99 <- ukb_pca[stat == "cum_prop"][order(pc), ][value > 0.99, ][1, pc]

#### mgi
mgi_phers_du_bm <- calculate_phers(
  pim        = mgi_pim,
  res        = ukb_cooccur[!(phecode %in% exclusionsX), ],
  method     = "pwide_sig",
  bonf_tests = ukb_m_0.99,
  phers_name = glue("phers0_t{time_threshold}_du_bm")
)
saveRDS(
  mgi_phers_du_bm,
  file = glue("results/mgi/{opt$mgi_version}/",
              "X{gsub('X', '', opt$outcome)}/phers/",
              "mgi_phers_t{opt$time_threshold}_du_bm.rds")
)

#### ukb
ukb_phers_du_bm <- calculate_phers(
  pim        = ukb_pim,
  res        = ukb_cooccur[!(phecode %in% exclusionsX), ],
  method     = "pwide_sig",
  bonf_tests = ukb_m_0.99,
  phers_name = glue("phers0_t{time_threshold}_du_bm")
)
saveRDS(
  ukb_phers_du_bm,
  file = glue("results/ukb/{opt$ukb_version}/",
              "X{gsub('X', '', opt$outcome)}/phers/",
              "ukb_phers_t{opt$time_threshold}_du_bm.rds")
)

## naive - top 50 hits
### using mgi data
#### mgi
mgi_phers_dm_h50 <- calculate_phers(
  pim        = mgi_pim,
  res        = mgi_cooccur[!(phecode %in% exclusionsX), ],
  method     = "tophits",
  tophits_n  = opt$tophits_n,
  phers_name = glue("phers0_t{time_threshold}_dm_h{opt$tophits_n}")
)
saveRDS(
  mgi_phers_du_bn,
  file = glue("results/mgi/{opt$mgi_version}/",
              "X{gsub('X', '', opt$outcome)}/phers/",
              "mgi_phers_t{opt$time_threshold}_dm_h50.rds")
)

#### ukb
ukb_phers_dm_h50 <- calculate_phers(
  pim        = ukb_pim,
  res        = mgi_cooccur[!(phecode %in% exclusionsX), ],
  method     = "tophits",
  tophits_n  = opt$tophits_n,
  phers_name = glue("phers0_t{time_threshold}_dm_h{opt$tophits_n}")
)
saveRDS(
  ukb_phers_dm_h50,
  file = glue("results/ukb/{opt$ukb_version}/",
              "X{gsub('X', '', opt$outcome)}/phers/",
              "ukb_phers_t{opt$time_threshold}_dm_h50.rds")
)

### using ukb data
#### mgi
mgi_phers_du_h50 <- calculate_phers(
  pim        = mgi_pim,
  res        = ukb_cooccur[!(phecode %in% exclusionsX), ],
  method     = "tophits",
  tophits_n  = opt$tophits_n,
  phers_name = glue("phers0_t{time_threshold}_du_h{opt$tophits_n}")
)
saveRDS(
  mgi_phers_du_h50,
  file = glue("results/mgi/{opt$mgi_version}/",
              "X{gsub('X', '', opt$outcome)}/phers/",
              "mgi_phers_t{opt$time_threshold}_du_bh50.rds")
)
#### ukb
ukb_phers_du_h50 <- calculate_phers(
  pim        = ukb_pim,
  res        = ukb_cooccur[!(phecode %in% exclusionsX), ],
  method     = "tophits",
  tophits_n  = opt$tophits_n,
  phers_name = glue("phers0_t{time_threshold}_du_h{opt$tophits_n}")
)
saveRDS(
  ukb_phers_du_h50,
  file = glue("results/ukb/{opt$ukb_version}/",
              "X{gsub('X', '', opt$outcome)}/phers/",
              "ukb_phers_t{opt$time_threshold}_du_h{opt$tophits_n}.rds")
)

# 6. aggregate phers -----------------------------------------------------------
mgi_phers <- Reduce(merge.data.table,
                    list(
                      mgi_pim[, .(id, case)],
                      mgi_phers_dm_bn$data,
                      mgi_phers_du_bn$data,
                      mgi_phers_dm_h50$data,
                      mgi_phers_du_h50$data,
                      mgi_phers_dm_bm$data,
                      mgi_phers_du_bm$data
                    ))
fwrite(x = mgi_phers,
       file = glue("results/mgi/{opt$mgi_version}/",
                   "X{gsub('X', '', opt$outcome)}/phers/",
                   "mgi_naive_phers_t{opt$time_threshold}.txt"))
cli_alert_info("MGI PheRS naive AUCs")
quick_naive_aucs(x = mgi_phers)

ukb_phers <- Reduce(merge.data.table,
                    list(
                      ukb_pim[, .(id, case)],
                      ukb_phers_dm_bn$data,
                      ukb_phers_du_bn$data,
                      ukb_phers_dm_h50$data,
                      ukb_phers_du_h50$data,
                      ukb_phers_dm_bm$data,
                      ukb_phers_du_bm$data
                    ))
fwrite(x = ukb_phers,
       file = glue("results/ukb/{opt$ukb_version}/",
                   "X{gsub('X', '', opt$outcome)}/phers/",
                   "ukb_naive_phers_t{opt$time_threshold}.txt"))
cli_alert_info("UKB PheRS naive AUCs")
quick_naive_aucs(x = ukb_phers)
