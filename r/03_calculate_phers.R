# construct phers
# author: max salvatore
# date:   20221108

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

set.seed(61787)

for (i in list.files("fn/")) source(paste0("fn/", i)) # load functions
source(glue("https://raw.githubusercontent.com/umich-cphds/",
            "createUKBphenome/master/scripts/function.expandPhecodes.r"))

# 2. specifications (specifies outcome) ----------------------------------------
mgi_version           <- "20210318"       # mgi phenome version
ukb_version           <- "20221117"       # ukb phenome version
outcome               <- "184.1"          # outcome phecode (157 - PanCan,
                                          # 155 - LivCan, 184.1 - OvCan)
# time_thresholds       <- c(0, 1, 2, 3, 5) # time thresholds
time_threshold <- 1
method                <- "tophits"        # 'tophits' or 'pwide_sig'

## pull file paths corresponding to the data version specified
file_paths <- get_files(mgi_version = mgi_version, ukb_version = ukb_version)

# 3. read data -----------------------------------------------------------------
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
## naive - pwide significant
### using MGI data
#### mgi
mgi_phers <- calculate_phers(
  pim        = mgi_pim,
  res        = mgi_cooccur[!(phecode %in% exclusionsX), ],
  method     = "pwide_sig",
  phers_name = glue("pwide_fm_t{time_threshold}_phers")
  )$data

#### ukb
ukb_phers <- calculate_phers(
    pim        = ukb_pim,
    res        = mgi_cooccur[!(phecode %in% exclusionsX), ],
    method     = "pwide_sig",
    phers_name = glue("pwide_fm_t{time_threshold}_phers")
    )$data

### using UKB data
#### mgi
mgi_phers <- merge.data.table(
  mgi_phers,
  calculate_phers(
  pim        = mgi_pim,
  res        = ukb_cooccur[!(phecode %in% exclusionsX), ],
  method     = "pwide_sig",
  phers_name = glue("pwide_fu_t{time_threshold}_phers")
)$data )

#### ukb
ukb_phers <- merge.data.table(
  ukb_phers,
  calculate_phers(
  pim        = ukb_pim,
  res        = ukb_cooccur[!(phecode %in% exclusionsX), ],
  method     = "pwide_sig",
  phers_name = glue("pwide_fu_t{time_threshold}_phers")
)$data )

## naive - # independent tests (i.e., # PCA explaining 99% of variance)
############# WHAT'S GOING ON BENEATH HERE #####################################

mgi_keepers <- names(mgi_pim)[!names(mgi_pim) %in%
                                c("id", "IID", "case",
                                  glue("X{gsub('X', '', outcome)}"),
                                  exclusionsX)]
mgi_tpca <- prcomp(mgi_pim[, ..mgi_keepers], center = FALSE, scale. = FALSE)
mgi_tpcs <- mgi_tpca$x[, 1:which.min(
  abs(summary(mgi_tpca)$importance[3, ] - 0.95)
  )]
max_mgi_pcs <- apply(abs(mgi_tpcs), 2, max)
mgi_tpcs_mod <- sweep(mgi_tpcs, MARGIN = 2, max_mgi_pcs, "/")
mgi_tpcs_mod <- mgi_tpcs_mod + 1

mgi_rotations <- data.table(mgi_tpca$rotation)
mgi_rotations <- mgi_rotations[,1:which.min(
  abs(summary(mgi_tpca)$importance[3, ] - 0.95)
  )]
#########
mgi_pca <- fread(glue("results/mgi/{mgi_version}/",
                      "mgi_pca_importance_{mgi_version}.txt"))
mgi_pca[, stat := c("sd", "prop", "cum_prop")]
mgi_pca <- melt(mgi_pca, id.vars = "stat")[, pc := as.numeric(
  gsub("PC", "", variable)
  )]
mgi_pca[stat == "cum_prop"][order(pc), ][value > 0.99, ][1, pc]

## naive - top 50 hits
### mgi
mgi_phers <- merge(
  mgi_phers,
  calculate_phers(
  pim        = mgi_pim,
  res        = mgi_cooccur[!(phecode %in% exclusionsX), ],
  method     = "tophits",
  tophits_n  = 50,
  phers_name = glue("th50_fm_t{time_threshold}_phers")
)$data) 
### ukb
ukb_phers <- merge.data.table(
  ukb_phers,
  calculate_phers(
    pim        = ukb_pim,
    res        = ukb_cooccur[!(phecode %in% exclusionsX), ],
    method     = "tophits",
    tophits_n  = 50,
    phers_name = glue("th50_fm_t{time_threshold}_phers")
  )$data )
