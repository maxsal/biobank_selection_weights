# construct phers
# author: max salvatore
# date:   20221108

# 1. libraries, functions, and options (outcome agnostic) ----------------------
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(caret)
  library(purrr)
  library(pROC)
  library(glue)
  library(logistf)
})

set.seed(61787)

for (i in list.files("fn/", full.names = TRUE)) source(i)

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
mgi_phers_dm_bn <- calculate_phers(
  pim        = mgi_pim,
  res        = mgi_cooccur[!(phecode %in% exclusionsX), ],
  method     = "pwide_sig",
  phers_name = glue("phers0_t{time_threshold}_dm_bn")
  )

#### ukb
ukb_phers_dm_bn <- calculate_phers(
    pim        = ukb_pim,
    res        = mgi_cooccur[!(phecode %in% exclusionsX), ],
    method     = "pwide_sig",
    phers_name = glue("phers0_t{time_threshold}_dm_bn")
    )$data

### using UKB data
#### mgi
mgi_phers_du_bn <- calculate_phers(
  pim        = mgi_pim,
  res        = ukb_cooccur[!(phecode %in% exclusionsX), ],
  method     = "pwide_sig",
  phers_name = glue("phers0_t{time_threshold}_du_bn")
)

#### ukb
ukb_phers_du_bn <- calculate_phers(
  pim        = ukb_pim,
  res        = ukb_cooccur[!(phecode %in% exclusionsX), ],
  method     = "pwide_sig",
  phers_name = glue("phers0_t{time_threshold}_du_bn")
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
#### ukb
ukb_phers_dm_bm <- calculate_phers(
  pim        = ukb_pim,
  res        = mgi_cooccur[!(phecode %in% exclusionsX), ],
  method     = "pwide_sig",
  bonf_tests = mgi_m_0.99,
  phers_name = glue("phers0_t{time_threshold}_dm_bm")
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
#### ukb
ukb_phers_du_bm <- calculate_phers(
  pim        = ukb_pim,
  res        = ukb_cooccur[!(phecode %in% exclusionsX), ],
  method     = "pwide_sig",
  bonf_tests = ukb_m_0.99,
  phers_name = glue("phers0_t{time_threshold}_du_bm")
)
## naive - top 50 hits
### using mgi data
#### mgi
mgi_phers_dm_h50 <- calculate_phers(
  pim        = mgi_pim,
  res        = mgi_cooccur[!(phecode %in% exclusionsX), ],
  method     = "tophits",
  tophits_n  = 50,
  phers_name = glue("phers0_t{time_threshold}_dm_h50")
)
#### ukb
ukb_phers_dm_h50 <- calculate_phers(
    pim        = ukb_pim,
    res        = mgi_cooccur[!(phecode %in% exclusionsX), ],
    method     = "tophits",
    tophits_n  = 50,
    phers_name = glue("phers0_t{time_threshold}_dm_h50")
  )

### using ukb data
#### mgi
mgi_phers_du_h50 <- calculate_phers(
    pim        = mgi_pim,
    res        = ukb_cooccur[!(phecode %in% exclusionsX), ],
    method     = "tophits",
    tophits_n  = 50,
    phers_name = glue("phers0_t{time_threshold}_du_h50")
  )
#### ukb
ukb_phers_du_h50 <- calculate_phers(
    pim        = ukb_pim,
    res        = ukb_cooccur[!(phecode %in% exclusionsX), ],
    method     = "tophits",
    tophits_n  = 50,
    phers_name = glue("phers0_t{time_threshold}_du_h50")
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

test_fn <- function(x) {
  
  x_arg <- deparse(substitute(x))
  phers_name <- x[["phers_name"]]
  
  disc <- ifelse(
    str_extract(phers_name, "(?<=_d)[:alpha:]") == "m",
    "mgi",
    "ukb")
  appl       <- substr(x_arg, 1, 3)
  phers_type <- ifelse(
    substr(gsub("_d[a-z]{1}_", "", str_extract(x_arg, "_d[a-z]{1}_[a-z]{1}.*")),
    1, 1) == "h", "top hits", "phenome wide significance")
  if (phers_type == "top hits") {
    bonf_tests <- NA
    top_hits <- as.numeric(gsub("_h", "",
                                str_extract(phers_name, "_h[0-9]{1,2}")))
  } else {
    bonf_tests <- ifelse(
     str_extract(phers_name, "_b[a-z]{1}") == "n",
     "no. of models", "estimated ind. tests (based on PCA @ 99%)"
    )
    top_hits <- NA
  }
  
  data.table(
    name = deparse(substitute(x)),
    discovery_cohort = dis,
    application_cohort = apple,
    phers_type = phers_type,
    n_phecodes = x[["n_phecodes"]],
    bonferroni_tests = bonf_tests,
    top_hits = top_hits,
    weights = x[["phecodes"]]
  )
}
test_fn(x = mgi_phers_dm_bn)
#########
out <- data.table()
for (i in names(mgi_phers)[grepl("phers", names(mgi_phers))]) {
  out <- rbindlist(list(
    out,
    data.table(
    phers = i,
    auc = pROC::roc(response = mgi_phers[["case"]],
              predictor = mgi_phers[[i]],
              family = binomial(),
              ci = TRUE)$auc
  )), fill = TRUE)
}
out[order(-auc)]

