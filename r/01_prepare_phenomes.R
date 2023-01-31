# prepare mgi and ukb data for time-restricted phecode-phecode phewas
# outputs: time-threshold phecode indicator matrices
# author:  max salvatore
# date:    20221207

# 1. libraries, functions, and options (outcome agnostic) ----------------------
options(stringsAsFactors = FALSE)

library(data.table)
library(MatchIt)
library(cli)
library(optparse)
library(glue)

set.seed(61787)

for (i in list.files("fn/")) source(paste0("fn/", i)) # load functions

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
  make_option("--nearest_matching_var", type = "character",
              default = "age_at_first_diagnosis,length_followup",
              help = glue("Matching variables by nearest  [default = ",
                          "age_at_first_diagnosis,length_followup]")),
  make_option("--exact_matching_var", type = "character", default = "female",
              help = "Matching variables by exact [default = female]"),
  make_option("--matching_caliper", type = "numeric", default = "0.25",
              help = "Matching caliper [default = 0.25]"),
  make_option("--matching_ratio", type = "numeric", default = "2",
              help = "Number of non-cases to match per case [default = 2]")
)

parser <- OptionParser(usage = "%prog [options]", option_list = option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

time_thresholds <- as.numeric(strsplit(opt$time_thresholds, ",")[[1]])

# 2. specifications (specifies outcome) --------------------------------------
time_thresholds <- as.numeric(strsplit(opt$time_thresholds, ",")[[1]])
nearest_matching_vars <- strsplit(opt$nearest_matching_var, ",")[[1]]
exact_matching_vars <- strsplit(opt$exact_matching_var, ",")[[1]]

# 3. extra preparations (outcome-specific) -------------------------------------
## confirm file structure for a given outcome exists - if not, create paths
### mgi
check_folder_structure(
  cohort          = "mgi",
  data_version    = opt$mgi_version,
  outcome_phecode = opt$outcome
  )

### ukb
check_folder_structure(
  cohort          = "ukb",
  data_version    = opt$ukb_version,
  outcome_phecode = opt$outcome
)

## pull file paths corresponding to the data version specified
file_paths <- get_files(mgi_version = opt$mgi_version,
                        ukb_version = opt$ukb_version)

# 4. read data -----------------------------------------------------------------
## mgi
cli_alert_info("loading mgi data...")
### phecode indicator matrix (PEDMASTER_0)
mgi_pim0 <- fread(file_paths[["mgi"]]$pim0_file)
setnames(mgi_pim0,
                     old = c("IID", paste0("X", opt$outcome)),
                     new = c("id", "outcome"))

### demographics
mgi_demo <- fread(file_paths[["mgi"]]$demo_file)[, .(
  id       = Deid_ID,
  age      = Age,
  alive    = as.numeric(AliveYN == "Y"),
  dead_dsb = Deceased_DaysSinceBirth,
  ethn     = EthnicityName,
  marital  = MaritalStatusCode,
  sex      = Sex,
  race     = RaceName)]

### icd-phecode data
mgi_full_phe <- get_full_icd_dsb_phecode(
  icd_file      = file_paths[["mgi"]]$icd9_file,
  more_icd_file = file_paths[["mgi"]]$icd10_file
)
mgi_first_phe <- mgi_full_phe[
  mgi_full_phe[, .I[which.min(dsb)], by = c("id", "phecode")]$V1
  ]

## ukb
cli_alert_info("loading ukb data...")
### phecode indicator matrix (PEDMASTER_0)
ukb_pim0 <- fread(file_paths[["ukb"]]$pim0_file)
setnames(ukb_pim0,
                     old = c("IID", paste0("X", gsub("X", "", opt$outcome))),
                     new = c("id", "outcome"))
### demographics
ukb_demo <- fread(file_paths[["ukb"]]$demo_file,
                              na.strings = c("", "NA", "."),
                              colClass = "character")
ukb_demo <- ukb_demo[, .(
  id   = as.character(id),
  dob  = as.Date(dob),
  age  = floor(as.numeric(as.Date("2022-11-17") - as.Date(dob)) / 365.25),
  ethn = ethnicity,
  sex)]
ukb_demo <- ukb_demo[complete.cases(ukb_demo), ]

### icd-phecode data
ukb_full_phe <- get_full_icd_dsb_phecode(
  icd_file      = file_paths[["ukb"]]$icd_phecode_file
)
ukb_first_phe <- ukb_full_phe[
  ukb_full_phe[, .I[which.min(dsb)], by = c("id", "phecode")]$V1
]

# 5. identify cases ------------------------------------------------------------
## mgi
mgi_case <- unique(mgi_first_phe[phecode == opt$outcome])
mgi_case_ids <- mgi_case[, unique(id)]

## ukb
ukb_case <- unique(ukb_first_phe[phecode == opt$outcome])
ukb_case_ids <- ukb_case[, unique(id)]

# 6. calculate diagnostic metrics ----------------------------------------------
## mgi
cli_alert_info("calculating diagnostic metrics in mgi...")
mgi_diag_metrics <- get_icd_phecode_metrics(
  full_phe_data = mgi_full_phe
)
mgi_matching_cov <- merge.data.table(
  mgi_demo[, .(id, age, female = as.numeric(sex == "F"))],
  mgi_diag_metrics,
  by = "id"
)[, case := fifelse(id %in% mgi_case_ids, 1, 0)]

## ukb
cli_alert_info("calculating diagnostic metrics in ukb...")
ukb_diag_metrics <- get_icd_phecode_metrics(
  full_phe_data = ukb_first_phe
)[, id := as.character(id)]
ukb_matching_cov <- merge.data.table(
  ukb_demo[, .(id, age, female = as.numeric(sex == "Female"))],
  ukb_diag_metrics,
  by = "id"
)[, case := fifelse(id %in% ukb_case_ids, 1, 0)]

# 7. perform matching ----------------------------------------------------------
## mgi
cli_alert_info(glue("performing 1:{opt$matching_ratio} case:non-case ",
                         "matching in mgi..."))
mgi_match_text <- glue("MatchIt::matchit(case ~ ",
                       "{glue_collapse(c(c(nearest_matching_vars, ",
                       "exact_matching_vars), sep = ' + ')}, ",
                       "data = mgi_matching_cov, calclosest = TRUE, ",
                       "mahvars = c({paste0(sapply(nearest_matching_vars, ",
                       "\(x) paste0('\\'', x, '\\'')), collapse = ', ')}), ",
                       "caliper = {opt$matching_caliper}, ",
                       "exact = c({paste0(sapply(exact_matching_vars, ",
                       "\(x) paste0('\\'', x, '\\'')), collapse = ', ')}), ",
                       "ratio = {opt$matching_ratio})")
mgi_match <- eval(parse(text = mgi_match_text))
mgi_matched <- MatchIt::match.data(mgi_match)
mgi_post_match_cov <- merge.data.table(
  mgi_matched,
  merge.data.table(
    mgi_case,
    mgi_matched,
    by = "id",
    all.x = TRUE
  )[, .(subclass, case_dsb = dsb)],
  by = "subclass"
)
for (i in time_thresholds) {
  mgi_post_match_cov[, (glue("t{i}_threshold")) :=
                       floor(case_dsb - (365.25 * i))]
}
for (i in time_thresholds) {
  mgi_post_match_cov[, (glue("t{i}_indicator")) :=
                       as.numeric(all(get(glue("t{i}_threshold")) > first_dsb)),
                     by = subclass]
}
if ( !dir.exists( glue("data/private/mgi/{opt$mgi_version}/",
                       "X{gsub('X','', opt$outcome)}/",
                       "time_restricted_phenomes/") ) ) {
  dir.create( glue("data/private/mgi/{opt$mgi_version}/", 
                   "X{gsub('X','', opt$outcome)}/time_restricted_phenomes/"),
              recursive = TRUE )
}
### save mgi matching data
fwrite(mgi_post_match_cov,
       glue("data/private/mgi/{opt$mgi_version}/",
            "X{gsub('X', '', opt$outcome)}/matched_covariates.txt"),
       sep = "\t")

## ukb
cli_alert_info(glue("performing 1:{opt$matching_ratio} case:non-case ",
                    "matching in ukb..."))
ukb_match_text <- glue("matchit(case ~ ",
                       "{glue_collapse(c(c(nearest_matching_vars, ",
                       "exact_matching_vars), sep = ' + ')}, ",
                       "data = ukb_matching_cov, calclosest = TRUE, ",
                       "mahvars = c({paste0(sapply(nearest_matching_vars,",
                       "\(x) paste0('\\'', x, '\\'')), collapse = ', ')}), ",
                       "caliper = {opt$matching_caliper}, ",
                       "exact = c({paste0(sapply(exact_matching_vars,",
                       "\(x) paste0('\\'', x, '\\'')), collapse = ', ')}), ",
                       "ratio = {opt$matching_ratio})")
ukb_match <- eval(parse(text = ukb_match_text))
ukb_matched <- match.data(ukb_match)
ukb_post_match_cov <- merge.data.table(
  ukb_matched,
  merge.data.table(
    ukb_case,
    ukb_matched,
    by = "id",
    all.x = TRUE
  )[, .(subclass, case_dsb = dsb)],
  by = "subclass"
)
for (i in time_thresholds) {
  ukb_post_match_cov[, (glue("t{i}_threshold")) :=
                       floor(case_dsb - (365.25 * i))]
}
for (i in time_thresholds) {
  ukb_post_match_cov[, (glue("t{i}_indicator")) :=
                       as.numeric(all(get(glue("t{i}_threshold")) > first_dsb)),
                     by = subclass]
}
if ( !dir.exists( glue("data/private/ukb/{opt$ukb_version}/",
                       "X{gsub('X','', opt$outcome)}/",
                       "time_restricted_phenomes/") ) ) {
  dir.create( glue("data/private/ukb/{opt$ukb_version}/",
                   "X{gsub('X','', opt$outcome)}/time_restricted_phenomes/"),
              recursive = TRUE )
}
### save ukb matching data
fwrite(ukb_post_match_cov,
       glue("data/private/ukb/{opt$ukb_version}/",
            "X{gsub('X', '', opt$outcome)}/matched_covariates.txt"),
       sep = "\t")

# 8. create time-restricted phenomes -------------------------------------------
## mgi
cli_alert_info("constructing time-restricted phenomes in mgi...")
mgi_matched_phe <- merge.data.table(
  mgi_first_phe[id %in% mgi_post_match_cov[, id]],
  mgi_post_match_cov,
  by = "id"
)
mgi_pims <- list()
for (i in seq_along(time_thresholds)) {
  mgi_pims[[i]] <- generate_restricted_phenome(phe_data    = mgi_matched_phe,
                                               threshold   = time_thresholds[i],
                                               cases       = mgi_case_ids,
                                               outcome_phe = opt$outcome)
}
names(mgi_pims) <- glue("t{time_thresholds}")
for (i in 1:length(mgi_pims)) {
  fwrite(
    x = mgi_pims[[i]],
    file = glue("data/private/mgi/{opt$mgi_version}/",
                "X{gsub('X', '', opt$outcome)}/time_restricted_phenomes/",
                "mgi_X{gsub('X', '', opt$outcome)}",
                "_{names(mgi_pims)[i]}_{opt$mgi_version}.txt"),
    sep = "\t"
  )
}

## ukb
cli_alert_info("constructing time-restricted phenomes in ukb...")
ukb_matched_phe <- merge.data.table(
  ukb_first_phe[id %in% ukb_post_match_cov[, id]],
  ukb_post_match_cov,
  by = "id"
)
ukb_pims <- list()
for (i in seq_along(time_thresholds)) {
  ukb_pims[[i]] <- generate_restricted_phenome(phe_data    = ukb_matched_phe,
                                               threshold   = time_thresholds[i],
                                               cases       = ukb_case_ids,
                                               outcome_phe = opt$outcome)
}
names(ukb_pims) <- glue("t{time_thresholds}")
for (i in 1:length(ukb_pims)) {
  fwrite(
    x = ukb_pims[[i]],
    file = glue("data/private/ukb/{opt$ukb_version}/",
                "X{gsub('X', '', opt$outcome)}/time_restricted_phenomes/",
                "ukb_X{gsub('X', '', opt$outcome)}",
                "_{names(ukb_pims)[i]}_{opt$ukb_version}.txt"),
    sep = "\t"
  )
}
