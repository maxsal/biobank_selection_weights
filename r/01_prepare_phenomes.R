# prepare mgi and ukb data for time-restricted phecode-phecode phewas
# outputs: time-threshold phecode indicator matrices
# author:  max salvatore
# date:    20221207

# 1. libraries, functions, and options -----------------------------------------
options(stringsAsFactors = FALSE)

library(data.table)
library(MatchIt)
library(cli)

set.seed(61787)

for (i in list.files("fn/")) source(paste0("fn/", i)) # load functions

# 2. specifications ------------------------------------------------------------
cohort          <- "ukb"
mgi_version     <- "20210318"       # mgi phenome version
ukb_version     <- "20221117"       # ukb phenome version
outcome         <- "157"            # outcome phecode (157 - PanCan, 155 - LivCan, 184.1 - OvCan)
time_thresholds <- c(0, 1, 2, 3, 5) # time thresholds
matching_ratio  <- 2                # number of noncases per case


# 3. extra preparations --------------------------------------------------------
## confirm file structure for a given outcome exists - if not, create paths
### mgi
check_folder_structure(
  cohort          = "mgi",
  data_version    = mgi_version,
  outcome_phecode = outcome
  )

### ukb
check_folder_structure(
  cohort          = "ukb",
  data_version    = ukb_version,
  outcome_phecode = outcome
)

## pull file paths corresponding to the data version specified
file_paths <- get_files()

# 4. read data -----------------------------------------------------------------
## mgi
cli::cli_alert_info("loading mgi data...")
### phecode indicator matrix (PEDMASTER_0)
mgi_pim0 <- fread(file_paths[["mgi"]]$pim0_file)
data.table::setnames(mgi_pim0,
                     old = c("IID", paste0("X", outcome)),
                     new = c("id", "outcome"))

### demographics
mgi_demo <- data.table::fread(file_paths[[cohort]]$demo_file)[, .(
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
cli::cli_alert_info("loading ukb data...")
### phecode indicator matrix (PEDMASTER_0)
ukb_pim0 <- fread(file_paths[["ukb"]]$pim0_file)
data.table::setnames(ukb_pim0,
                     old = c("IID", paste0("X", outcome)),
                     new = c("id", "outcome"))
### demographics
ukb_demo <- data.table::fread(file_paths[["ukb"]]$demo_file,
                              na.strings = c(""))[, .(
                                id,
                                ethn = ethnicity,
                                sex)]

### icd-phecode data
ukb_first_phe <- get_full_icd_dsb_phecode(
  icd_file      = file_paths[["ukb"]]$icd_phecode_file
)

# 5. identify cases -------------------------------------------------------------
## mgi
mgi_case <- unique(mgi_first_phe[phecode == outcome])
mgi_case_ids <- mgi_case[, unique(id)]

## ukb
ukb_case <- unique(ukb_first_phe[phecode == outcome])
ukb_case_ids <- ukb_case[, unique(id)]

# 6. calculate diagnostic metrics -----------------------------------------------
## mgi
cli::cli_alert_info("calculating diagnostic metrics in mgi...")
mgi_diag_metrics <- get_icd_phecode_metrics(
  full_phe_data = mgi_full_phe
)
mgi_matching_cov <- merge.data.table(
  mgi_demo[, .(id, age, female = as.numeric(sex == "F"))],
  mgi_diag_metrics,
  by = "id"
)[, case := fifelse(id %in% mgi_case_ids, 1, 0)]

## ukb
cli::cli_alert_info("calculating diagnostic metrics in ukb...")
ukb_diag_metrics <- get_icd_phecode_metrics(
  full_phe_data = ukb_first_phe
)
ukb_matching_cov <- merge.data.table(
  ukb_demo[, .(id, age, female = as.numeric(sex == "F"))]
)[, case := fifelse(id %in% ukb_case_ids, 1, 0)]

# 7. perform matching ----------------------------------------------------------
## mgi
cli::cli_alert_info("performing 1:{matching_ratio} case:non-case matching in mgi...")
mgi_match_text <- glue::glue("MatchIt::matchit(case ~ ",
                             "{paste0(c(nearest_matching_vars, exact_matching_vars), collapse = ' + ')}, ", 
                             "data = mgi_matching_cov, calclosest = TRUE, ",
                             "mahvars = c({paste0(sapply(nearest_matching_vars, function(x) paste0('\\'', x, '\\'')), collapse = ', ')}), ",
                             "caliper = {matching_caliper}, ",
                             "exact = c({paste0(sapply(exact_matching_vars, function(x) paste0('\\'', x, '\\'')), collapse = ', ')}), ",
                             "ratio = {matching_ratio})")
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
  mgi_post_match_cov[, (glue::glue("t{i}_threshold")) := floor(case_dsb - (365.25 *i))]
}
for (i in time_thresholds) {
  mgi_post_match_cov[, (glue::glue("t{i}_indicator")) := as.numeric(all(get(glue::glue("t{i}_threshold")) > first_dsb)), by = subclass]
}
if ( !dir.exists( glue::glue("data/private/mgi/{mgi_version}/X{gsub('X','', outcome)}/time_restricted_phenomes/") ) ) {
  dir.create( glue::glue("data/private/mgi/{mgi_version}/X{gsub('X','', outcome)}/time_restricted_phenomes/"), recursive = TRUE )
}
### save mgi matching data
data.table::fwrite(mgi_post_match_cov,
                   glue::glue("data/private/mgi/{mgi_version}/",
                              "X{gsub('X', '', outcome)}/matched_covariates.txt"),
                   sep = "\t")

## ukb
cli::cli_alert_info("performing 1:{matching_ratio} case:non-case matching in ukb...")
ukb_match_text <- glue::glue("MatchIt::matchit(case ~ ",
                             "{paste0(c(nearest_matching_vars, exact_matching_vars), collapse = ' + ')}, ", 
                             "data = ukb_matching_cov, calclosest = TRUE, ",
                             "mahvars = c({paste0(sapply(nearest_matching_vars, function(x) paste0('\\'', x, '\\'')), collapse = ', ')}), ",
                             "caliper = {matching_caliper}, ",
                             "exact = c({paste0(sapply(exact_matching_vars, function(x) paste0('\\'', x, '\\'')), collapse = ', ')}), ",
                             "ratio = {matching_ratio})")
ukb_match <- eval(parse(text = ukb_match_text))
ukb_matched <- MatchIt::match.data(ukb_match)
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
  ukb_post_match_cov[, (glue::glue("t{i}_threshold")) := floor(case_dsb - (365.25 *i))]
}
for (i in time_thresholds) {
  ukb_post_match_cov[, (glue::glue("t{i}_indicator")) := as.numeric(all(get(glue::glue("t{i}_threshold")) > first_dsb)), by = subclass]
}
if ( !dir.exists( glue::glue("data/private/ukb/{ukb_version}/X{gsub('X','', outcome)}/time_restricted_phenomes/") ) ) {
  dir.create( glue::glue("data/private/ukb/{ukb_version}/X{gsub('X','', outcome)}/time_restricted_phenomes/"), recursive = TRUE )
}
### save ukb matching data
data.table::fwrite(ukb_post_match_cov,
                   glue::glue("data/private/ukb/{ukb_version}/",
                              "X{gsub('X', '', outcome)}/matched_covariates.txt"),
                   sep = "\t")

# create time-restricted phenomes ----------------------------------------------
## mgi
cli::cli_alert_info("constructing time-restricted phenomes in mgi...")
mgi_matched_phe <- data.table::merge.data.table(
  mgi_first_phe[id %in% mgi_post_match_cov[, id]],
  mgi_post_match_cov,
  by = "id"
)
mgi_pims <- list()
for (i in seq_along(time_thresholds)) {
  mgi_pims[[i]] <- generate_restricted_phenome(phe_data    = mgi_matched_phe,
                                               threshold   = time_thresholds[i],
                                               cases       = mgi_case_ids,
                                               outcome_phe = outcome)
}
names(mgi_pims) <- glue::glue("t{time_thresholds}")
for (i in 1:length(mgi_pims)) {
  data.table::fwrite(
    x = mgi_pims[[i]],
    file = glue::glue("data/private/mgi/{version}/X{gsub('X', '', outcome)}",
                      "/time_restricted_phenomes/mgi_X{gsub('X', '', outcome)}",
                      "_{names(tr_pims)[i]}_{version}.txt"),
    sep = "\t"
  )
}

## ukb
cli::cli_alert_info("constructing time-restricted phenomes in ukb...")
ukb_matched_phe <- data.table::merge.data.table(
  ukb_first_phe[id %in% ukb_post_match_cov[, id]],
  ukb_post_match_cov,
  by = "id"
)
ukb_pims <- list()
for (i in seq_along(time_thresholds)) {
  ukb_pims[[i]] <- generate_restricted_phenome(phe_data    = ukb_matched_phe,
                                               threshold   = time_thresholds[i],
                                               cases       = ukb_case_ids,
                                               outcome_phe = outcome)
}
names(ukb_pims) <- glue::glue("t{time_thresholds}")
for (i in 1:length(ukb_pims)) {
  data.table::fwrite(
    x = ukb_pims[[i]],
    file = glue::glue("data/private/ukb/{version}/X{gsub('X', '', outcome)}",
                      "/time_restricted_phenomes/ukb_X{gsub('X', '', outcome)}",
                      "_{names(tr_pims)[i]}_{version}.txt"),
    sep = "\t"
  )
}




