### libraries ------------------------------------------------------------------
options(stringsAsFactors = FALSE)
library(data.table)
library(purrr)
library(fst)
library(MatchIt)
library(optparse)

set.seed(61787)

### parser ---------------------------------------------------------------------
option_list <- list(
  make_option("--cohort", type="character", default="",
              help="cohort - 'mgi' or 'ukb'"),
  make_option("--mgi_version", type="character", default="20210318",
              help="data version for mgi data (e.g., '20210318')"),
  make_option("--ukb_version", type="character", default="20221117",
              help="data version for ukb data (e.g., '20221117')"),
  make_option("--matching_ratio", type="character", default="2",
              help="number of noncases to select per case (e.g., '2')"),
  make_option("--matching_covariates", type="character", default="age_at_first_diagnosis, length_followup, female",
              help="covariates on which to match cases to non-cases (e.g., 'age_at_first_diagnosis, length_followup, female')"),
  make_option("--outcome", type="character", default="",
              help="outcome phecode (e.g., '157')"),
  make_option("--time_thresholds", type="character", default="0, 1, 2, 3, 5",
              help="one or multiple time-thresholds (in years)")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options

if(opt$cohort == "" | opt$outcome == "") {
  print("Not all parameters present")
}

### specs ----------------------------------------------------------------------
# cohort <- "mgi"
cohort <- "ukb"
# if (opt$mgi_version == "") { mgi_version     <- "20210318" } # mgi data version
mgi_version     <- "20210318"       # mgi phenome version
ukb_version     <- "20221117"       # ukb phenome version
outcome         <- "157"            # outcome phecode (157 - PanCan, 155 - LivCan, 184.1 - OvCan)
time_thresholds <- c(0, 1, 2, 3, 5) # time thresholds
matching_ratio  <- 2                # number of noncases per case

### load functions -------------------------------------------------------------
for (i in list.files("fn/")) {source(paste0("fn/", i))}

### prep -----------------------------------------------------------------------

## confirm file structure for a given outcome exists - if not, create paths
check_folder_structure(cohort = cohort, data_version = ifelse(tolower(cohort) == "mgi", mgi_version, ukb_version), outcome_phecode = outcome)

## pull file paths corresponding to the data version specified
file_paths <- get_files()

### load data ------------------------------------------------------------------

## phecode indicator matrix (PEDMASTER_0)
pim0 <- fread(file_paths[[cohort]]$pim0_file)
setnames(pim0, old = c("IID", paste0("X", outcome)), new = c("id", "outcome"))

## demographics data
if (tolower(cohort) == "mgi") {
  demo <- data.table::fread(file_paths[[cohort]]$demo_file)[, .(id = Deid_ID, age = Age, alive = as.numeric(AliveYN == "Y"), dead_dsb = Deceased_DaysSinceBirth, ethn = EthnicityName, marital = MaritalStatusCode, sex = Sex, race = RaceName)]
}
if (tolower(cohort) == "ukb") {
  demo <- data.table::fread(file_paths[[cohort]]$demo_file, na.strings = c(""))[, .(id, ethn = ethnicity, sex)]
}

## phecode data (first occurrence of each phecode)
## see 'get_full_icd' for all time-stamps
first_by_phe <- get_first_phecode(cohort      = cohort,
                                  mgi_version = mgi_version,
                                  ukb_file = file_paths[[cohort]]$icd_phecode_file)
first_by_phe[, dsb := as.numeric(dsb)]

## identify cases
first_case <- unique(first_by_phe[phecode == outcome])
case_ids   <- first_case[, unique(id)]

### do matching ----------------------------------------------------------------
# get_diag_metrics is a function that contains variables regarding length of
# follow-up, number of encounters, density of encounters, and other potentially
# relevant variables
diag_metrics <- get_diag_metrics(cohort = cohort, icd_phecode_data = first_by_phe, data_version = mgi_version)

## merge diagnostic metrics into demographic table
matching_cov <- merge.data.table(
  demo[, .(id, age, female = as.numeric(sex == "F"))],
  diag_metrics,
  by = "id"
)[, case := fifelse(id %in% case_ids, 1, 0)][]

## perform matching
matched <- lets_match(
  dataset      = "matching_cov",
  nearest_vars = c("age_at_first_diagnosis", "length_followup"),
  exact_vars   = c("female"),
  calip        = 0.25,
  r            = matching_ratio
)

sum_matched <- summary(matched)
my_matched  <- match.data(matched)

post_match_cov <- merge.data.table(
  my_matched,
  merge.data.table(
    first_case,
    my_matched,
    by = "id",
    all.x = TRUE
    )[, .(subclass, case_dsb = dsb)],
  by = "subclass"
)

# create DSB-based thresholds for thresholding phenome data
for (i in time_thresholds) {
  post_match_cov[, (paste0("t", i, "_threshold")) := floor(case_dsb - (365.25 * i))]
}

## create indicators on whether to keep a matched group at each threshold
# if the threshold precedes or occurs on the same date as the first recorded
# diagnosis of any patient in the matched group,
# remove the entire group (0), otherwise keep it (1)
for (i in time_thresholds) {
  post_match_cov[, (paste0("t", .x, "_indicator")) := as.numeric(all(get(paste0("t", .x, "_threshold")) > first_dsb)), by = subclass]
}

# n at each time threshold
map_df(post_match_cov[, (paste0("t", time_thresholds, "_indicator")), with = FALSE], ~sum(.x))
# n_cases at each time threshold
map_df(post_match_cov[case == 1][, (paste0("t", time_thresholds, "_indicator")), with = FALSE], ~sum(.x))

if (!dir.exists( paste0("data/", mgi_version, "/processed/X", gsub("X", "", outcome), "/"))) {
  dir.create( paste0("data/", mgi_version, "/processed/X", gsub("X", "", outcome), "/time_restricted_phenomes/"), recursive = TRUE)
}

### save matched covariate data
data.table::fwrite(post_match_cov,
       paste0("data/", mgi_version, "/processed/X", gsub("X", "", outcome), "/matched_covariates.txt"),
       sep = "\t")

### create restricted phenomes -------------------------------------------------
matched_icd <- data.table::merge.data.table(
  first_by_phe[id %in% post_match_cov[, id]],
  post_match_cov,
  by.x = "id",
  by.y = "id"
)

generate_restricted_phenome <- function(icd_data, threshold, cases, outcome_phe) {
  
  if (!grepl("X", outcome_phe)) {
    outcome_phe <- paste0("X", outcome_phe)
  }
  
  out <- data.table::dcast(
    unique(icd_data[get(paste0("t", threshold, "_indicator")) == 1][dsb < get(paste0("t", threshold, "_threshold"))][, .(id, phecode = paste0("X", phecode))]),
    id ~ phecode,
    value.var = "phecode",
    fun.aggregate = length,
    fill = 0
  )
  
  out[, case := data.table::fifelse(id %in% cases, 1, 0)]
  if (outcome_phe %in% names(out)) {
    out[, c(outcome_phe) := NULL]
  }
  
  return(out)
  
}


# list object of time-restricted phecode indicator matrices
tr_pims <- map(time_thresholds,
    ~generate_restricted_phenome(icd_data = matched_icd, threshold = .x, cases = case_ids, outcome_phe = outcome))
names(tr_pims) <- paste0("t", time_thresholds)


# save phecode indicator matrices ----------------------------------------------
output_tsv <- function(data, names, ver = NULL, outc = NULL){ 
  if (is.null(ver)) {stop("Missing version (`ver`) argument. Check is `version` is specified in environment.")}
  if (is.null(outc)) {stop("Missing outcome (`outc`) arguement. Check if `outcome` is specified in environmnet.")}
  fwrite(data, paste0("./data/", ver, "/processed/X", gsub("X", "", outc), "/time_restricted_phenomes/mgi_X", gsub("X", "", outc), "_", names, "_", ver, ".txt"), sep = "\t")
}

list(
  data = tr_pims,
  names = names(tr_pims)
) |>
  purrr::pmap(~output_tsv(data = .x, names = .y, ver = mgi_version, outc = outcome))
