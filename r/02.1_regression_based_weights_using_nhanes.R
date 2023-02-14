# simplex regression-based weights
# requires: mgi data; nhanes wave info
# outputs:  nhanes and uncorrected cancer weights
# author:   max salvatore
# date:     20230203

# libraries --------------------------------------------------------------------
library(haven)
library(survey)
library(dplyr)
library(pracma)
library(simplexreg)
library(data.table)
library(glue)
library(fst)
library(optparse)

# optparse list ---
option_list <- list(
  make_option("--cohort_version", type = "character", default = "20220822",
              help = "Cohort version in /net/junglebook/magic_data/EHRdata/ [default = '20220822']"),
  make_option("--mgi_cohort", type = "character", default = "comb",
              help = "Specific MGI cohort [default = %default]"),
  make_option("--nhanes_wave_letter", type = "character",
              default = "J",
              help = glue("NHANES data prefix corresponding to wave ",
                          "[default = %default]")),
  make_option("--nhanes_wave_years", type = "character",
              default = "2017-2018",
              help = glue("NHANES wave years corresponding to wave ",
                          "[default = %default]")),
  make_option("--nhanes_survey_names", type = "character",
              default = "DEMO,BMX,SMQ,DIQ,MCQ",
              help = glue("NHANES wave years corresponding to wave ",
                          "[default = %default]"))
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args   <- parse_args(parser, positional_arguments = 0)
opt    <- args$options

data_path <- glue("data/private/mgi/{opt$cohort_version}/")

lapply(list.files("fn", full.names = TRUE), source) |> # load functions
  invisible()



# load data --------------------------------------------------------------------
mgi <- read_fst(glue("{data_path}data_{opt$cohort_version}_{opt$mgi_cohort}.fst"),
               as.data.table = TRUE)

nhanes_datasets <- unlist(strsplit(opt$nhanes_survey_names, ","))
nhanes_merged <- download_nhanes_data(
  wave_letter = opt$nhanes_wave_letter,
  wave_years  = opt$nhanes_wave_years,
  datasets    = nhanes_datasets
)

keep_vars <- c("SEQN", "RIAGENDR", "WTINT2YR", "RIDAGEYR", "RIDRETH1", "MCQ220",
               "BMXBMI", "SMQ040", "SMQ020", "DIQ010", "MCQ160C", "WTMEC2YR",
               "SDMVSTRA", "SDMVPSU")

if ("WTMECPRP" %in% names(nhanes_merged)) {
  setnames(nhanes_merged,
           "WTMECPRP",
           "WTMEC2YR")
}
if ("WTINTPRP" %in% names(nhanes_merged)) {
  setnames(nhanes_merged,
           "WTINTPRP",
           "WTINT2YR")
}

nhanes_merged <- nhanes_merged[, ..keep_vars]

prepped_nhanes <- prepare_nhanes_data(
  nhanes_data = nhanes_merged,
  mec_wt_var = "WTMEC2YR")

stacked <- rbindlist(list(
  prepped_nhanes,
  mgi[][, dataset := "MGI"]
), use.names = TRUE, fill = TRUE)

# weight estimation function ---------------------------------------------------
# ADAPTED FROM: /net/junglebook/home/kundur/EHR/Processed Code/Weighted_using_lauren_code_bb.R
lauren_nhanes <- function(stacked_data) {
  
  stacked_data[dataset == "NHANES", weight_nhanes := .N * weight_nhanes / sum(weight_nhanes, na.rm = TRUE)]

  selection_NHANES_NOCAN <- simplexreg(samp_nhanes ~ as.numeric(age_cat == 5) + as.numeric(age_cat == 6) + cad + diabetes + smoking_current + smoking_former + bmi_under + bmi_overweight + bmi_obese + nhanes_nhw, data = stacked_data[dataset == 'NHANES', ])
  
  mgiselect_NOCAN <- glm(as.numeric(dataset == "MGI") ~ as.numeric(age_cat == 5) + as.numeric(age_cat == 6) + diabetes + cad + bmi_under + bmi_overweight + bmi_obese + smoking_current + smoking_former + nhanes_nhw, data = stacked_data, family = quasibinomial())
  
  p_Sext <- predict(selection_NHANES_NOCAN, newdata = stacked_data[dataset == "MGI", ], type = "response")[, 1]
  p_MGI  <- predict(mgiselect_NOCAN, newdata = stacked_data[dataset == "MGI", ], type = "response")
  temp <- rep(0, times = length(p_MGI))
  temp[which(rownames(data.frame(p_Sext)) %in% rownames(data.frame(p_MGI)) == T)] <- p_Sext
  temp[which(rownames(data.frame(p_Sext)) %in% rownames(data.frame(p_MGI)) == F)] <- NA
  p_Sext <- temp
  p_Sext[which(p_Sext == 0)] <- 1.921e-05
  SELECT_NHANES_NOCAN <- p_Sext * (p_MGI / (1 - p_MGI))
  
  chopr <- function(x) {
    quant2.5  <- quantile(x, probs = 0.025, na.rm = TRUE)
    quant97.5 <- quantile(x, probs = 0.975, na.rm = TRUE)
    x[x < quant2.5]  <- quant2.5
    x[x > quant97.5] <- quant97.5
    return(x)
  }
  
  SELECT_NHANES_NOCAN <- chopr(SELECT_NHANES_NOCAN)
  WEIGHT_NHANES_NOCAN <- 1 / SELECT_NHANES_NOCAN
  WEIGHT_NHANES_NOCAN <- stacked_data[dataset == "MGI", .N] * WEIGHT_NHANES_NOCAN / sum(WEIGHT_NHANES_NOCAN, na.rm = TRUE)

    ## With Cancer
  NHANES_cancer_model <- glm(cancer ~ as.numeric(age_cat == 5) + as.numeric(age_cat == 6) + diabetes + cad + bmi_under + bmi_overweight + bmi_obese + smoking_current + smoking_former + nhanes_nhw, data = stacked_data[dataset == "NHANES", ], weights = weight_nhanes, family = quasibinomial())
  
  modelCAN_NHANES <- predict(NHANES_cancer_model, type = "response", newdata = stacked_data)
  
  MGI_cancer_model <- glm(cancer ~ as.numeric(age_cat == 5) + as.numeric(age_cat == 6) + diabetes + cad + bmi_under + bmi_overweight + bmi_obese + smoking_current + smoking_former + nhanes_nhw, data = stacked_data[dataset == "MGI", ], family = quasibinomial())
  
  modelCAN_MGI <- predict(MGI_cancer_model, type = "response", newdata = stacked_data)
  denom <- ifelse(stacked_data[, cancer] == 1, modelCAN_MGI, 1 - modelCAN_MGI)
  num <- ifelse(stacked_data[, cancer] == 1, modelCAN_NHANES, 1 - modelCAN_NHANES)
  cancer_factor <- (num[stacked_data[, dataset] == "MGI"] / denom[stacked_data[, dataset] == "MGI"])
  cancer_factor <- chopr(cancer_factor)
  CANCER_NHANES_UNCORRECTED = cancer_factor * WEIGHT_NHANES_NOCAN
  CANCER_NHANES_UNCORRECTED  = stacked_data[dataset == "MGI", .N] * CANCER_NHANES_UNCORRECTED / sum(CANCER_NHANES_UNCORRECTED, na.rm = TRUE)
  
  data.table(
    "WEIGHT_NHANES_NOCAN"       = WEIGHT_NHANES_NOCAN,
    "CANCER_NHANES_UNCORRECTED" = CANCER_NHANES_UNCORRECTED
    )

}

estimated_weights <- lauren_nhanes(stacked_data = stacked)

merged <- cbind(mgi, estimated_weights)

# estimate cancer~female log(OR) -----------------------------------------------
m0 <- glm(cancer~as.numeric(Sex == "F"), family = quasibinomial(), data = merged)
m1 <- glm(cancer~as.numeric(Sex == "F"), family = quasibinomial(), data = merged, weights = WEIGHT_NHANES_NOCAN)
m2 <- glm(cancer~as.numeric(Sex == "F"), family = quasibinomial(), data = merged, weights = CANCER_NHANES_UNCORRECTED)

extractr <- function(x, weight_name) {
  suppressMessages({y <- confint(x)})
  data.table(
    "weights" = weight_name,
    "est"     = coef(x)[[2]],
    "lower"   = y[2, 1],
    "upper"   = y[2, 2],
    "var"     = diag(summary(x)$cov.scaled)[[2]]
  )
}

# summarize estimates cancer~female log(OR)
rbindlist(list(
  extractr(x = m0, weight_name = "None"),
  extractr(x = m1, weight_name = "No cancer"),
  extractr(x = m2, weight_name = "Cancer (uncorrected)")
)) |>
  fwrite(file = glue("{data_path}cancer_female_logor_est_{opt$cohort_version}_{opt$mgi_cohort}.csv"))

## poststratification weights --------------------------------------------------
poststratification <- function(
    mgi_data,
    last_entry_age_var = "AgeLastEntry",
    cancer_var         = "cancer",
    chd_var            = "cad",
    smoke_var          = "smoking_current",
    diabetes_var       = "diabetes"
    ) {
  
  Nobs      <- nrow(mgi_data)
  age_var_2 <- as.numeric(mgi_data[[last_entry_age_var]])
  ages      <- seq(18, 85, 1)
  which_between <- function(vec, mat) {
    between(x = vec, lower = mat[["lower"]], upper = mat[["upper"]])
  }
  
  # cancer prevalence by age (US) 
  # https://seer.cancer.gov/csr/1975_2016/results_merged/topic_prevalence.pdf
  cancer_prevalence <- age_grp_table(
    lower_ages   = c(0, 10, 20, 30, 40, 50, 60, 70, 80),
    num_vec      = c(0.0899, 0.2023, 0.3922, 0.8989, 2.1532, 4.9326, 10.4420, 18.3168, 21.5939) / 100,
    num_var_name = "prevalence"
  )
  
  cancer_func_pop <- stepfun(x = cancer_prevalence[["lower"]], y = c(0, cancer_prevalence[["prevalence"]]), right = FALSE)
  b <- aggregate(as.formula(paste0(cancer_var, " ~ ", last_entry_age_var)), FUN = mean, data = mgi_data)
  cancer_func_mgi <- stepfun(x = b[, 1], y = c(0, b[, 2]), right = FALSE)
  
  # diabetes prevalence by age (US)
  # https://www.cdc.gov/diabetes/pdfs/data/statistics/national-diabetes-statistics-report.pdf 
  # youngest group will have no people, no number provided in documentation
  diabetes_prevalence <- age_grp_table(
    lower_ages   = c(0, 18, 45, 65),
    num_vec      = c(0.01, 0.030, 0.138, 0.214),
    num_var_name = "prevalence"
  )
  
  diabetes_func_pop <- stepfun(x = diabetes_prevalence[["lower"]], y = c(0, diabetes_prevalence[["prevalence"]]), right = FALSE)
  d_b <- aggregate(as.formula(paste0(diabetes_var, " ~ ", last_entry_age_var)), FUN = mean, data = mgi_data)
  diabetes_func_mgi <- stepfun(x = d_b[, 1], y = c(0, d_b[, 2]), right = FALSE)
  
  # chd prevalence by age (us)
  # https://www.cdc.gov/nchs/fastats/heart-disease.htm
  # youngest group will have no people, no number provided in documentation
  chd_prevalence <- age_grp_table(
    lower_ages   = c(0,18,45,65,75),
    num_vec      = c(0.01, 0.01, 0.060, 0.155, 0.239),
    num_var_name = "prevalence"
  )
  
  chd_func_pop <- stepfun(x = chd_prevalence[["lower"]], y = c(0, chd_prevalence[["prevalence"]]), right = FALSE)
  chd_b <- aggregate(as.formula(paste0(chd_var, " ~ ", last_entry_age_var)), FUN = mean, data = mgi_data)
  chd_func_mgi <- stepfun(x = chd_b[, 1], y = c(0, chd_b[, 2]), right = FALSE)

  # smoking prevalence by age (us)
  # https://www.cdc.gov/nchs/fastats/heart-disease.htm
  # youngest group will have no people, no number provided in documentation
  smoke_prevalence <- age_grp_table(
    lower_ages   = c(0,18,25,45,65),
    num_vec      = c(0.01, 0.074, 0.141, 0.149, 0.09),
    num_var_name = "prevalence"
  )
  
  smoke_func_pop <- stepfun(x = smoke_prevalence[["lower"]], y = c(0, smoke_prevalence[["prevalence"]]), right = FALSE)
  smoke_b        <- aggregate(as.formula(paste0(smoke_var, " ~ ", last_entry_age_var)), FUN = mean, data = mgi_data)
  smoke_func_mgi <- stepfun(x = smoke_b[, 1], y = c(0, smoke_b[, 2]), right = FALSE)
  
  # age distribution (us)
  # https://www.census.gov/data/tables/2000/dec/phc-t-09.html
  total_population  <- 281421906
  male_population   <- 138053563
  female_population <- total_population - male_population
  
  low_ages <- seq(0, 85, 5)
  age_counts_male <- age_grp_table(
    lower_ages   = low_ages,
    num_vec      = c(9810733,10523277,10520197,10391004,9687814,9798760,10321769,11318696,
                11129102,9889506,8607724,6508729,5136627,4400362,3902912,3044456,1834897,1226998),
    num_var_name = "counts"
  )
  age_counts_female <- age_grp_table(
    lower_ages   = low_ages,
    num_vec      = c(9365065,10026228,10007875,9828886,9276187,9582576,10188619,11387968,11312761,10202898,
                     8977824,6960508,5668820,5133183,4954529,4371357,3110470,3012589),
    num_var_name = "counts"
  )
  age_prevalence <- age_grp_table(
    lower_ages   = low_ages,
    num_vec      = (age_counts_male[["counts"]] + age_counts_female[["counts"]]) / total_population,
    num_var_name = "prevalence"
  )
  age_func_pop <- stepfun(x = age_prevalence[["lower"]], y = c(0, age_prevalence[["prevalence"]]), right = FALSE)
  age_func_mgi <- stepfun(x = age_prevalence[["lower"]], y = as.numeric(c(0, table(cut(age_var_2, breaks = c(0, age_prevalence[["upper"]]), labels = age_prevalence[["group"]])) / Nobs)), right = FALSE)
  
  # with cancer
  population <- fifelse(mgi_data[[cancer_var]] == 1, cancer_func_pop(age_var_2), 1 - cancer_func_pop(age_var_2))
  population <- population * fifelse(mgi_data[[diabetes_var]] == 1, diabetes_func_pop(age_var_2), 1 - diabetes_func_pop(age_var_2))
  population <- population * fifelse(mgi_data[[chd_var]] == 1, chd_func_pop(age_var_2), 1 - chd_func_pop(age_var_2))
  population <- population * fifelse(mgi_data[[smoke_var]] == 1, smoke_func_pop(age_var_2), 1 - smoke_func_pop(age_var_2))
  population <- population * age_func_pop(age_var_2)

  mgi <- fifelse(mgi_data[[cancer_var]] == 1, cancer_func_mgi(age_var_2), 1 - cancer_func_mgi(age_var_2))  
  mgi <- mgi * fifelse(mgi_data[[diabetes_var]] == 1, diabetes_func_mgi(age_var_2), 1 - diabetes_func_mgi(age_var_2))
  mgi <- mgi * fifelse(mgi_data[[chd_var]] == 1, chd_func_mgi(age_var_2), 1 - chd_func_mgi(age_var_2))
  mgi <- mgi * fifelse(mgi_data[[smoke_var]] == 1, smoke_func_mgi(age_var_2), 1 - smoke_func_mgi(age_var_2))
  mgi <- mgi * age_func_mgi(age_var_2)
  
  post <- population / mgi
  post <- (Nobs * post) / sum(post, na.rm = TRUE)
  
  return(
    data.table(
      cancer_post_uncorrected = post
    )
  )

}

poststratification(mgi_data = mgi)

tmp <- cbind(mgi_data, post)

glm(cancer~as.numeric(Sex == "F"), family = quasibinomial(), data = tmp) |> coef()
glm(cancer~as.numeric(Sex == "F"), family = quasibinomial(), data = tmp, weights = post) |> coef()

# save -------------------------------------------------------------------------
write_fst(
  x = merged[, .(id = DeID_PatientID, no_cancer_weights = WEIGHT_NHANES_NOCAN, cancer_indirect_weights = CANCER_NHANES_UNCORRECTED)],
  path = glue("{data_path}weights_{opt$cohort_version}_{opt$mgi_cohort}.fst")
  )

cli_alert_success("script success! see {.path {data_path}} and suffix {.emph {opt$mgi_cohort}}")
