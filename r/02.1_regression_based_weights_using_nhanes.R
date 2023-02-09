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

isource <- function(fn_path) {
  lapply(list.files(fn_path, full.names = TRUE), source) |> invisible()
}
isource("fn")


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
    quant2.5 <- quantile(x, probs = 0.025, na.rm = TRUE)
    quant97.5 <- quantile(x, probs = 0.975, na.rm = TRUE)
    x[x < quant2.5] <- quant2.5
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

# save -------------------------------------------------------------------------
write_fst(
  x = merged[, .(id = DeID_PatientID, no_cancer_weights = WEIGHT_NHANES_NOCAN, cancer_indirect_weights = CANCER_NHANES_UNCORRECTED)],
  path = glue("{data_path}weights_{opt$cohort_version}_{opt$mgi_cohort}.fst")
  )

cli_alert_success("script success! see {.path {data_path}} and suffix {.emph {opt$mgi_cohort}}")
