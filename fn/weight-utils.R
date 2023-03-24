# calculate IPW from stacked data
# ADAPTED FROM: /net/junglebook/home/kundur/EHR/Processed Code/Weighted_using_lauren_code_bb.R
ipw <- function(stacked_data, dataset_name = "MGI", id_var = "id") {
  
  stacked_data[dataset == "NHANES", weight_nhanes := .N * weight_nhanes / sum(weight_nhanes, na.rm = TRUE)]
  
  selection_NHANES_NOCAN <- simplexreg(samp_nhanes ~ as.numeric(age_cat == 5) + as.numeric(age_cat == 6) + cad + diabetes + smoking_current + smoking_former + bmi_under + bmi_overweight + bmi_obese + nhanes_nhw + female, data = stacked_data[dataset == 'NHANES', ])
  
  mgiselect_NOCAN <- glm(as.numeric(dataset == dataset_name) ~ as.numeric(age_cat == 5) + as.numeric(age_cat == 6) + diabetes + cad + bmi_under + bmi_overweight + bmi_obese + smoking_current + smoking_former + nhanes_nhw + female, data = stacked_data, family = quasibinomial())
  
  p_Sext <- predict(selection_NHANES_NOCAN, newdata = stacked_data[dataset == dataset_name, ], type = "response")[, 1]
  p_MGI  <- predict(mgiselect_NOCAN, newdata = stacked_data[dataset == dataset_name, ], type = "response")
  temp <- rep(0, times = length(p_MGI))
  temp[which(rownames(data.frame(p_Sext)) %in% rownames(data.frame(p_MGI)) == T)] <- p_Sext
  temp[which(rownames(data.frame(p_Sext)) %in% rownames(data.frame(p_MGI)) == F)] <- NA
  p_Sext <- temp
  p_Sext[which(p_Sext == 0)] <- 1.921e-05
  SELECT_NHANES_NOCAN <- p_Sext * (p_MGI / (1 - p_MGI))
  
  SELECT_NHANES_NOCAN <- chopr(SELECT_NHANES_NOCAN)
  WEIGHT_NHANES_NOCAN <- 1 / SELECT_NHANES_NOCAN
  WEIGHT_NHANES_NOCAN <- stacked_data[dataset == dataset_name, .N] * WEIGHT_NHANES_NOCAN / sum(WEIGHT_NHANES_NOCAN, na.rm = TRUE)
  
  ## With Cancer
  NHANES_cancer_model <- glm(cancer ~ as.numeric(age_cat == 5) + as.numeric(age_cat == 6) + diabetes + cad + bmi_under + bmi_overweight + bmi_obese + smoking_current + smoking_former + nhanes_nhw + female, data = stacked_data[dataset == "NHANES", ], weights = weight_nhanes, family = quasibinomial())
  
  modelCAN_NHANES <- predict(NHANES_cancer_model, type = "response", newdata = stacked_data)
  
  MGI_cancer_model <- glm(cancer ~ as.numeric(age_cat == 5) + as.numeric(age_cat == 6) + diabetes + cad + bmi_under + bmi_overweight + bmi_obese + smoking_current + smoking_former + nhanes_nhw + female, data = stacked_data[dataset == dataset_name, ], family = quasibinomial())
  
  modelCAN_MGI <- predict(MGI_cancer_model, type = "response", newdata = stacked_data)
  denom <- ifelse(stacked_data[, cancer] == 1, modelCAN_MGI, 1 - modelCAN_MGI)
  num <- ifelse(stacked_data[, cancer] == 1, modelCAN_NHANES, 1 - modelCAN_NHANES)
  cancer_factor <- (num[stacked_data[, dataset] == dataset_name] / denom[stacked_data[, dataset] == dataset_name])
  cancer_factor <- chopr(cancer_factor)
  CANCER_NHANES_UNCORRECTED = cancer_factor * WEIGHT_NHANES_NOCAN
  CANCER_NHANES_UNCORRECTED  = stacked_data[dataset == dataset_name, .N] * CANCER_NHANES_UNCORRECTED / sum(CANCER_NHANES_UNCORRECTED, na.rm = TRUE)
  
  data.table(
    "id"                  = stacked_data[dataset == dataset_name, ][[id_var]],
    "no_cancer_ipw"       = WEIGHT_NHANES_NOCAN,
    "cancer_ipw"          = CANCER_NHANES_UNCORRECTED
  )
  
}

# calculate poststratification weights from MGI
# ADAPTED FROM: /net/junglebook/home/kundur/EHR/Processed Code/Weighted_using_lauren_code_bb.R
poststratification <- function(
    mgi_data,
    id_var             = "id",
    last_entry_age_var = "AgeLastEntry",
    cancer_var         = "cancer",
    chd_var            = "cad",
    smoke_var          = "smoking_current",
    diabetes_var       = "diabetes",
    female_var         = "female",
    chop               = FALSE,
    use_female         = FALSE
) {
  
  if (use_female == FALSE) {
    cli_alert_warning("Estimating poststratification weights without `female` variable. IPW estimation includes female by default.")
  }
  
  Nobs    <- nrow(mgi_data)
  age_vec <- as.numeric(mgi_data[[last_entry_age_var]])
  which_between <- function(vec, mat) {
    which(data.table::between(x = vec, lower = mat[["lower"]], upper = mat[["upper"]]))
  }
  
  # cancer prevalence by age (US) 
  # https://seer.cancer.gov/csr/1975_2016/results_merged/topic_prevalence.pdf
  cancer_prevalence <- age_grp_table(
    lower_ages   = c(0, 10, 20, 30, 40, 50, 60, 70, 80),
    num_vec      = c(0.0899, 0.2023, 0.3922, 0.8989, 2.1532, 4.9326, 10.4420, 18.3168, 21.5939) / 100,
    num_var_name = "prevalence"
  )
  
  cancer_func_pop <- stepfun(x = cancer_prevalence[["lower"]], y = c(0, cancer_prevalence[["prevalence"]]), right = FALSE)
  b               <- aggregate(as.formula(paste0(cancer_var, " ~ ", last_entry_age_var)), FUN = mean, data = mgi_data)
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
  d_b               <- aggregate(as.formula(paste0(diabetes_var, " ~ ", last_entry_age_var)), FUN = mean, data = mgi_data)
  diabetes_func_mgi <- stepfun(x = d_b[, 1], y = c(0, d_b[, 2]), right = FALSE)
  
  # chd prevalence by age (us)
  # https://www.cdc.gov/nchs/fastats/heart-disease.htm
  # youngest group will have no people, no number provided in documentation
  chd_prevalence <- age_grp_table(
    lower_ages   = c(0, 18, 45, 65, 75),
    num_vec      = c(0.01, 0.01, 0.060, 0.155, 0.239),
    num_var_name = "prevalence"
  )
  
  chd_func_pop <- stepfun(x = chd_prevalence[["lower"]], y = c(0, chd_prevalence[["prevalence"]]), right = FALSE)
  chd_b        <- aggregate(as.formula(paste0(chd_var, " ~ ", last_entry_age_var)), FUN = mean, data = mgi_data)
  chd_func_mgi <- stepfun(x = chd_b[, 1], y = c(0, chd_b[, 2]), right = FALSE)
  
  # smoking prevalence by age (us)
  # https://www.cdc.gov/nchs/fastats/heart-disease.htm
  # youngest group will have no people, no number provided in documentation
  smoke_prevalence <- age_grp_table(
    lower_ages   = c(0, 18, 25, 45, 65),
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
    upper_offset = 0,
    num_vec      = (age_counts_male[["counts"]] + age_counts_female[["counts"]]) / total_population,
    num_var_name = "prevalence"
  )
  
  age_vals      <- age_prevalence[, group]
  mgi_age_group <- age_vals[as.vector(unlist(apply(as.matrix(age_vec), 1, FUN = which_between, mat = age_prevalence)))]
  age_func_pop  <- stepfun(x = age_prevalence[["lower"]], y = c(0, age_prevalence[["prevalence"]]), right = FALSE)
  age_func_mgi  <- stepfun(x = age_prevalence[["lower"]], y = as.numeric(c(0, table(cut(age_vec, breaks = c(0, age_prevalence[["upper"]]), labels = age_prevalence[["group"]], right = FALSE), useNA = "ifany") / Nobs)), right = FALSE)
  age_func_mgi2 <- stepfun(x = age_prevalence[["lower"]], y = c(0, as.numeric(table(factor(mgi_age_group, levels = age_vals))) / Nobs), right = FALSE)
  
  # sex distribution by age (us)
  # same source as above for age
  if (use_female == TRUE) {
  female_prevalence <- age_grp_table(
    lower_ages   = low_ages,
    upper_offset = 0,
    num_vec      = age_counts_female[["counts"]] / (age_counts_male[["counts"]] + age_counts_female[["counts"]]),
    num_var_name = "prevalence"
  )
  female_func_pop <- stepfun(x = female_prevalence[["lower"]], y = c(0, female_prevalence[["prevalence"]]), right = FALSE)
  female_b        <- aggregate(as.formula(paste0(female_var, " ~ ", last_entry_age_var)), FUN = mean, data = mgi_data)
  female_func_mgi <- stepfun(x = female_b[, 1], y = c(0, female_b[, 2]), right = FALSE)
  }
  
  # without cancer
  population <- fifelse(mgi_data[[diabetes_var]] == 1, diabetes_func_pop(age_vec), 1 - diabetes_func_pop(age_vec))
  population <- population * fifelse(mgi_data[[chd_var]] == 1, chd_func_pop(age_vec), 1 - chd_func_pop(age_vec))
  population <- population * fifelse(mgi_data[[smoke_var]] == 1, smoke_func_pop(age_vec), 1 - smoke_func_pop(age_vec))
  if (use_female == TRUE) {
    population <- population * fifelse(mgi_data[[female_var]] == 1, female_func_pop(age_vec), 1 - female_func_pop(age_vec))
    }
  population <- population * age_func_pop(age_vec)
  
  mgi <- fifelse(mgi_data[[diabetes_var]] == 1, diabetes_func_mgi(age_vec), 1 - diabetes_func_mgi(age_vec))
  mgi <- mgi * fifelse(mgi_data[[chd_var]] == 1, chd_func_mgi(age_vec), 1 - chd_func_mgi(age_vec))
  mgi <- mgi * fifelse(mgi_data[[smoke_var]] == 1, smoke_func_mgi(age_vec), 1 - smoke_func_mgi(age_vec))
  if (use_female == TRUE) {
    mgi <- mgi * fifelse(mgi_data[[female_var]] == 1, female_func_mgi(age_vec), 1 - female_func_mgi(age_vec))
  }
  mgi <- mgi * age_func_mgi(age_vec)
  
  if (chop == TRUE) {
    post_no_cancer <- chopr(population / mgi)
  } else {
    post_no_cancer <- population / mgi
  }
  post_no_cancer <- (Nobs * post_no_cancer) / sum(post_no_cancer, na.rm = TRUE)
  
  # with cancer
  population <- population * fifelse(mgi_data[[cancer_var]] == 1, cancer_func_pop(age_vec), 1 - cancer_func_pop(age_vec))
  mgi        <- mgi * fifelse(mgi_data[[cancer_var]] == 1, cancer_func_mgi(age_vec), 1 - cancer_func_mgi(age_vec))
  
  if (chop == TRUE) {
    post_cancer <- chopr(population / mgi)
  } else {
    post_cancer <- population / mgi
  }
  post_cancer <- (Nobs * post_cancer) / sum(post_cancer, na.rm = TRUE)
  
  return(
    data.table(
      id              = mgi_data[[id_var]],
      no_cancer_postw = post_no_cancer,
      cancer_postw    = post_cancer
    )
  )
  
}