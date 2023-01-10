library(data.table)

prepare_nhanes_data <- function(
    nhanes_data,
    age_var        = "RIDAGEYR",
    mec_wt_var     = "WTMEC2YR",
    sex_var        = "RIAGENDR",
    cancer_var     = "MCQ220",
    chd_var        = "MCQ160C",
    smoke_ev_var   = "SMQ020",
    smoke_curr_var = "SMQ040",
    diabetes_var   = "DIQ010",
    race_eth_var   = "RIDRETH1"
) {
  
  merged <- nhanes_data
  
  merged <- merged[get(age_var) >= 12, ]
  merged <- merged[get(mec_wt_var) != 0, ]
  merged[, SAMP_NHANES := 1 / get(mec_wt_var)]
  merged[, NHANES_MEC_WT := length(get(mec_wt_var)) * get(mec_wt_var) /
           sum(get(mec_wt_var))]
  merged[, SEX := fcase(
    get(sex_var) == 1, "Male",
    get(sex_var) == 2, "Female"
  )]
  
  yes_no_vars <- c(cancer_var, chd_var, smoke_ev_var)
  char_vars <- c(smoke_curr_var, diabetes_var, race_eth_var)
  yes_no_recode <- function(x) {
    fcase(
      x == "1", "Yes",
      x == "2", "No",
      x %in% c("7", "9"), NA_character_
    )
  }
  merged[, (c(yes_no_vars, char_vars)) :=
           lapply(.SD, as.character), .SDcols = c(yes_no_vars, char_vars)]
  merged[, (yes_no_vars) := lapply(.SD, yes_no_recode), .SDcols = yes_no_vars]
  
  merged[, SMOKING_CURRENT := fcase(
    get(smoke_curr_var) %in% c("1", "2"), "Yes",
    get(smoke_curr_var) == "3", "No",
    get(smoke_curr_var) %in% c("7", "9"), NA_character_
  )]
  merged[, DIABETES := fcase(
    get(diabetes_var) == "1", "Yes",
    get(diabetes_var) %in% c("2", "3"), "No",
    get(diabetes_var) %in% c("7", "9"), NA_character_
  )]
  merged[BMXBMI > 86.2, BMXBMI := NA]
  merged[, RACE_ETH := fcase(
    get(race_eth_var) %in% c("1", "2"), "Hispanic",
    get(race_eth_var) == "3", "Non-Hispanic White",
    get(race_eth_var) == "4", "Non-Hispanic Black",
    get(race_eth_var) == "5", "Other"
  )]
  merged[, NHANES_NHW := as.numeric(RACE_ETH == "Non-Hispanic White")]
  
  merged[, Smoking := fcase(
    get(smoke_ev_var) == "No", "Never",
    get(smoke_ev_var) == "Yes" & SMOKING_CURRENT == "No", "Former",
    get(smoke_ev_var) == "Yes" & SMOKING_CURRENT == "Yes", "Current"
  )]
  
  merged[, NHANES_AGECAT := fcase(
    between(get(age_var), 0, 5), 1,
    between(get(age_var), 6, 11), 2,
    between(get(age_var), 12, 19), 3,
    between(get(age_var), 20, 39), 4,
    between(get(age_var), 40, 59), 5,
    between(get(age_var), 60, 150), 6
  )]
  
  merged[, NHANES_BMICAT := fcase(
    between(BMXBMI, 0, 18.499), 1,
    between(BMXBMI, 18.5, 24.999), 2,
    between(BMXBMI, 25.0, 29.999), 3,
    between(BMXBMI, 30, 120), 4
  )]
  
  tidy_table <- data.table(
    dataset          = rep("NHANES", nrow(merged)),
    age              = merged[, get(age_var)],
    female           = as.numeric(merged[, SEX] == "Female"),
    age_cat          = merged[, NHANES_AGECAT],
    nhw              = merged[, NHANES_NHW],
    cancer           = as.numeric(merged[, get(cancer_var)] == "Yes"),
    cancer_missing   = as.numeric(is.na(merged[, get(cancer_var)])),
    diabetes         = as.numeric(merged[, DIABETES] == "Yes"),
    diabetes_missing = as.numeric(is.na(merged[, DIABETES])),
    chd              = as.numeric(merged[, get(chd_var)] == "Yes"),
    chd_missing      = as.numeric(is.na(merged[, get(chd_var)])),
    bmi_under        = as.numeric(merged[, NHANES_BMICAT] == 1),
    bmi_overweight   = as.numeric(merged[, NHANES_BMICAT] == 3),
    bmi_obese        = as.numeric(merged[, NHANES_BMICAT] == 4),
    bmi_missing      = as.numeric(is.na(merged[, NHANES_BMICAT])),
    smoking_current  = as.numeric(merged[, Smoking] == "Current"),
    smoking_former   = as.numeric(merged[, Smoking] == "Former"),
    smoking_missing  = as.numeric(is.na(merged[, Smoking])),
    race_hispanic    = as.numeric(merged[, RACE_ETH] == "Hispanic"),
    race_black       = as.numeric(merged[, RACE_ETH] == "Non-Hispanic Black"),
    race_other       = as.numeric(merged[, RACE_ETH] == "Other"),
    samp_nhanes      = merged[, SAMP_NHANES],
    weight_nhanes    = merged[, NHANES_MEC_WT],
    bmi              = merged[, BMXBMI]
  )
  
  tidy_table[, weight_nhanes := length(weight_nhanes) * weight_nhanes /
               sum(weight_nhanes)][]
  
}
