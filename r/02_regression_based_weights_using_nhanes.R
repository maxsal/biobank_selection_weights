# beta regression-based weights
# requires: mgi pim, demo, smk, and bmi files; nhanes data
# outputs:  nhanes and uncorrected cancer weights
# author:   max salvatore
# date:     20221220

# 1. libraries, functions, and options (outcome agnostic) ----------------------
options(stringsAsFactors = FALSE)

library(data.table)
library(MatchIt)
library(cli)
library(optparse)
library(glue)
library(betareg)

set.seed(61787)

for (i in list.files("fn/")) source(paste0("fn/", i)) # load functions
cancer_phecodes <- fread("data/public/cancer_phecodes.txt",
                         colClasses = "character")[[1]]

# optparse list ----------------------------------------------------------------
option_list <- list(
  make_option("--outcome", type = "character", default = "",
              help = "Outcome phecode"),
  make_option("--mgi_version", type = "character", default = "20210318",
              help = "Version of MGI data [default = 20210318]"),
  make_option("--ukb_version", type = "character", default = "20221117",
              help = "Version of UKB data [default = 20221117]"),
  make_option("--time_threshold", type = "character", default = "0",
              help = glue("Time threshold for the phenome data ",
                          "[default = 0]")),
  make_option("--nhanes_data_path", type = "character",
              default = "data/public/nhanes/",
              help = glue("Relative path to directory with NHANES data ",
                          "[default = data/public/nhanes/]")),
  make_option("--nhanes_data_prefix", type = "character",
              default = "P",
              help = glue("NHANES data prefix corresponding to wave ",
                          "[default = P]"))
)

parser <- OptionParser(usage = "%prog [options]", option_list = option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

# 2. specifications ------------------------------------------------------------
nhanes_data_path   <- opt$nhanes_data_path
nhanes_data_prefix <- op$nhanes_data_prefix
time_thresholds    <- as.numeric(opt$time_threshold)

## pull file paths corresponding to the data version specified
file_paths <- get_files(mgi_version = opt$mgi_version,
                        ukb_version = opt$ukb_version)

# 3. read data -----------------------------------------------------------------
## nhanes
cli_alert_info("loading nhanes data...")

nhanes_datasets  <- c("DEMO", "BMX", "SMQ", "DIQ", "MCQ")
nhanes_data_list <- list()
for (i in seq_along(nhanes_datasets)) {
  nhanes_data_list[[i]] <- read_xpt(glue("{nhanes_data_path}",
                                         "{nhanes_data_prefix}_",
                                         "{nhanes_datasets[i]}.XPT")) |>
    as.data.table()
}

nhanes_merged <- Reduce(\(x, y) {
  merge.data.table(x = x, y = y, by = "SEQN", all = TRUE)
},
data_list
)[, c("SEQN", "RIAGENDR", "WTINTPRP", "RIDAGEYR", "RIDRETH1", "MCQ220",
      "BMXBMI", "SMQ040", "SMQ020", "DIQ010", "MCQ160C", "WTMECPRP")]

prepped_nhanes <- prepare_nhanes_data(nhanes_data = nhanes_merged)

## mgi
cli_alert_info("loading mgi data...")
### phecode indicator matrix (PEDMASTER_0)
  mgi_pim0 <- fread(file_paths[["mgi"]]$pim0_file)
  setnames(mgi_pim0,
           old = c("IID"),
           new = c("id"))
  
  Nobs <- length(unlist(mgi_pim0[, 1]))
  
  mgi_pim0[, cancer := rowSums(.SD, na.rm = TRUE),
           .SDcols = glue("X{cancer_phecodes}")]
  mgi_pim0[, cancer := fifelse(cancer >= 1, 1, 0)]

### demographics data
  mgi_demo <- fread(file_paths[["mgi"]]$demo_file)[, .(
    id       = Deid_ID,
    age      = Age,
    ethn     = EthnicityName,
    sex      = Sex,
    race     = RaceName)]
  mgi_demo[, `:=`(
    female = as.numeric(sex == "F"),
    age_cat = fcase(
      between(age, 0, 5), 1,
      between(age, 6, 11), 2,
      between(age, 12, 19), 3,
      between(age, 20, 39), 4,
      between(age, 40, 59), 5,
      between(age, 60, 150), 6
    ),
    nhw = as.numeric(ethn == "Non-Hispanic or Latino" & race == "Caucasian"),
    race_black = as.numeric(ethn == "Non-Hispanic or Latino" &
                              race == "African American"),
    race_hispanic = as.numeric(ethn == "Hispanic or Latino"),
    race_other = as.numeric(ethn != "Hispanic or Latino" &
                              !(race %in% c("African American", "Caucasian")))
  )]

### smoking data
  mgi_smk <- fread(file_paths[["mgi"]]["smk_file"],
                   select = c("DeID_PatientID", "DaysSinceBirth",
                              "SmokingStatus"))
  setnames(mgi_smk, c("DeID_PatientID", "DaysSinceBirth", "SmokingStatus"),
           c("id", "dsb", "smoking"))
  mgi_smk <- mgi_smk[ mgi_smk[, .I[which.max(dsb)], by = "id"][, V1] ][
    , !c("dsb")]
  
  mgi_smk[, `:=`(
    smoking_current = as.numeric(smoking == "Current"),
    smoking_former = as.numeric(smoking == "Former")
  )]

### bmi data
  mgi_bmi <- fread(file_paths[["mgi"]]["bmi_file"],
                   select = c("DeID_PatientID", "DaysSinceBirth", "BMI"))
  setnames(mgi_bmi, c("DeID_PatientID", "DaysSinceBirth", "BMI"),
           c("id", "dsb", "bmi"))
  mgi_bmi <- mgi_bmi[, !c("dsb")][, .(bmi = median(as.numeric(bmi),
                                                   na.rm = TRUE)), by = "id"]
  mgi_bmi[, `:=`(
    bmi_cat = fcase(
      between(bmi, 0, 18.499), 1,        # underweight
      between(bmi, 18.5, 24.999), 2,     # "normal"
      between(bmi, 25.0, 29.999), 3,     # overweight
      between(bmi, 30, 120), 4)          # obese
  )][, `:=`(
    bmi_under = as.numeric(bmi_cat == 1),
    bmi_overweight = as.numeric(bmi_cat == 3),
    bmi_obese = as.numeric(bmi_cat == 4)
  )]
  
# 4. prepare and stack data ----------------------------------------------------
cli_alert_info("prepping and stacking data...")
mgi_dis  <- mgi_pim0[, .(id, chd = X411.4, diabetes = X250, cancer)]
mgi_covs <- mgi_demo[, .(id, age, female, age_cat, nhw, race_black,
                         race_hispanic, race_other)]

mgi_merged <- Reduce(\(x, y) {
  merge.data.table(x = x, y = y, by = "id", all = TRUE)
},
list(mgi_covs, mgi_dis, mgi_smk, mgi_bmi)
)[, !c("smoking")]

stacked <- rbindlist(list(
  prepped_nhanes,
  mgi_merged[, !c("id")][, dataset := "MGI"]
), use.names = TRUE, fill = TRUE)

# 5. calculate weights ---------------------------------------------------------
cli_alert_info("calculating weights...")
## without cancer
select_nhanes_nocan <- betareg(samp_nhanes ~ as.numeric(age_cat %in% c(5, 6))
                               + chd + diabetes + bmi_under + bmi_overweight +
                                 + bmi_obese + nhw +
                                 as.numeric(smoking_current + smoking_former),
                               data = stacked[dataset == "NHANES", ])

mgiselect_nocan <- glm(as.numeric(dataset == "MGI") ~
                         as.numeric(age_cat %in% c(5, 6)) + chd + diabetes +
                         bmi_under + bmi_overweight + bmi_obese +
                         smoking_current + smoking_former + nhw,
                       data = stacked, family = quasibinomial())

p_Sext <- predict(select_nhanes_nocan, newdata = stacked[dataset == "MGI", ],
                  type = "response")
p_MGI <- predict(mgiselect_nocan, newdata = stacked[dataset == "MGI", ],
                 type = "response")

# tmp <- rep(0, times = length(p_MGI))
# tmp[which(rownames(p_Sext) %in% rownames(data.frame(p_MGI)) == TRUE)] <- p_Sext
# tmp[which(rownames(p_Sext) %in% rownames(data.frame(p_MGI)) == TRUE)] <- NA
# p_Sext <- tmp
# p_Sext[which(p_Sext == 0)] <- 1.921e-05

select_nhanes_nocan <- p_Sext * (p_MGI / (1 - p_MGI))
weights_nhanes_nocan <- 1/select_nhanes_nocan
weights_nhanes_nocan <- Nobs * weights_nhanes_nocan /
  sum(weights_nhanes_nocan, na.rm = TRUE)

## with cancer
nhanes_cancer_model <- glm(cancer ~ as.numeric(age_cat %in% c(5, 6)) + diabetes
                           + chd + bmi_under + bmi_overweight + bmi_obese +
                             smoking_current + smoking_former + nhw,
                           data = stacked[dataset == "NHANES", ],
                           weights = stacked[dataset == "NHANES", weight_nhanes],
                           family = quasibinomial())
mcan_nhanes <- predict(nhanes_cancer_model, type = "response", newdata = stacked)

mgi_cancer_model <- glm(cancer ~ as.numeric(age_cat %in% c(5, 6)) + diabetes
                        + chd + bmi_under + bmi_overweight + bmi_obese +
                          smoking_current + smoking_former + nhw,
                        data = stacked[dataset == "MGI", ],
                        family = quasibinomial())
mcan_mgi <- predict(mgi_cancer_model, type = "response", newdata = stacked)

denom <- ifelse(stacked[, cancer] == 1, mcan_mgi, 1 - mcan_mgi)
num   <- ifelse(stacked[, cancer] == 1, mcan_nhanes, 1 - mcan_nhanes)

cancer_nhanes_uncorrected <- ((num[stacked[, dataset] == "MGI"]) /
                                (denom[stacked[, dataset] == "MGI"])) *
  weights_nhanes_nocan
cancer_nhanes_uncorrected <- Nobs * cancer_nhanes_uncorrected /
  sim(cancer_nhanes_uncorrected, na.rm = TRUE)

# 6. output --------------------------------------------------------------------
data.table(
  id                        = stacked[dataset == "MGI", id],
  weight_nhanes_nocan       = weights_nhanes_nocan,
  cancer_nhanes_uncorrected = cancer_nhanes_uncorrected
)