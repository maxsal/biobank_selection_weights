# beta regression-based weights
# requires: mgi pim, demo, smk, and bmi files; nhanes data
# outputs:  nhanes and uncorrected cancer weights
# author:   max salvatore
# date:     20230110

# 1. libraries, functions, and options (outcome agnostic) ----------------------
options(stringsAsFactors = FALSE)

library(data.table)
# library(MatchIt)
library(cli)
library(optparse)
library(glue)
library(betareg)
library(haven)
library(survey)
library(ggplot2)
library(colorblindr)
library(simplexreg)

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
  make_option("--nhanes_wave_letter", type = "character",
              default = "J",
              help = glue("NHANES data prefix corresponding to wave ",
                          "[default = J]")),
  make_option("--nhanes_wave_years", type = "character",
              default = "2017-2018",
              help = glue("NHANES wave years corresponding to wave ",
                          "[default = 2017-2018]"))
)

parser <- OptionParser(usage = "%prog [options]", option_list = option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

# 2. specifications ------------------------------------------------------------
time_thresholds    <- as.numeric(opt$time_threshold)

## pull file paths corresponding to the data version specified
file_paths <- get_files(mgi_version = opt$mgi_version,
                        ukb_version = opt$ukb_version)

# 3. read data -----------------------------------------------------------------
## nhanes
cli_alert_info("loading nhanes data...")

nhanes_datasets  <- c("DEMO", "BMX", "SMQ", "DIQ", "MCQ")
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
  mgi_smk <- fread(file = file_paths[["mgi"]][["smk_file"]],
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
  mgi_bmi <- fread(file_paths[["mgi"]][["bmi_file"]],
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
  mgi_merged[][, dataset := "MGI"]
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

select_nhanes_nocan <- p_Sext * (p_MGI / (1 - p_MGI))
weights_nhanes_nocan <- 1/select_nhanes_nocan

weights_nhanes_nocan <- Nobs * weights_nhanes_nocan /
  sum(weights_nhanes_nocan, na.rm = TRUE)

## with cancer - as exposure?
select_nhanes_can <- betareg(samp_nhanes ~ as.numeric(age_cat %in% c(5, 6))
                             + chd + diabetes + bmi_under + bmi_overweight +
                               + bmi_obese + nhw + cancer + 
                               as.numeric(smoking_current + smoking_former),
                             data = stacked[dataset == "NHANES", ])
mgiselect_can <- glm(as.numeric(dataset == "MGI") ~
                         as.numeric(age_cat %in% c(5, 6)) + chd + diabetes +
                         bmi_under + bmi_overweight + bmi_obese + cancer +
                         smoking_current + smoking_former + nhw,
                       data = stacked, family = quasibinomial())

p_Sext_can <- predict(select_nhanes_can, newdata = stacked[dataset == "MGI", ],
                      type = "response")
p_MGI_can <- predict(mgiselect_can, newdata = stacked[dataset == "MGI", ],
                 type = "response")

select_nhanes_can <- p_Sext_can * (p_MGI_can / (1 - p_MGI_can))
weights_nhanes_can <- 1/select_nhanes_can

weights_nhanes_can <- Nobs * weights_nhanes_can /
  sum(weights_nhanes_can, na.rm = TRUE)

## with cancer - separate
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
  sum(cancer_nhanes_uncorrected, na.rm = TRUE)

# 6. output --------------------------------------------------------------------
weight_data <- data.table(
  id                     = stacked[dataset == "MGI", id],
  weights_no_cancer      = weights_nhanes_nocan,
  weights_cancer_direct   = weights_nhanes_can,
  weights_cancer_indirect = cancer_nhanes_uncorrected
)

# truncate weights at 10
weight_data[weights_no_cancer > 10]       <- 10
weight_data[weights_cancer_direct > 10]   <- 10
weight_data[weights_cancer_indirect > 10] <- 10

# save weights
fwrite(weight_data,
       glue("results/mgi/{opt$mgi_version}/cancer_weights.txt"))

# transform to long data for plotting
weight_data_long <- melt(weight_data, id.vars = c("id"))
weight_data_long[variable == "weights_no_cancer", name := "No cancer"]
weight_data_long[variable == "weights_cancer_direct", name := "Cancer - direct"]
weight_data_long[variable == "weights_cancer_indirect", name := "Cancer - indirect"]
weight_data_long[, name := factor(name, levels = c("No cancer", "Cancer - indirect", "Cancer - direct"))]

# 7. summary -------------------------------------------------------------------
mgi_sum <- merge.data.table(mgi_merged, weight_data)
mgi_sum <- mgi_sum[!is.na(weights_no_cancer), ]
table(is.na(mgi_sum[, weights_no_cancer]))
table(is.na(mgi_sum[, weights_cancer_indirect]))

## describe weight distribution
vw_plot <- weight_data_long |>
  ggplot(aes(name, value)) +
  geom_violin(aes(fill = name),
              draw_quantiles = c(0.25, 0.5, 0.75)) +
  labs(
    title   = "Estimated selection adjustment weights",
    x       = "",
    y       = "Estimated selection weights",
    caption = glue("All weights were trimmed at 10. Horizontal black lines correspond to the 25th, 50ths, and 75th quantiles.\n",
                   "'No cancer' weights were calculated using a Beta regression model adjusted for age, coronary heart diease,\ndiabetes, BMI, smoking, and non-Hispanic White.\n",
                   "'Cancer - indirect' are the 'no cancer' weights multiplied by a ratio of predicted cancer in MGI and NHANES.\n",
                   "'Cancer - direct' uses the same model as 'cnao cancer' but additionally adjusts for cancer.car")
  ) +
  scale_fill_OkabeIto() +
  theme_minimal() +
  theme(
    plot.title      = element_text(hjust = 0, face = "bold"),
    plot.caption    = element_text(hjust = 0),
    legend.position = "none"
  )
ggsave(
 filename = glue("results/mgi/{opt$mgi_version}/cancer_weights_violin_plot.pdf"),
 plot = vw_plot,
 height = 5, width = 7,
 device = cairo_pdf
)

## survey designs
### CHECK DESIGN
mgi_dsn <- svydesign(ids = ~1,
                     weights = ~weights_no_cancer,
                     data = mgi_sum)
mgi_can_dsn <- svydesign(ids = ~1,
                         weights = ~weights_cancer_direct,
                         data = mgi_sum)
mgi_ncan_cor_dsn <- svydesign(ids = ~1,
                         weights = ~weights_cancer_indirect,
                         data = mgi_sum)
### CHECK DESIGN
nhanes_design <- svydesign(data    = nhanes_merged,
                           id      = ~SDMVPSU,
                           strata  = ~SDMVSTRA,
                           weights = ~WTMEC2YR,
                           nest    = TRUE)

## compare means
age_sum <- describe_var_weights(
  mgi_data = mgi_sum,
  nhanes_data = nhanes_merged,
  mgi_var_name = c("age"),
  nhanes_var_name = "RIDAGEYR",
  designs = list(
    "mgi_dsn" = mgi_dsn,
    "mgi_can_dsn" = mgi_can_dsn,
    "mgi_ncan_cor_dsn" = mgi_ncan_cor_dsn,
    "nhanes_design" = nhanes_design)
)

female_sum <- describe_var_weights(
  mgi_data = mgi_sum,
  nhanes_data = nhanes_merged,
  mgi_var_name = "female",
  nhanes_var_name = "RIAGENDR",
  nhanes_var_val = "2",
  designs = list(
    "mgi_dsn" = mgi_dsn,
    "mgi_can_dsn" = mgi_can_dsn,
    "mgi_ncan_cor_dsn" = mgi_ncan_cor_dsn,
    "nhanes_design" = nhanes_design)
)

nhw_sum <- describe_var_weights(
  mgi_data = mgi_sum,
  nhanes_data = nhanes_merged,
  mgi_var_name = "nhw",
  nhanes_var_name = "RIDRETH1",
  nhanes_var_val = "3",
  designs = list(
    "mgi_dsn" = mgi_dsn,
    "mgi_can_dsn" = mgi_can_dsn,
    "mgi_ncan_cor_dsn" = mgi_ncan_cor_dsn,
    "nhanes_design" = nhanes_design)
)

bmi_sum <- describe_var_weights(
  mgi_data = mgi_sum,
  nhanes_data = nhanes_merged,
  mgi_var_name = "bmi",
  nhanes_var_name = "BMXBMI",
  designs = list(
    "mgi_dsn" = mgi_dsn,
    "mgi_can_dsn" = mgi_can_dsn,
    "mgi_ncan_cor_dsn" = mgi_ncan_cor_dsn,
    "nhanes_design" = nhanes_design)
)

cancer_sum <- describe_var_weights(
  mgi_data = mgi_sum,
  nhanes_data = nhanes_merged,
  mgi_var_name = "cancer",
  nhanes_var_name = "MCQ220",
  nhanes_var_val = "1",
  designs = list(
    "mgi_dsn" = mgi_dsn,
    "mgi_can_dsn" = mgi_can_dsn,
    "mgi_ncan_cor_dsn" = mgi_ncan_cor_dsn,
    "nhanes_design" = nhanes_design)
)

summary_results <- rbindlist(list(
  age_sum, female_sum, nhw_sum, bmi_sum, cancer_sum
))

fwrite(
  summary_results,
  glue("results/mgi/{opt$mgi_version}/cancer_weighted_summary_means.txt")
)
