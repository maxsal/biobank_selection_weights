# aim one analysis pipeline
library(glue)

mgi_version     <- "20220822"
ukb_version     <- "20221117"
outcome         <- "157"
time_thresholds <- c(0, 0.5, 1, 2, 3, 5)
use_geno_pcs    <- TRUE

# phase 0 scripts --------------------------------------------------------------
## prepare mgi data
system(glue("/usr/bin/time -v -o logs/00_prepare_mgi_data.txt Rscript r/",
            "00_prepare_mgi_data.R --mgi_version={mgi_version}"))

# phase 1 scripts --------------------------------------------------------------
## conduct a PCA analysis in MGI and UKB
system(glue("/usr/bin/time -v -o logs/01_pca_analysis.txt Rscript r/",
            "01_pca_analysis.R"))

## calculate partial correlations in MGI
system(glue("/usr/bin/time -v -o logs/01_mgi_partial_correlations.txt Rscript r/",
            "01_mgi_partial_correlations.R --mgi_version={mgi_version} ",
            "--use_geno={use_geno_pcs}"))

## prepare time-restricted phenomes in MGI and UKB
system(glue("/usr/bin/time -v -o logs/01_prepare_phenomes.txt Rscript r/",
            "01_prepare_phenomes.R --mgi_version={mgi_version} --ukb_version={ukb_version} ",
            "--time_thresholds={paste0(time_thresholds, collapse = ',')} --outcome={outcome}"))

## estimate ipw and poststratification weights in MGI
system(glue("/usr/bin/time -v -o logs/01_estimate_weights.txt Rscript r/",
            "01_estimate_weights.R --cohort_version={mgi_version}"))

## phase 2 scripts -------------------------------------------------------------
for (i in seq_along(time_thresholds)) {
  cli::cli_alert_info("running random forest at {time_thresholds[i]}")
  system(glue("/usr/bin/time -v -o logs/02_random_forest.txt Rscript r/",
              "02_random_forest.R --mgi_version={mgi_version} --ukb_version={ukb_version} ",
              "--time_threshold={time_thresholds[i]} --outcome={outcome}"))
}

for (i in seq_along(time_thresholds)) {
  cli::cli_alert_info("running SuperLearner at {time_thresholds[i]}")
  system(glue("/usr/bin/time -v -o logs/02_super_learner.txt Rscript r/",
              "02_super_learner.R --mgi_version={mgi_version} --ukb_version={ukb_version} ",
              "--time_threshold={time_thresholds[i]} --outcome={outcome}"))
}

## phase 3 scripts -------------------------------------------------------------
for (i in time_thresholds) {
  system(glue("/usr/bin/time -v -o logs/03_calculate_naive_phers.txt Rscript r/",
              "03_calculate_naive_phers.R --outcome={outcome} --time_threshold={i}"))
}

## phase 4 scripts -------------------------------------------------------------
for (i in time_threshold) {
  system(glue("/usr/bin/time -v -o logs/04_evaluate_phers_new.txt Rscript r/",
              "04_evaluate_phers_new.R --outcome={outcome} --time_threshold={i}"))
}


