# Calculate simplex regression-based inverse probability weights using NHANES
# data and postratificatiod-based weights using Census, CDC, and SEER data in
# MGI
# author:   max salvatore
# date:     20230809

# libraries --------------------------------------------------------------------
ms::libri(
  maxsal/ms, haven, survey, dplyr, pracma, simplexreg, data.table,
  glue, cli, qs, optparse, ggExtra
)

# optparse list ---
option_list <- list(
  make_option("--cohort_version",
    type = "character", default = "20220822",
    help = "Cohort version in /net/junglebook/magic_data/EHRdata/ [default = '20220822']"
  ),
  make_option("--mgi_cohort",
    type = "character", default = "comb",
    help = "Specific MGI cohort [default = %default]"
  ),
  make_option("--nhanes_wave_letter",
    type = "character",
    default = "J",
    help = glue(
      "NHANES data prefix corresponding to wave ",
      "[default = %default]"
    )
  ),
  make_option("--nhanes_wave_years",
    type = "character",
    default = "2017-2018",
    help = glue(
      "NHANES wave years corresponding to wave ",
      "[default = %default]"
    )
  ),
  make_option("--nhanes_survey_names",
    type = "character",
    default = "DEMO,BMX,SMQ,DIQ,MCQ,DPQ,BPQ,HIV",
    help = glue(
      "NHANES wave years corresponding to wave ",
      "[default = %default]"
    )
  )
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args   <- parse_args(parser, positional_arguments = 0)
opt    <- args$options
print(opt)

data_path <- glue("data/private/mgi/{opt$cohort_version}/")

for (i in list.files("fn", full.names = TRUE)) source(i)

# load data --------------------------------------------------------------------
cli_alert("loading data...")
# MGI
mgi <- read_qs(glue("{data_path}datax_{opt$cohort_version}_{opt$mgi_cohort}.qs"))[, `:=` (
  smoking_current = as.numeric(SmokingStatus == "Current"),
  smoking_former  = as.numeric(SmokingStatus == "Former")
  )]
setnames(mgi, c("DeID_PatientID", "age_at_last_diagnosisx"), c("id", "age"), skip_absent = TRUE)

mgi[, bmi_cat := relevel(factor(fcase(
  bmi_cat == 1, "Underweight",
  bmi_cat == 2, "Normal weight",
  bmi_cat == 3, "Overweight",
  bmi_cat == 4, "Obese",
  default = "Unknown"
)), ref = "Normal weight")]
mgi[, c("cancer", "diabetes", "cad") := NULL]
setnames(mgi, c("cancerx", "cadx", "diabetesx"), c("cancer", "cad", "diabetes"))

mgi[, race_eth := fcase(
  RaceEthnicity == "Caucasian / Non-Hispanic", "NH White",
  RaceEthnicity == "African American / Non-Hispanic", "NH Black",
  RaceEthnicity %in% c(
      "Asian / Hispanic", "Caucasian / Hispanic", "Native American / Hispanic", "African American / Hispanic"
  ), "Hispanic",
  RaceEthnicity == "Asian / Hispanic", "NH Asian",
  default = "Other/Unknown"
)]

# NHIS
prepped_nhis <- fread("https://raw.githubusercontent.com/maxsal/public_data/main/prepped_nhis.csv")[, dataset := "NHIS"]
setnames(prepped_nhis, "chd", "cad")
prepped_nhis[, samp_WTFA_A := 1/WTFA_A]
prepped_nhis[, bmi_cat := relevel(factor(bmi_cat), ref = "Normal weight")]

# stack
stacked <- rbindlist(list(
  prepped_nhis,
  mgi[][, dataset := "MGI"]
), use.names = TRUE, fill = TRUE)

# estimate ipw and postratification weights ------------------------------------
cli_alert("estimating ipw weights...")

ip_weights_list <- list(
  "selection" = c("as.numeric(age>= 50)", "female", "nhw", "hypertension", "diabetes", "cancer", "anxiety", "depression", "bmi_cat")
)

ip_weights <- list()
cli_progress_bar(name = "estimating ipw weights", total = length(names(ip_weights_list)))
for (i in seq_along(names(ip_weights_list))) {
  tmp_weights <- ipw(
    stacked_data       = stacked,
    weight_outcome_var = "WTFA_A",
    samp_var           = "samp_WTFA_A",
    external_dataset   = "NHIS",
    covs               = ip_weights_list[[names(ip_weights_list)[i]]]
  )
  setnames(tmp_weights, "ip_weight", paste0("ip_", names(ip_weights_list)[i]))
  ip_weights[[i]] <- tmp_weights
  cli_progress_update()
}
cli_progress_done()

ipws <- ms::merge_list(ip_weights, join_fn = dplyr::full_join)

cli_alert("estimating poststratification weights...")
post_weights_list <- list(
  "selection" = c("female", "nhw", "hypertension", "diabetes", "cancer", "anxiety", "depression", "bmi_cat")
)

post_weights <- list()
cli_progress_bar(name = "estimating poststratification weights", total = length(names(post_weights_list)))
for (i in seq_along(names(post_weights_list))) {
  tmp_weights <- psw(
    int_data       = mgi,
    ext_data       = prepped_nhis,
    psu_var        = "PPSU",
    strata_var     = "PSTRAT",
    weight_var     = "WTFA_A",
    chop           = TRUE,
    age_bin        = TRUE,
    age_bin_breaks = c(0, 50, 150),
    not_num_vars   = c("bmi_cat", "race_eth"),
    covs           = post_weights_list[[names(post_weights_list)[i]]]
  )
  setnames(tmp_weights, "ps_weight", paste0("ps_", names(post_weights_list)[i]))
  post_weights[[i]] <- tmp_weights
  cli_progress_update()
}
cli_progress_done()

posts <- merge_list(post_weights, join_fn = dplyr::full_join)

weights <- merge.data.table(ipws, posts, by = "id")
merged  <- merge_list(list(mgi, weights), by = "id", join_fn = dplyr::full_join)

# estimate cancer~female log(OR) -----------------------------------------------
cli_alert("estimating cancer~female log(OR) as sanity check...")

extract_estimates <- function(
  outcome,
  exposure,
  data,
  weights = NULL
) {
  if (is.null(weights)) {
    mod <- glm(
      as.formula(paste0(outcome, " ~ ", exposure)),
      family  = "quasibinomial",
      data    = data)
  } else {
    tmp_design <- svydesign(
      id      = ~1,
      weights = ~get(weights),
      data    = data[!is.na(get(weights)), ]
    )
    mod <- svyglm(
      as.formula(paste0(outcome, " ~ ", exposure)),
      family = "quasibinomial",
      design = tmp_design
    )
  }

  est     <- summary(mod)$coef[exposure, 1:2]
  beta_lo <- est[[1]] - qnorm(0.975) * est[[2]]
  beta_hi <- est[[1]] + qnorm(0.975) * est[[2]]
  rare <- ifelse(sum(mod$data[[outcome]], na.rm = TRUE) /
    length(na.omit(mod$data[[outcome]])) >= 0.15, FALSE, TRUE)
  evals <- suppressMessages(evalues.OR(
    est  = exp(est[[1]]),
    lo   = exp(beta_lo),
    hi   = exp(beta_hi),
    rare = rare
  ))
  data.table(
    outcome  = outcome,
    exposure = exposure,
    weights  = ifelse(is.null(weights), "None", weights),
    beta     = est[[1]],
    beta_lo  = beta_lo,
    beta_hi  = beta_hi,
    or       = exp(est[[1]]),
    or_lo    = exp(beta_lo),
    or_hi    = exp(beta_hi),
    eval     = evals[2, 1],
    eval_lo  = evals[2, 2],
    eval_hi  = evals[2, 3]
  )
}

est_out <- list()
weight_vars <- names(weights)[names(weights) != "id"]
cli_progress_bar(name = "estimating cancer~female log(OR)", total = length(weight_vars))
for (i in seq_along(weight_vars)) {
  est_out[[i]] <- extract_estimates(
    outcome  = "cancer",
    exposure = 'as.numeric(Sex == "F")',
    data     = merged,
    weights  = weight_vars[i]
  )
  cli_progress_update()
}
cli_progress_done()

est_out[[i+1]] <- extract_estimates(
        outcome  = "cancer",
        exposure = 'as.numeric(Sex == "F")',
        data     = merged
      )

est_out <- rbindlist(est_out)

fwrite(
  x    = est_out,
  file = glue("{data_path}cancerx_female_logor_est_{opt$cohort_version}_{opt$mgi_cohort}.txt"),
  sep  = "\t"
)

# save -------------------------------------------------------------------------
save_qs(
  x    = weights,
  file = glue("{data_path}weightsx_{opt$cohort_version}_{opt$mgi_cohort}.qs")
)

scat_plot <- weights |>
  ggplot(aes(x = ip_selection, y = ps_selection)) +
  geom_point(size = 1, alpha = 0.5) +
  # geom_smooth(formula = "y ~ x", method = "loess", se = FALSE) +
  labs(title = str_wrap("MGI weight scatterplot with marginal distibution", width = 50), x = "IP weight", y = "PS weight") +
  theme_ms()

marg_plot <- ggMarginal(scat_plot, type = "violin", margins = "both")

ggsave(
  plot     = marg_plot,
  filename = glue("results/mgi/{opt$cohort_version}/mgi_weights_dist.png"),
  width    = 8,
  height   = 8,
  units    = "in",
  dpi      = 300
)

mgi_weight_violin_plot <- weights[, .(`IP weight` = ip_selection, `PS weight` = ps_selection)] |>
  melt() |>
  ggplot(aes(x = variable, y = value, fill = variable)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  labs(x = "Weight type", y = "Weight value", title = "MGI weight distribution") +
  scale_fill_ms() +
  theme_ms() +
  theme(
    legend.position = "none"
  )

ggsave(
  plot = mgi_weight_violin_plot,
  filename = glue("results/mgi/{opt$cohort_version}/mgi_weights_violin.pdf"),
  width = 6, height = 6, device = cairo_pdf
)

cli_alert_success("script success! see {.path {data_path}} and suffix {opt$mgi_cohort} 🥳")
