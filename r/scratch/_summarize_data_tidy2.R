compare_ors <- function(
    data,
    outcome,
    exposure,
    weights = c("ip_simple", "ip_simple_f", "ip_selection", "ip_selection_f", 
                "ip_cancer", "ip_depression", "ip_cad", "ip_diabetes", "ip_hypertension",
                "ps_selection", "ps_selection_f", "ps_nhw", "ps_nhw_f", "ps_depression", "ps_diabetes", "ps_hypertension"
    )) {
  out <- list()
  pb <- txtProgressBar(max = length(weights) + 1, width = 50, style = 3)
  mod <- suppressWarnings(suppressMessages(glm(as.formula(paste0(outcome, " ~ ", exposure)), data = data, family = "binomial")))
  est <- summary(mod)$coef[exposure, 1:2]
  out[[1]] <- data.table(
    outcome  = outcome,
    exposure = exposure,
    weights  = "None",
    beta_est = est[[1]],
    beta_lo  = est[[1]] - qnorm(0.975) * est[[2]],
    beta_hi  = est[[1]] + qnorm(0.975) * est[[2]],
    or_est   = exp(est[[1]]),
    or_lo    = exp(est[[1]] - qnorm(0.975) * est[[2]]),
    or_hi    = exp(est[[1]] - qnorm(0.975) * est[[2]])
  )
  setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
  for (w in seq_along(weights)) {
    mod      <- suppressWarnings(suppressMessages(glm(as.formula(paste0(outcome, " ~ ", exposure)), data = data, family = "binomial", weights = merged[[weights[w]]])))
    est <- summary(mod)$coef[exposure, 1:2]
    out[[w + 1]] <- data.table(
      outcome  = outcome,
      exposure = exposure,
      weights  = weights[w],
      beta_est = est[[1]],
      beta_lo  = est[[1]] - qnorm(0.975) * est[[2]],
      beta_hi  = est[[1]] + qnorm(0.975) * est[[2]],
      or_est   = exp(est[[1]]),
      or_lo    = exp(est[[1]] - qnorm(0.975) * est[[2]]),
      or_hi    = exp(est[[1]] - qnorm(0.975) * est[[2]])
    )
    setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
  }
  close(pb)
  do.call(rbind, out)
}

cancer_female       <- compare_ors(data = merged, outcome = "cancer", exposure = "female")
diabetes_female     <- compare_ors(data = merged, outcome = "diabetes", exposure = "female")
cad_female          <- compare_ors(data = merged, outcome = "cad", exposure = "female")
depression_female   <- compare_ors(data = merged, outcome = "depression", exposure = "female")

stacked <- rbindlist(list(
  cancer_female, diabetes_female, cad_female, depression_female
))

fwrite(stacked, "bin/weighted_ors.txt", sep = "\t")

#####
lqsum(
  mgi_cov,
  vars = c("age_at_last_diagnosis", "age_verbose", "sex", "race_eth", "bmi", "bmi_verbose",
    "cancer", "diabetes", "cad", "anxiety", "depression", "SmokingStatus")
)
