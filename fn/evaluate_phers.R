# evaluate the performance of a phers in a dataset with outcome (case),
# predictor (phers), and covariates
# author: max salvatore
# date:   20230118

evaluate_phers <- function(
    merged_data,
    covars,
    outcome,
    phers_name,
    n_phecodes,
    calc_r2 = TRUE,
    pctile_or = TRUE
) {
  # initialize
  tmp_max <- NULL
  tmp_lars <- NULL
  
  tmp_max$trait     <- gsub("X", "", outcome)
  tmp_max$predictor <- phers_name
  tmp_max$phecodes  <- n_phecodes
  tmp_lars$trait     <- gsub("X", "", outcome)
  tmp_lars$predictor <- phers_name
  tmp_lars$phecodes  <- tmp_max$phecodes
  
  cov_check <- apply(merged_data[, ..covars], 2, \(x) length(unique(x)))
  x_cov <- which(cov_check > 1)
  add_cov <- ifelse(length(x_cov) > 0,
                    paste0("+", names(cov_check[x_cov]), collapse = "", sep = ""),
                    "")
  
  cli_alert_info("getting 'max' diagnostics....")
  max_m1 <- glm(glue("case ~ {phers_name} {add_cov}"), data = merged_data, family = binomial())
  vars <- diag(vcov(max_m1))
  
  # error_catcher <- try(-pchisq((max_m1$coefficient[[phers_name]]^2/vars[[phers_name]]), 1, lower.tail = FALSE, log = TRUE) / log(10))
  # if(inherits(t, "try-error")) {
  #   stop("Problem calculating tmp_max$log10p")
  # }
  max_log10p <- -pchisq((max_m1$coefficient[[phers_name]]^2/vars[[phers_name]]), 1,
                        lower.tail = FALSE, log = TRUE)/log(10)
  
  if (coef(summary(max_m1))[phers_name, 4] == 0) {
    tmp_max$p_val <- log10toP(max_log10p)
  } else {
    tmp_max$p_val <- coef(summary(max_m1))[phers_name, 4]
  }
  tmp_max$beta   <- max_m1$coefficient[[phers_name]]
  tmp_max$sebeta <- sqrt(diag(vcov(max_m1))[[phers_name]])
  tmp_max$or     <- exp(max_m1$coefficient[[phers_name]])
  tmp_confint    <- suppressMessages(confint(max_m1))
  tmp_max$or_lo  <- exp(tmp_confint[phers_name, 1])
  tmp_max$or_hi  <- exp(tmp_confint[phers_name, 2])
  tmp_max$or_print  <- glue("{format(round(tmp_max$or, 2), nsmall = 2)} (",
                            "{format(round(tmp_max$or_lo, 2), nsmall = 2)}, ",
                            "{format(round(tmp_max$or_hi, 2), nsmall = 2)})")
  tmp_max$log10p <- max_log10p
  
  max_auc <- suppressMessages(
    roc(response = merged_data[["case"]], predictor = merged_data[[phers_name]],
        family = binomial(), ci = TRUE)
  )
  tmp_max$auc <- max_auc$auc
  tmp_max$auc_lo <- max_auc$ci[1]
  tmp_max$auc_hi <- max_auc$ci[3]
  tmp_max$auc_print <- glue(
    "{format(round(tmp_max$auc, 3), nsmall = 3)} (",
    "{format(round(tmp_max$auc_lo, 3), nsmall = 3)}, ",
    "{format(round(tmp_max$auc_hi, 3), nsmall = 3)})"
  )
  
  max_m2 <- glm(glue("case ~ {phers_name}"), data = merged_data, family = binomial())
  if (calc_r2 == TRUE) {
    TOGGLE = (class(max_m2)[1] == "lm" | class(max_m2)[1] == "glm")
    if (!TOGGLE) {stop("Function written for `lm` and `glm` objects only")}
    null_m2 <- update(max_m2, ~ 1)
    N_m2    <- nobs(max_m2)
    m_m2    <- suppressWarnings(logLik(max_m2, REML = FALSE))[1]
    n_m2    <- suppressWarnings(logLik(null_m2, REML = FALSE))[1]
    
    tmp_max$r2_mcfadden      <- signif(1 - m_m2/n_m2, digits = 6)
    cs_m2                <- 1 - exp(-2/N_m2 * (m_m2 - n_m2))
    tmp_max$r2_cox_and_snell <- signif(cs_m2, digits = 6)
    tmp_max$r2_nagelkerke    <- signif(cs_m2 / (1 - exp((2 / N_m2) * n_m2)),
                                       digits = 6)
  }
  
  max_m2_hl           <- hoslem.test(max_m2$y, fitted(max_m2), g = 10)
  tmp_max$hl_chisq    <- max_m2_hl$statistic[["X-squared"]]
  tmp_max$hl_p        <- max_m2_hl$p.value
  tmp_max$brier_score <- BrierScore(max_m2)
  
  cli_alert_info("performing cross validation....")
  # cross validation with 5 random folds ---
  G                <- 5
  merged_data$folds   <- createFolds(y = 1:nrow(merged_data), k = G, list = F)
  linearpredictors <- list()
  cli_progress_bar(name = "cross validation")
  for (g in 1:G) {
    testdata  <- merged_data[merged_data$folds == g, ]
    traindata <- merged_data[merged_data$folds != g, ]
    
    # fit
    formula.fit <- glue("case ~ {phers_name}")
    fit.train   <- logistf(as.formula(formula.fit), data = traindata,
                           control = logistf.control(maxit = 1000))
    
    # prediction
    betas <- coef(fit.train)
    X     <- model.matrix(as.formula(formula.fit), data = testdata)
    
    testdata$fittedPredictions <- as.numeric(1 / (1 + exp(-X %*% betas)))
    form_predict <- glue("logistf(case ~ fittedPredictions, ",
                         "family = binomial(logit), data = testdata, ",
                         "ci = TRUE, ",
                         "control = logistf.control(maxit = 1000))")
    pred_test               <- eval(parse(text = form_predict))
    linearpredictors[[g]] <- data.table("id" = testdata$id,
                                        "linearpredictor" = 
                                          pred_test$linear.predictors)
    cli_progress_update()
  }
  
  merged_data <- merge.data.table(
    merged_data, rbindlist(linearpredictors),
    by = "id"
  )
  
  cli_alert_success("cross validataion complete. getting 'lars' diagnostics....")
  
  lars_m1            <- logistf(glue("case ~ {phers_name} {add_cov}"), data = merged_data)
  lars_sebeta        <- sqrt(diag(vcov(lars_m1)))
  names(lars_sebeta) <- lars_m1$terms
  if (glue("{phers_name}{1}") %in% lars_m1$terms) {
    predictor <- glue("{phers_name}{1}")
  }
  vars        <- diag(lars_m1$var)
  names(vars) <- lars_m1$terms
  lars_log10p <- -pchisq((lars_m1$coefficient[[phers_name]]^2/vars[[phers_name]]), 1,
                         lower.tail = FALSE, log = TRUE)/log(10)
  
  tmp_lars$p_val  <- log10toP(lars_log10p)
  tmp_lars$beta   <- as.numeric(lars_m1$coefficient[[phers_name]])
  tmp_lars$sebeta <- as.numeric(lars_sebeta[[phers_name]])
  tmp_lars$or       <- exp(lars_m1$coefficient[[phers_name]])
  tmp_confint       <- suppressMessages(confint(lars_m1))
  tmp_lars$or_lo    <- exp(tmp_confint[phers_name, 1])
  tmp_lars$or_hi    <- exp(tmp_confint[phers_name, 2])
  tmp_lars$or_print <- glue("{format(round(tmp_lars$or, 2), nsmall = 2)} (",
                            "{format(round(tmp_lars$or_lo, 2), nsmall = 2)}, ",
                            "{format(round(tmp_lars$or_hi, 2), nsmall = 2)})")
  tmp_lars$log10p <- lars_log10p
  
  lars_auc <- roc(
    response  = merged_data[["case"]],
    predictor = merged_data[["linearpredictor"]],
    family    = binomial(),
    ci        = TRUE
  )
  
  tmp_lars$auc    <- lars_auc$auc
  tmp_lars$auc_lo <- lars_auc$ci[1]
  tmp_lars$auc_hi <- lars_auc$ci[3]
  tmp_lars$auc_print <- glue(
    "{format(round(tmp_lars$auc, 3), nsmall = 3)} (",
    "{format(round(tmp_lars$auc_lo, 3), nsmall = 3)}, ",
    "{format(round(tmp_lars$auc_hi, 3), nsmall = 3)})"
  )
  
  lars_m2 <- glm("case ~ linearpredictor", data = merged_data, family = binomial)
  if (calc_r2 == TRUE) {
    TOGGLE = (class(lars_m2)[1] == "lm" | class(lars_m2)[1] == "glm")
    if (!TOGGLE) {stop("Function written for `lm` and `glm` objects only")}
    null_m2 <- update(lars_m2, ~ 1)
    N_m2    <- nobs(lars_m2)
    m_m2    <- suppressWarnings(logLik(lars_m2, REML = FALSE))[1]
    n_m2    <- suppressWarnings(logLik(null_m2, REML = FALSE))[1]
    
    tmp_lars$r2_mcfadden      <- signif(1 - m_m2/n_m2, digits = 6)
    cs_m2                <- 1 - exp(-2/N_m2 * (m_m2 - n_m2))
    tmp_lars$r2_cox_and_snell <- signif(cs_m2, digits = 6)
    tmp_lars$r2_nagelkerke    <- signif(cs_m2 / (1 - exp((2 / N_m2) * n_m2)),
                                        digits = 6)
  }
  
  lars_m2_hl           <- hoslem.test(lars_m2$y, fitted(lars_m2), g = 10)
  tmp_lars$hl_chisq    <- lars_m2_hl$statistic[["X-squared"]]
  tmp_lars$hl_p        <- lars_m2_hl$p.value
  tmp_lars$brier_score <- BrierScore(lars_m2)
  
  if (pctile_or == TRUE) {
    # OR of extreme PRS
    # TOP 1, 2, 5, 10 and 25%
    print(glue("Calculating OR of extreme PheRS:"))
    for (prob in c(0.01, 0.02, 0.05, 0.1, 0.25)) {
      print(glue("Percentile = {prob}..."))
      tmp_lars <- c(tmp_lars, getTopEffects(prob, pred = phers_name,
                                            cov = covars, dat = merged_data))
    }
    
    # determine percentile with at least 80% power
    #   (or the most powered percentiles)
    ptest    <- sapply(seq(0.005, 0.5, by = 0.005),
                       \(x) {getTopPower(prob = x, dat = merged_data,
                                         pred = phers_name)})
    presults <- data.table('Top' = seq(0.005, 0.5, by = 0.005), t(ptest))
    maxp     <- presults$Top[which(presults$h == max(presults$h,
                                                     na.rm = TRUE))[1]]
    
    p80 <- presults$Top[which(presults$power >= 0.80)]
    underpowered <- F
    if (length(p80) == 0) {
      p80 <- maxp
      underpowered <- T
    } else {
      p80 <- p80[1]
    }
    
    prs.p80        <- getTopEffects(p80, pred = phers_name,
                                    cov = covars, dat = merged_data)
    names(prs.p80) <- gsub(p80,"MinPower80", names(prs.p80))
    
    tmp_lars <- c(tmp_lars, "MinPower80" = p80, prs.p80,
                  "Top_Underpowered" = underpowered)
  }
  
  return(
    data.table(
      "stat" = names(tmp_lars),
      "max"  = c(unlist(tmp_max), rep(NA, length(unlist(tmp_lars)) - length(unlist(tmp_max)))),
      "lars" = unlist(tmp_lars)
    )
  )
  
}
