# functions to assist in evaluation of cooccurence and phers results
# author: max salvatore
# date:   20221212

suppressPackageStartupMessages({
  library(data.table)
  library(progress)
  library(glue)
  library(logistf)
  library(ResourceSelection)
  library(DescTools)
  library(pROC)
})

# log10toP: converts log10 p-values to decimal format --------------------------
log10toP <- function(log10P) {
  log10P <- abs(as.numeric(log10P))
  if (is.na(log10P)) {
    return(NA)
  }
  if (log10P > 300) {
    part1 <- log10P %/% 100 * 100
    part2 <- log10P - part1
    if (part2 != 0) {
      P <- format(signif(10^-part2, 3), scientific = T)
      P <- paste(as.numeric(gsub("e-.+", "", P)), "e-",
                 as.numeric(gsub(".+-", "", P), sep = "") + part1, sep = "")
    } else {
      P <- paste("1e-", part1, sep = "")
    }
  } else {
    P <- signif(10^-log10P, 3)
  }
  return(as.character(P))
}

# getTopEffects (from Lars) ----------------------------------------------------
getTopEffects <- function(prob, pred, cov, dat) {
  add_cov <- ifelse(length(cov) == 0, "",
                    paste0(" + ", paste0(cov, collapse = " + ")))
  riskBin <- glue::glue("Top_{prob}")
  dat[[riskBin]] <- ifelse(
    dat[[pred]] >= quantile(dat[case == 0, ][[pred]], probs = 1 - prob),
    1,
    0)
  tmp_vars <- c("case", riskBin, cov)
  fitAsso2 <- logistf(
    as.formula(paste0("case ~ ", riskBin, add_cov)),
    data = na.omit(dat[, ..tmp_vars]),
    control = logistf.control(maxit = 1000))
  getValues(fitAsso2, riskBin)
}

# getTopPower (from Lars) ------------------------------------------------------
getTopPower <- function(prob, dat, pred) {
  riskBin               <- glue("Top_{prob}")
  dat[[riskBin]] <- ifelse(
    dat[[pred]] >= quantile(dat[case == 0, ][[pred]], probs = 1 - prob),
    1,
    0)
  ptable           <- table('Top' = dat[[riskBin]], 'trait' = dat[["case"]])
  if (nrow(ptable) != 2) return(c('power' = NA,'h' = NA))
  ntop             <- sum(ptable[, 2])
  nrest            <- sum(ptable[, 1])
  ptop             <- ptable[2, 2] / ntop
  prest            <- ptable[2, 1] / nrest
  h                <- pwr::ES.h(ptop, prest)
  unlist(pwr::pwr.2p2n.test(h         = h,
                            n1        = ntop,
                            n2        = nrest,
                            sig.level = 0.05,
                            power     = NULL)[c("power", "h")])
}

# getValues (from Lars) --------------------------------------------------------
getValues <- function(fitAsso, predictor) {
  sebeta        <- sqrt(diag(vcov(fitAsso)))
  names(sebeta) <- fitAsso$terms
  BETA          <- round(fitAsso$coefficient[predictor],4)
  SEBETA        <- round(sebeta[predictor],4)
  
  vars        <- diag(fitAsso$var)
  names(vars) <- fitAsso$terms
  LOG10P      <- -pchisq((fitAsso$coefficient[predictor]^2/vars[predictor]),
                         1, lower.tail = F,log = T) / log(10)
  P           <- log10toP(LOG10P)
  
  OR         <- round(exp(fitAsso$coefficient[predictor]),4)
  CI1        <- round(exp(fitAsso$ci.lower[predictor]),4)
  CI2	       <- round(exp(fitAsso$ci.upper[predictor]),4)
  out        <- c(BETA,SEBETA,P,OR,CI1,CI2,signif(LOG10P,4))
  names(out) <- paste0(predictor, "_" ,
                       c("BETA", "SEBETA", "P", "OR", "CI1", "CI2", "LOG10P"))
  return(out)
}

# evaluate_phers ---------------------------------------------------------------
evaluate_phers <- function(
    eval_pim,
    predictor = "phers",
    phers_phes,
    covariates = c("age_at_threshold", "female"),
    out_phe,
    calc_r2 = TRUE,
    pctile_or = TRUE
) {
  
  # initialize output list object
  out <- list()
  
  out$trait      <- out_phe
  out$predictor  <- predictor
  out$phecodes   <- paste0(phers_phes, collapse = ", ")
  out$n_phecodes <- length(phers_phes)
  
  add_cov <- ifelse(length(covariates) == 0, "",
                    paste0(" + ", paste0(covariates, collapse = " + ")))
  
  # m1: adjusted GLM model
  message("m1: adjusted GLM model")
  m1   <- stats::glm(paste0("case ~ ", predictor, add_cov),
                     family = binomial(), data = eval_pim)
  vars <- diag(vcov(m1))
  
  if (coef(summary(m1))[predictor, 4] == 0) {
    out$p_value <- log10toP(-pchisq((m1$coefficient[predictor]^2 / 
                                       vars[predictor]),
                                    1,
                                    lower.tail = F,
                                    log = T) / log(10))
  } else {
    out$p_value <- coef(summary(m1))[predictor, 4]
  }
  
  m1_confint <- confint(m1)
  
  # basic performance as predictor
  out$beta   <- m1$coefficient[[predictor]]
  out$sebeta <- sqrt(diag(vcov(m1))[[predictor]])
  out$or     <- exp(m1$coefficient[[predictor]])
  out$or_ci  <- glue::glue("{round(exp(m1_confint[predictor, 1]), 5)}, ",
                           "{round(exp(m1_confint[predictor, 2]), 5)}")
  out$log10p <- -pchisq((m1$coefficient[[predictor]] ^ 2 / vars[[predictor]]),
                        1,
                        lower.tail = FALSE,
                        log = TRUE)/log(10)
  
  # auc
  message("m1 AUC")
  m1_auc     <- pROC::roc(response = eval_pim[["case"]],
                          predictor = eval_pim[[predictor]],
                          family = binomial(),
                          ci = TRUE)
  out$auc    <- m1_auc$auc
  out$auc_ci <- paste0(m1_auc$ci[c(1, 3)], collapse = ", ")
  
  # m2: unadjusted GLM model
  message("m2: unadjusted GLM model")
  m2      <- glm("case ~ phers", family = binomial(), data = eval_pim)
  
  ### calculate R2
  if (calc_r2 == TRUE) {
    TOGGLE = (class(m2)[1] == "lm" | class(m2)[1] == "glm")
    if (!TOGGLE) {stop("Function written for `lm` and `glm` objects only")}
    null_m2 <- update(m2, ~ 1)
    N_m2    <- nobs(m2)
    m_m2    <- suppressWarnings(logLik(m2, REML = FALSE))[1]
    n_m2    <- suppressWarnings(logLik(null_m2, REML = FALSE))[1]
    
    out$r2_mcfadden      <- signif(1 - m_m2/n_m2, digits = 6)
    cs_m2                <- 1 - exp(-2/N_m2 * (m_m2 - n_m2))
    out$r2_cox_and_snell <- signif(cs_m2, digits = 6)
    out$r2_nagelkerke    <- signif(cs_m2 / (1 - exp((2 / N_m2) * n_m2)),
                                   digits = 6)
  }
  
  m2_hl           <- hoslem.test(m2$y, fitted(m2), g = 10)
  out$hl_chisq    <- m2_hl$statistic[["X-squared"]]
  out$hl_p        <- m2_hl$p.value
  out$brier_score <- BrierScore(m2)
  
  ### OR CALCULATION
  ### cross-validation
  print(glue("Performing cross validation..."))
  # cross validation with 5 random folds ---
  G                <- 5
  eval_pim$folds   <- createFolds(y = 1:nrow(eval_pim), k = G, list = F)
  linearpredictors <- list()
  
  for (g in 1:G) {
    testdata  <- eval_pim[eval_pim$folds == g, ]
    traindata <- eval_pim[eval_pim$folds != g, ]
    
    # fit
    formula.fit <- glue("case ~ {predictor}")
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
  }
  
  eval_pim <- merge.data.table(
    eval_pim, rbindlist(linearpredictors),
    by = "id"
  )
  
  print(glue("Cross validation complete. Sit tight..."))
  
  if (pctile_or == TRUE) {
    # OR of extreme PRS
    # TOP 1, 2, 5, 10 and 25%
    print(glue("Calculating OR of extreme PheRS:"))
    for (prob in c(0.01, 0.02, 0.05, 0.1, 0.25)) {
      print(glue("Percentile = {prob}..."))
      out <- c(out, getTopEffects(prob, pred = predictor,
                                  cov = covariates, dat = eval_pim))
    }
    
    # determine percentile with at least 80% power
    #   (or the most powered percentiles)
    ptest    <- sapply(seq(0.005, 0.5, by = 0.005),
                       \(x) {getTopPower(prob = x, dat = eval_pim,
                                         pred = predictor)})
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
    
    prs.p80        <- getTopEffects(p80, pred = predictor,
                                    cov = covariates, dat = eval_pim)
    names(prs.p80) <- gsub(p80,"MinPower80", names(prs.p80))
    
    out <- c(out, "MinPower80" = p80, prs.p80,
             "Top_Underpowered" = underpowered)
  }
  
  return(data.table(
    "stat"  = names(unlist(out)),
    "value" = unlist(out)
  ))
  
}

# calculate_phers --------------------------------------------------------------
calculate_phers <- function(
    pim,                 # phecode indicator matrix
    res,                 # results matrix - phecode, beta, p-value
    method,              # tophits or pwide_sig
    tophits_n    = 50,   # n of hits to select based on p-value
    bonf_tests   = NULL, # denominator of bonferroni correction
    phers_name   = NULL,
    reverse_code = FALSE # reverse code negatives to keep phers above 0?
    ) {
  
  # check that phers_name argument is specified
  if (is.null(phers_name)) {
    stop("'phers_name' argument is not specified!")
  }
  
  # check that columns exist in res data
  res_data_cols <- c("phecode", "beta", "p_value")
  if (!any(res_data_cols %in% names(res))) {
    stop(glue("Missing columns in 'res': ",
              "{glue_collapse(res_data_cols[!(res_data_cols %in% names(res))]",
              ", sep = ', ')}"))
  }
  
  ## select significant hits
  if (method == "tophits") {
    phers_hits <- res[order(p_value)][1:tophits_n]
  }
  if (method == "pwide_sig") {
    if (is.null(bonf_tests)) {
      phers_hits <- res[p_value < 0.05/nrow(pim)]
    } else {
      phers_hits <- res[p_value < 0.05/bonf_tests]
    }
  }
  
  out <- data.table::copy(pim)
  
  out[, phers := 0]
  message("calculating phers...")
  pb <- progress_bar$new(total = length(phers_hits[, phecode]),
                         format = "[:bar] :percent eta: :eta")
  
  if (reverse_code == FALSE) {
    for (i in phers_hits[, phecode]) {
      out[, phers := phers + (phers_hits[phecode == i, beta] * get(i))]
      pb$tick()
    }
  } else {
  for (i in phers_hits[, phecode]) {
      out[, phers := phers + (
        phers_hits[phecode == i, beta] * get(i) * 
          as.numeric(phers_hits[phecode == i, beta] > 0)) - (
          phers_hits[phecode == i, beta] * (1 - get(i)) * 
            as.numeric(phers_hits[phecode == i, beta] < 0)
        )]
      pb$tick()
    }
  }
  
  setnames(out, old = "phers", new = phers_name)
  keep_cols <- c("id", phers_name)

  return(list(
    n_phecodes = length(phers_hits[, phecode]),
    method     = method,
    phecodes   = phers_hits,
    data       = out[, ..keep_cols],
    phers_name = phers_name
    ))
  
}

