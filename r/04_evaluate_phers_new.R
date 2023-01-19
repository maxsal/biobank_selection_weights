# evaluate a phers
# author: max salvatore
# date:   20230118

# 1. libraries, functions, and options (outcome agnostic) ----------------------
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(caret)
  library(purrr)
  library(progress)
  library(pROC)
  library(glue)
  library(logistf)
  library(cli)
  library(optparse)
})

set.seed(61787)

for (i in list.files("fn/", full.names = TRUE)) source(i) # load functions

# optparse list ----------------------------------------------------------------
option_list <- list(
  make_option("--outcome", type = "character", default = "157",
              help = "Outcome phecode [default = 157]"),
  make_option("--mgi_version", type = "character", default = "20210318",
              help = "Version of MGI data [default = 20210318]"),
  make_option("--ukb_version", type = "character", default = "20221117",
              help = "Version of UKB data [default = 20221117]"),
  make_option("--time_threshold", type = "numeric", default = "0",
              help = glue("Time threshold for the phenome data ",
                          "[default = 0]")),
  make_option("--tophits_n", type = "numeric", default = "50",
              help = glue("Number of top hits to use in top hits PheRS ",
                          "[default = 50]")),
  make_option("--pctile_or", type = "logical", default = "TRUE",
              help = glue("Perform percentile-based OR diagnostics ",
                          "which is relatively time consuming ",
                          "[default = TRUE]"))
)

####!!! ADD arguments specifying covariates

parser <- OptionParser(usage="%prog [options]", option_list = option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

# 2. specifications (specifies outcome) ----------------------------------------
mgi_version    <- opt$mgi_version  # mgi phenome version
ukb_version    <- opt$ukb_version  # ukb phenome version
outcome        <- opt$outcome      # outcome phecode
time_threshold <- opt$time_threshold

## extract file paths
file_paths <- get_files(mgi_version = opt$mgi_version,
                        ukb_version = opt$ukb_version)

# read data --------------------------------------------------------------------
## mgi
### phers
#### naive
mgi_naive_phers <- fread(glue("results/mgi/{opt$mgi_version}/",
                              "X{gsub('X', '', opt$outcome)}/phers/",
                              "mgi_naive_phers_t{opt$time_threshold}.txt"))

mgi_naive_phers_info <- list()
for (i in grep(
  list.files(glue("results/mgi/{opt$mgi_version}/",
                "X{gsub('X', '', opt$outcome)}/phers"),
             full.names = TRUE),
  pattern = glue("mgi_phers_t{opt$time_threshold}"),
  value = TRUE
)) {
  mgi_naive_phers_info[[i]] <- readRDS(i)
}

### covariates
mgi_covariates <- fread(glue("data/private/mgi/{opt$mgi_version}/",
                             "X{gsub('X', '', opt$outcome)}/",
                             "matched_covariates.txt"))[
                               id %in% mgi_naive_phers[, id]
                               ][
                                 , age_at_threshold := round(get(glue("t{opt$time_threshold}_threshold")) / 365.25, 1)
                               ]

mgi_merged <- merge.data.table(
  mgi_naive_phers,
  mgi_covariates[, !c("case")],
  by = "id"
)

## ukb
### phers
ukb_naive_phers <- fread(glue("results/ukb/{opt$ukb_version}/",
                              "X{gsub('X', '', opt$outcome)}/phers/",
                              "ukb_naive_phers_t{opt$time_threshold}.txt"))

ukb_naive_phers_info <- list()
for (i in grep(
  list.files(glue("results/ukb/{opt$ukb_version}/",
                  "X{gsub('X', '', opt$outcome)}/phers"),
             full.names = TRUE),
  pattern = glue("ukb_phers_t{opt$time_threshold}"),
  value = TRUE
)) {
  ukb_naive_phers_info[[i]] <- readRDS(i)
}

### covariates
ukb_covariates <- fread(glue("data/private/ukb/{opt$ukb_version}/",
                             "X{gsub('X', '', opt$outcome)}/",
                             "matched_covariates.txt"))[
                               id %in% ukb_naive_phers[, id]
                             ][
                               , age_at_threshold := round(get(glue("t{opt$time_threshold}_threshold")) / 365.25, 1)
                             ]

ukb_merged <- merge.data.table(
  ukb_naive_phers,
  ukb_covariates[, !c("case")],
  by = "id"
)

## phenome
pheinfo <- fread("data/public/Phecode_Definitions_FullTable_Modified.txt",
                 colClasses = "character")

## evaluate

tmp_eval_fn <- function(
  merged_data,
  covars,
  outcome,
  phers_name,
  phers_info_list,
  calc_r2 = TRUE,
  pctile_or = TRUE
  ) {
  # initialize
  tmp_max <- NULL
  tmp_lars <- NULL
  
  tmp_max$trait     <- gsub("X", "", outcome)
  tmp_max$predictor <- phers_name
  print(grep(
    x       = names(phers_info_list),
    pattern = phers_name,
    value   = TRUE
  ))
  tmp_max$phecodes  <- phers_info_list[[grep(
    x       = names(phers_info_list),
    pattern = phers_name,
    value   = TRUE
  )]][["n_phecodes"]]
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

# setnames(mgi_merged, "phers0_t0_dm_bn", "phers_t0_dm_bn")

tmp_eval_fn(
  merged_data     = mgi_merged,
  covars          = c("age_at_threshold", "female"),
  outcome         = opt$outcome,
  phers_name      = "phers_t0_dm_bn",
  phers_info_list = mgi_naive_phers_info,
  pctile_or       = opt$pctile_or
)


