suppressPackageStartupMessages({
  require(SPAtest)
  require(logistf)
})

### quick_mod -----------
# function for obtaining beta and p-values using either SPAtest or logistf
# SPAtest is *much* faster but logistf allows for inclusion of weights
quick_cooccur_mod <- function(
    dat,
    covs       = c("age_at_threshold", "female", "length_followup"),
    ex_code    = "X157",
    mod_type   = "glm",
    weight_var = NULL
) {
  
  ### SPAtest
  if (mod_type == "SPAtest") {
    mod <- ScoreTest_SPA(
      genos       = t(dat[[ex_code]]),
      pheno       = dat[["case"]],
      cov         = dat[, ..covs],
      method      = "fastSPA",
      alpha       = as.numeric(0.05),
      beta.out    = TRUE,
      beta.Cutoff = 1
    )
    
    data.table(
      phecode = ex_code,
      beta    = mod$beta,
      se_beta = mod$SEbeta,
      p_value = mod$p.value
    )
  } 
  
  ### logistf
  if (mod_type == "logistf") {
    if (!is.null(weight_var)) {
      wgts <- dat[[weight_var]]
      mod <- logistf(as.formula(paste0("case ~ ", ex_code, " + ", paste0(covs, collapse = " + "))),
                     data    = dat,
                     weights = wgts)
    } else {
      mod <- logistf(paste0("case ~ ", ex_code, " + ", paste0(covs, collapse = " + ")), data = dat)
    }
    
    data.table(
      phecode = ex_code,
      beta    = mod$coefficients[[ex_code]],
      se_beta = sqrt(diag(vcov(mod)))[[ex_code]],
      p_value = mod$prob[[ex_code]],
      log10p  = log10(mod$prob[[ex_code]])
    )
    
  }
  
  ### glm
  if (mod_type == "glm") {
    if (!is.null(weight_var)) {
      wgts <- dat[[weight_var]]
      mod <- glm(as.formula(paste0("case ~ ", ex_code, " + ", paste0(covs, collapse = " + "))),
                 data    = dat,
                 weights = wgts,
                 family  = binomial())
    } else {
      mod <- glm(paste0("case ~ ", ex_code, " + ", paste0(covs, collapse = " + ")), data = dat, family = binomial())
    }
    
    data.table(
      phecode = ex_code,
      beta    = coef(mod)[[ex_code]],
      se_beta = coef(summary(mod))[ex_code, 2],
      p_value = coef(summary(mod))[ex_code, 4],
      log10p  = log10(coef(summary(mod))[ex_code, 4])
    )
    
  }
  
}

### output_mgi_cooccur_results -----------
# a function for quickly iterating cooccur analyses over a phecode indicator
# matrix (developed using MGI data)
# provides options for using SPAtest and logistf
# SPAtest is *much* faster but logistf allows for inclusion of weights
output_cooccurrence_results <- function(
    pim_data,
    cov_data,
    covariates,
    t_thresh,
    all_phecodes = paste0("X", pheinfo[, phecode]),
    model_type   = "logistf",
    w_data       = NULL,
    w_var        = NULL,
    ncore        = parallel::detectCores() / 4,
    parallel     = TRUE
    ) {
  
  # 1. identify analytic phecodes
  possible_phecodes    <- names(pim_data)[names(pim_data) %in% all_phecodes]
  phecodes_to_consider <- melt(pim_data[, ..possible_phecodes][, lapply(.SD, \(x) sum(x, na.rm = TRUE))], variable.name = "phecode", value.name = "n", id.vars = character())[n >= 10, as.character(phecode)]
  
  # 2. merge covariates
  if (is.null(w_data) | is.null(w_var)) {
    merged <- merge.data.table(
      pim_data[, !c("case")],
      cov_data,
      by = "id",
      all.x = TRUE
    )[, age_at_threshold := round(get(paste0("t", t_thresh, "_threshold")) / 365.25, 1)][]
  } else {
    merged <- Reduce(f = \(x, y) merge.data.table(x, y, by = "id", all.x = TRUE),
                     x = list(pim_data[, !c("case")], cov_data, w_data)
                     )[, age_at_threshold := round(get(paste0("t", t_thresh, "_threshold")) / 365.25, 1)][]
  }
  
  # 3. run analyses
  if (parallel == FALSE) {
    out <- list()
    cli_progress_bar(name = glue("t = {t_thresh} threshold"),
                     total = length(phecodes_to_consider))
    for (i in seq_along(phecodes_to_consider)) {
      out[[i]] <- quick_cooccur_mod(
        dat        = merged,
        covs       = covariates,
        ex_code    = phecodes_to_consider[i],
        mod_type   = model_type,
        weight_var = w_var
      )
      cli_progress_update()
    }
    out <- rbindlist(out)
  }
  
  if (parallel == TRUE) {
    require(doMC)
    registerDoMC(cores = ncore)
    columns <- phecodes_to_consider
    cols    <- seq_along(phecodes_to_consider)
    output  <- foreach(i = cols) %dopar% {
      quick_cooccur_mod(
        dat        = merged,
        covs       = covariates,
        ex_code    = phecodes_to_consider[i],
        mod_type   = model_type,
        weight_var = w_var
      )
    }
    out <- rbindlist(output, use.names = TRUE, fill = TRUE)
  }
  
  return(out)
  
}

