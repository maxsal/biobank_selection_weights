### quick_mod -----------
# function for obtaining beta and p-values using either SPAtest or logistf
# SPAtest is *much* faster but logistf allows for inclusion of weights
quick_cooccur_mod <- function(dat, covs = c("age_at_threshold", "female", "length_followup"), ex_code = "X157", mod_type = "logistf") {
  
  ### SPAtest
  if (mod_type == "SPAtest") {
    mod <- SPAtest::ScoreTest_SPA(
      genos       = t(dat[[ex_code]]),
      pheno       = dat[["case"]],
      cov         = dat[, ..covs],
      method      = "fastSPA",
      alpha       = as.numeric(0.05),
      beta.out    = TRUE,
      beta.Cutoff = 1
    )
    
    data.table::data.table(
      phecode = ex_code,
      beta    = mod$beta,
      se_beta = mod$SEbeta,
      p_value = mod$p.value
    )
  } 
  
  ### logistf
  if (mod_type == "logistf") {
    mod <- logistf::logistf(paste0("case ~ ", ex_code, " + ", paste0(covs, collapse = " + ")), data = dat)
    
    data.table::data.table(
      phecode = ex_code,
      beta    = mod$coefficients[[ex_code]],
      se_beta = sqrt(diag(vcov(mod)))[[ex_code]],
      p_value = mod$prob[[ex_code]]
    )
  }
  
}

### output_mgi_cooccur_results -----------
# a function for quickly iterating cooccur analyses over a phecode indicator
# matrix (developed using MGI data)
# provides options for using SPAtest and logistf
# SPAtest is *much* faster but logistf allows for inclusion of weights
output_cooccurrence_results <- function(pim_data, cov_data, covariates, t_thresh, all_phecodes = paste0("X", pheinfo[, phecode]), model_type = "logistf") {
  
  # 1. identify analytic phecodes
  possible_phecodes    <- names(pim_data)[names(pim_data) %in% all_phecodes]
  phecodes_to_consider <- melt(pim_data[, ..possible_phecodes][, lapply(.SD, \(x) sum(x, na.rm = TRUE))], variable.name = "phecode", value.name = "n", id.vars = character())[n >= 10, as.character(phecode)]
  
  # 2. merge covariates
  merged <- merge.data.table(
    pim_data[, !c("case")],
    cov_data,
    by = "id",
    all.x = TRUE
  )[, age_at_threshold := round(get(paste0("t", t_thresh, "_threshold")) / 365.25, 1)][]
  
  # 3. run analyses
  out <- list()
  pb <- progress::progress_bar$new(total = length(phecodes_to_consider),
                                   format = ":what [:bar] :percent eta: :eta")
  for (i in seq_along(phecodes_to_consider)) {
    out[[i]] <- quick_cooccur_mod(
      dat      = merged,
      covs     = covariates,
      ex_code  = phecodes_to_consider[i],
      mod_type = model_type
    )
    pb$tick(tokens = list(what = glue::glue("t = {t_thresh} threshold")))
  }
  
  out <- rbindlist(out)
  
  return(out)
  
}
