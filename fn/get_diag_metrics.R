# getting diagnostic record metrics
get_diag_metrics <- function(cohort, icd_phecode_data, data_version = NULL, force = FALSE) {
  
  if (tolower(cohort) == "mgi") {
    if( file.exists( paste0("data/", data_version, "/processed/", cohort, "_diag_metrics.txt") ) & force == FALSE ){
      cli::cli_alert_info("congrats! {cohort}_diag_metrics already exists")
      out <- data.table::fread(paste0("data/", data_version, "/processed/", cohort, "diag_metrics.txt"))
    } else {
      first_phe <- icd_phecode_data[ icd_phecode_data[, .I[which.min(DaysSinceBirth)], by = "IID"]$V1 ][, .(id = IID, first_dsb = as.numeric(DaysSinceBirth))]
      last_phe <- icd_phecode_data[ icd_phecode_data[, .I[which.max(DaysSinceBirth)], by = "IID"]$V1 ][, .(id = IID, last_dsb = as.numeric(DaysSinceBirth))]
      n_encounters <- unique(icd_phecode_data[, .(id = IID, dsb = DaysSinceBirth)])[, .N, by = "id"][, .(id, n_encounters = N)]
      
      out <- Reduce(merge.data.table, list(first_phe, last_phe, n_encounters))
      out[, `:=` (
        age_at_first_diagnosis = round(first_dsb / 365.25, 1),
        age_at_last_diagnosis  = round(last_dsb / 365.25, 1)
      )][, `:=` (
        length_followup = round((age_at_last_diagnosis - age_at_first_diagnosis)),
        encounters_per_year = round(n_encounters / (age_at_last_diagnosis - age_at_first_diagnosis), 1)
        )]
      out[is.infinite(encounters_per_year), encounters_per_year := NA]
      data.table::fwrite(x = out, file = paste0("data/", data_version, "/processed/", cohort, "_diag_metrics.txt"))
    }
  }
  
  if (tolower(cohort) == "ukb") {
    if( file.exists( paste0("data/", data_version, "/processed/mgi_diag_metrics.txt") ) & force == FALSE ) {
        cli::cli_alert_info("congrats! {cohort}_diag_metrics already exists")
        out <- data.table::fread(paste0("data/", data_version, "/processed/", cohort, "diag_metrics.txt"))
    } else {
      first_phe          <- icd_phecode_data[ icd_phecode_data[, .I[which.min(dsb)], by = "id"]$V1 ][, .(id, first_dsb = dsb)]
      last_phe           <- icd_phecode_data[ icd_phecode_data[, .I[which.max(dsb)], by = "id"]$V1 ][, .(id, first_dsb = dsb)]
      n_unique_diagnoses <- unique(icd_phecode_data[, .(id, phecode)])[, .N, by = "id"][, .(id, n_unique_diagnoses = N)]
      
      out <- Reduce(merge.data.table, list(first_phe, last_phe, n_unique_diagnoses))
      out[, `:=` (
        age_at_first_diagnosis = round(first_dsb / 365.25, 1),
        age_at_last_diagnosis = round(last_dsb/ 365.25, 1)
      )][, `:=` (
        length_followup = round((age_at_last_diagnosis - age_at_first_diagnosis))
      )]
      data.table::fwrite(x = out, file = paste0("data/", data_version, "/processed/", cohort, "_diag_metrics.txt"))
      }
  }
  
  return(out[])
  
}
