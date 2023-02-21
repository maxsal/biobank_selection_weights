suppressPackageStartupMessages({
  require(data.table)
  require(glue)
})
generate_restricted_phenome <- function(phe_data, threshold, cases, outcome_phe) {
  
  if (!grepl("X", outcome_phe)) {
    outcome_phe <- glue::glue("X{outcome_phe}")
  }
  
  out <- data.table::dcast(
    unique(phe_data[get(glue::glue("t{threshold}_indicator")) == 1, ][dsb < get(glue::glue("t{threshold}_threshold")), ][, .(id, phecode = paste0("X", phecode))]),
    id ~ phecode,
    value.var     = "phecode",
    fun.aggregate = length,
    fill          = 0
  )
  
  out[, case := data.table::fifelse(id %in% cases, 1, 0)]
  if (outcome_phe %in% names(out)) {
    out[, c(outcome_phe) := NULL]
  }
  
  return(out)
  
}