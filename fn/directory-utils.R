## function that checks for presence of folders corresponding to data date
## and, optionally, outcome phecode and makes missing folders if absent
check_folder_structure <- function(cohort = "mgi", data_version, outcome_phecode) {
  
  created <- 0
  
  # data -----------------------------------------------------------------------
  if (!dir.exists( glue::glue("data/private/{cohort}/{data_version}/X{gsub('X', '', outcome_phecode)}") )) {
    dir.create(glue::glue("data/private/{cohort}/{data_version}/X{gsub('X', '', outcome_phecode)}"), recursive = TRUE)
    cli::cli_alert_info("data/private/{cohort}/{data_version}/X{gsub('X', '', outcome_phecode)}/ folder created")
    created <- created + 1
  }
  
  if ( !dir.exists( glue::glue("results/{cohort}/{data_version}/") ) ) {
    dir.create(glue::glue("results/{cohort}/{data_version}/"), recursive = TRUE)
    cli::cli_alert_info("results/{cohort}/{data_version}/ folder created")
    created <- created + 1
  }
  
  # results --------------------------------------------------------------------
  if (!is.null(outcome_phecode)) {
    if ( !dir.exists( glue::glue("results/{cohort}/{data_version}/X{outcome_phecode}/") ) ) {
      dir.create(glue::glue("results/{cohort}/{data_version}/X{outcome_phecode}/"), recursive = TRUE)
      cli::cli_alert_info("results/{cohort}/{data_version}/X{outcome_phecode}/ folder created")
      created <- created + 1
    }
  }
  
  # summary --------------------------------------------------------------------
  if (created == 0) {
    if (is.null(outcome_phecode)) {
      cli::cli_alert_success("file structure already exists for {cohort} version {data_version}")
    } else {
      cli::cli_alert_success("file structure already exists for {cohort} version {data_version} and outcome phecode {outcome_phecode}")
    }
  }
  
}
