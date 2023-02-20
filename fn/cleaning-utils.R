# replace NA with value --------------------------------------------------------
## inspiration from: https://stackoverflow.com/questions/7235657/fastest-way-to-replace-nas-in-a-large-data-table
replace_missing = function(data, cols = NULL, new_value = 0) {
  if (!is.null(cols)) {
    if ( !all(cols %in% names(data)) ) {
      missing_cols <- col[!(cols %in% names(data))]
      cli::cli_alert_warning("Some names in cols argument not present in data set: {paste0(missing_cols, collapse = ', ')}")
    } 
    for (i in cols) { data[is.na(get(i)), (i) := new_value] }
  } else {
    cli::cli_alert_info("Replacing missing in *all* columns with: {new_value}")
    for (i in names(data)) { data[is.na(get(i)), (i) := new_value] }
  }
}

# replace an existing value with NA --------------------------------------------
make_missing = function(data, cols = NULL, old_value = "Unknown") {
  if (!is.null(cols)) {
    if ( !all(cols %in% names(data)) ) {
      missing_cols <- col[!(cols %in% names(data))]
      cli::cli_alert_warning("Some names in cols argument not present in data set: {paste0(missing_cols, collapse = ', ')}")
    } 
    for (i in cols) { data[get(i) == old_value, (i) := NA] }
  } else {
    cli::cli_alert_info("Replacing '{old_value}' in *all* columns with missing (NA)")
    for (i in names(data)) { data[get(i) == old_value, (i) := NA] }
  }
}

# quickly generate an age group prevalence table
age_grp_table <- function(
    lower_ages,
    upper_char = "+",
    upper_val = 150,
    num_vec,
    num_var_name = "prevalences"
) {
  out <- data.table(
    group      = paste0(lower_ages, c(paste0("-", lower_ages[-1] - 1), upper_char)),
    lower      = lower_ages,
    upper      = c(lower_ages[-1] - 1, upper_val),
    num_var    = num_vec
  )
  setnames(out, "num_var", num_var_name)
  return(out)
}

## helper for truncating probabilities
chopr <- function(x) {
  quant2.5  <- quantile(x, probs = 0.025, na.rm = TRUE)
  quant97.5 <- quantile(x, probs = 0.975, na.rm = TRUE)
  x[x < quant2.5]  <- quant2.5
  x[x > quant97.5] <- quant97.5
  return(x)
}
