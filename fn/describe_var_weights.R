# summarize unweighted and weighted MGI and NHANES variables by weighting design
# author: max salvatore
# date:   20230110

describe_var_weights <- function(
    mgi_data,
    nhanes_data,
    mgi_var_name,
    nhanes_var_name,
    nhanes_var_val = NULL,
    designs
) {

  unweighted <- data.table("mean" = mean(mgi_data[[mgi_var_name]], na.rm = TRUE))[, weight := "unweighted"]
  
  out <- list()
  for (i in 1:length(designs)) {
    if (grepl("nhanes", names(designs)[i])) {
      if (!is.null(nhanes_var_val)) {
        out[[i]] <- as.data.table(
          svymean(as.numeric(nhanes_data[[nhanes_var_name]] == nhanes_var_val),
                  design = designs[[i]],
                  na.rm = TRUE)
        )[, weight := names(designs)[i]]
      } else {
        out[[i]] <- as.data.table(
          svymean(nhanes_data[[nhanes_var_name]],
                  design = designs[[i]],
                  na.rm = TRUE)
        )[, weight := names(designs)[i]]
      }
    } else {
      out[[i]] <- as.data.table(
        svymean(mgi_data[[mgi_var_name]], design = designs[[i]])
      )[, weight := names(designs)[i]]
    }
  }
  
  rbindlist(
    list(unweighted, rbindlist(out)),
    fill = TRUE, use.names = TRUE
  )[, variable := mgi_var_name][]
}

