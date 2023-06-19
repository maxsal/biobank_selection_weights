# libraries and such -----------------------------------------------------------
suppressPackageStartupMessages({
    library(ms)
    library(qs)
    library(data.table)
    library(tidyverse)
    library(survey)
    library(stringr)
})

# helper functions -------------------------------------------------------------
source("fn/download_nhanes_data.R")
source("fn/prepare_nhanes_data.R")

svyprop <- function(strata_var, prop_var = NULL, svydesign, value_var = "N") {
    if (!is.null(prop_var)) {
        tmp_form <- as.formula(paste0("~", strata_var, "+", prop_var))
        tmp <- svytable(tmp_form, svydesign) |>
            as.data.table() |>
            dcast(as.formula(paste0(strata_var, "~", prop_var)), value.var = value_var)
        tmp[, prop := ifelse(`0`+`1` == 0, 0, `1`/(`0`+`1`))]
        keep <- c(strata_var, "prop")
        out <- tmp[, ..keep]
        setnames(out, "prop", paste0(prop_var, "_prop"))
    } else {
        tmp_form  <- as.formula(paste0("~", strata_var))
        tmp       <- svymean(tmp_form, svydesign, na.rm = TRUE)
        row_names <- gsub(strata_var, "", rownames(as.data.frame(tmp)))
        out       <- as.data.table(tmp)[, .(prop = mean)][, tmp := row_names]
        setnames(out, c("tmp", "prop"), c(strata_var, paste0(strata_var, "_prop")))
    }
    return(out[])
}

str_extract_to_number <- function(x, pattern) {
    str_extract(x, pattern) |>
        gsub(pattern, "\\1", x = _) |>
        as.numeric()
}

# poststrat function -----------------------------------------------------------
poststrat_nhanes <- function() {

    # load and prep nhanes data ------------------------------------------------
    nhanes <- download_nhanes_data(datasets = c("DEMO", "BMX", "SMQ", "DIQ", "MCQ", "BPQ", "DPQ"))
    phanes <- prepare_nhanes_data(nhanes)
    phanes[, `:=` (
        smoking_ever = ifelse(smoking_former == 1 | smoking_current == 1, 1, 0),
        age_bin      = cut(age, c(seq(0, 80, by = 10), 150), right = FALSE)
    )]


}