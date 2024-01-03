# libraries and such -----------------------------------------------------------
suppressPackageStartupMessages({
    library(ms)
    library(qs)
    library(data.table)
    library(tidyverse)
    library(survey)
    library(stringr)
})

source("fn/download_nhanes_data.R")
source("fn/prepare_nhanes_data.R")

# load data --------------------------------------------------------------------
nhanes <- download_nhanes_data(
    datasets = c("DEMO", "BMX", "SMQ", "DIQ", "MCQ", "BPQ", "DPQ")
)
phanes <- prepare_nhanes_data(nhanes)
phanes[, `:=` (
    smoking_ever = ifelse(smoking_former == 1 | smoking_current == 1, 1, 0),
    age_bin      = cut(age, c(seq(0, 80, by = 10), 150), right = FALSE)
)]

nhanes_design <- svydesign(
    id      = ~psu_nhanes,
    strata  = ~strata_nhanes,
    weights = ~weight_nhanes,
    nest    = TRUE,
    data    = phanes
)

# functions --------------------------------------------------------------------
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

# get proportions --------------------------------------------------------------
props <- list()
for (i in c("diabetes", "cad", "cancer", "bmi_obese", "nhanes_nhw","smoking_current")) {
    props[[i]] <- svyprop("age_bin", i, nhanes_design)
}
props[[length(props) + 1]] <- svyprop("age_bin", svydesign = nhanes_design)

cprop <- Reduce(\(x, y) merge.data.table(x, y, by = "age_bin"), props)

cprop[, `:=` (
    lower = str_extract_to_number(age_bin, "\\[(\\d+),"),
    upper = str_extract_to_number(age_bin, ",(\\d+)\\)")
)]

prop_vars <- grep("_prop", names(cprop), value = TRUE)

# create step functions --------------------------------------------------------
func_pops <- list()
for (i in seq_along(prop_vars)) {
    func_pops[[i]] <- stepfun(
        x = cprop[, lower],
        y = c(0, cprop[[prop_vars[i]]]),
        right = FALSE
    )
    names(func_pops)[i] <- gsub("_prop", "", prop_vars[i])
}

# get internal data step functions ---------------------------------------------
get_func_int <- function(out_var, x_var, FUN = mean, data) {
    tmp <- aggregate(
        as.formula(
            paste0(
                out_var, "~", x_var
            )
        ),
        FUN = FUN, data = data
    )
    stepfun(
        x = tmp[, 1],
        y = c(0, tmp[, 2],
        right = FALSE)
    )
}

func_pops_internal <- list()
