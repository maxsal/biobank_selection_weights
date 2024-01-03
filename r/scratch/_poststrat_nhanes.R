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
poststrat_nhanes <- function(
    int_data,
    id_var             = "id",
    last_entry_age_var = "AgeLastEntry",
    cancer_var         = "cancer",      
    cad_var            = "cad",         
    smoke_var          = "smoker",      
    diabetes_var       = "diabetes",    
    depression_var     = "depression",  
    female_var         = "female",      
    hypertension_var   = "hypertension",
    nhw_var            = "nhw",
    covs               = c("cad", "smoker", "diabetes", "female"),
    chop               = TRUE
) {

    # prep ---------------------------------------------------------------------
    if (!all(covs %in% c("cad", "smoker", "diabetes", "female", "depression", "hypertension", "cancer", "nhw"))) {
        stop("function only supports covs 'cad', 'smoker', 'diabetes', 'female', 'depression', 'hypertension', 'cancer', and 'nhw' right now")
    }

    Nobs    <- nrow(int_data)
    age_vec <- as.numeric(int_data[[last_entry_age_var]])
    which_between <- function(vec, mat) {
        which(data.table::between(
            x = vec, lower = mat[["lower"]],
            upper = mat[["upper"]]
        ))
    }

    # load and prep nhanes data ------------------------------------------------
    nhanes <- download_nhanes_data(datasets = c("DEMO", "BMX", "SMQ", "DIQ", "MCQ", "BPQ", "DPQ"))
    phanes <- prepare_nhanes_data(nhanes)
    phanes[, `:=` (
        smoker = ifelse(smoking_former == 1 | smoking_current == 1, 1, 0),
        age_bin      = cut(age, c(seq(0, 80, by = 10), 150), right = FALSE)
    )]
    setnames(phanes, "nhanes_nhw", "nhw")

    nhanes_design <- svydesign(
        id      = ~psu_nhanes,
        strata  = ~strata_nhanes,
        weights = ~weight_nhanes,
        nest    = TRUE,
        data    = phanes
    )

    # get nhanes proportions ---------------------------------------------------
    props <- list()
    for (i in c("female", "diabetes", "cad", "cancer", "bmi_obese", "depression", "hypertension", "nhw", "smoker")) {
        props[[i]] <- svyprop("age_bin", i, nhanes_design)
    }
    props[[length(props) + 1]] <- svyprop("age_bin", svydesign = nhanes_design)

    cprop <- Reduce(\(x, y) merge.data.table(x, y, by = "age_bin"), props)

    cprop[, `:=`(
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
            x     = tmp[, 1],
            y     = c(0, tmp[, 2]),
            right = FALSE
        )
    }

    func_pops_internal <- list()
    for (i in seq_along(covs)) {
        func_pops_internal[[i]] <- get_func_int(
            out_var = covs[i],
            x_var   = last_entry_age_var,
            data    = int_data
        )
        names(func_pops_internal)[i] <- gsub("_prop", "", covs[i])
    }
    func_pops_internal[[last_entry_age_var]] <- stepfun(
        x = seq(0, 85, 5),
        y = as.numeric(c(
            0,
            table(
                cut(
                    age_vec,
                    c(seq(0, 85, by = 5), 150),
                    right = FALSE
                ),
                useNA = "ifany"
            ) / Nobs)
        ), right = FALSE)
    
    # population probabilities ---------------------------------------------------
    population <- func_pops[["age_bin"]](age_vec)
    for (i in seq_along(covs)) {
        population <- population * 
            fifelse(
                int_data[[covs[i]]] == 1,
                func_pops[[covs[i]]](age_vec),
                1 - func_pops[[covs[i]]](age_vec)
            )
    }

    # internal probabilities -----------------------------------------------------
    internal <- func_pops_internal[[last_entry_age_var]](age_vec)
    for (i in seq_along(covs)) {
        internal <- internal * 
            fifelse(
                int_data[[covs[i]]] == 1,
                func_pops_internal[[covs[i]]](age_vec),
                1 - func_pops_internal[[covs[i]]](age_vec)
            )
    }

    # process --------------------------------------------------------------------
    if (chop == TRUE) {
        poststrat_weights <- chopr(population / internal)
    } else {
        poststrat_weights <- population / internal
    }

    poststrat_weights <- (Nobs * poststrat_weights) / sum(poststrat_weights, na.rm = TRUE)

    return(
        data.table(
            id        = int_data[[id_var]],
            ps_weight = poststrat_weights
        )
    )

}

setnames(mgi, "nhanes_nhw", "nhw")

tmp <- poststrat_nhanes(int_data = mgi)
tmp <- poststrat_nhanes(int_data = mgi, covs = c("nhw", "cad", "diabetes" "smoker"))

tmp2 <- merge.data.table(
    mgi,
    tmp,
    by = "id"
)

glm(cancer ~ female, data = tmp2, family = "binomial")
glm(cancer ~ female, data = tmp2, family = "binomial", weights = ps_weight)

tmp_dsn <- svydesign(
    id      = ~1,
    weight  = ~ps_weight,
    data    = tmp2[!is.na(ps_weight), ]
)

mean(tmp2[, female], na.rm = TRUE)
svymean(~female, tmp_dsn)
svymean(~female, nhanes_design)

mean(tmp2[, nhanes_nhw], na.rm = TRUE)
svymean(~nhanes_nhw, tmp_dsn)
svymean(~nhw, nhanes_design)[[1]]

tmp_fn <- function(data, var) {
    data.table(
        "variable"  = var,
        "raw"       = mean(data[[var]], na.rm = TRUE),
        "poststrat" = svymean(as.formula(paste0("~", var)), svydesign(~1, ~ps_weight, data = data[!is.na(ps_weight), ]))[[1]],
        "nhanes"    = svymean(as.formula(paste0("~", var)), nhanes_design, na.rm = TRUE)[[1]]
    )
}

tmp_fn(data = tmp2, var = "cad")
tmp_fn(data = tmp2, var = "smoker")
tmp_fn(data = tmp2, var = "diabetes")
tmp_fn(data = tmp2, var = "female")

lapply(
    c("nhw", "cad", "cancer", "smoker", "diabetes", "female"),
    \(x) tmp_fn(var = x, data = tmp2)
) |>
rbindlist()

mean(tmp2[, AgeLastEntry], na.rm = TRUE)
svymean(~AgeLastEntry, tmp_dsn)
svymean(~age, nhanes_design)[[1]]

svytable(~female+cancer+nhw, tmp_dsn)

svytable(~female, nhanes_design)

tmp_table <- svytable(~nhw + cad + cancer + smoker + diabetes + female, nhanes_design)
dt_table <- as.data.table(tmp_table)
prop_table <- as.data.table(prop.table(tmp_table))
tmp_table <- tmp_table[N > 0]

prop.table(tmp_table)

mgi_design <- svydesign(
    id      = ~1,
    data    = mgi[(!(is.na(nhw)) | is.na(cad) | is.na(cancer) | is.na(smoker) | is.na(diabetes) | is.na(female)), ]
)
mgi_post_design <- postStratify(mgi_design, ~ nhw + cad + cancer + smoker + diabetes + female,
                                tmp_table)


mgi_prop_table <- table(female = mgi[, female], cancer = mgi[, cancer], nhw = mgi[, nhw], cad = mgi[, cad]) |>
    prop.table() |>
    as.data.table()

num_cols <- c("female", "cancer", "nhw", "cad")

tmp2 <- merge.data.table(
    mgi,
    mgi_prop_table[, (num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols],
    by = c("female", "cancer", "nhw", "cad")
)

mgi_prop_table[, length(unique(N))]
tmp2[, length(unique(N))]
