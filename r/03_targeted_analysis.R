# Conduct targeted hypertension ~ sex analysis in MGI and UKB
# author:  max salvatore
# date:    20230718

ms::libri(
    ms, data.table, qs, glue, cli, optparse
)

set.seed(61787)

for (i in list.files("fn/", full.names = TRUE)) source(i)

# optparse list ----
option_list <- list(
  make_option("--mgi_version",
    type = "character", default = "20220822",
    help = "MGI cohort version in /net/junglebook/magic_data/EHRdata/ [default = %default]"
  ),
  make_option("--ukb_version",
    type = "character", default = "20221117",
    help = "UKB cohort version in /net/junglebook/magic_data/EHRdata/ [default = %default]"
  ),
  make_option("--mgi_weights",
    type = "character", default = "ps_nhw_f,ip_selection_f",
    help = "Name of weight variable to use for MGI [default = %default]"
  )
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

file_paths <- get_files(
  mgi_version = opt$mgi_version,
  ukb_version = opt$ukb_version
)

# data
## mgi
### demo
mgi_demo <- read_qs(file_paths[["mgi"]][["cov_processed_file"]])
setnames(mgi_demo, "DeID_PatientID", "id")

mgi_phecode <- get(load(file_paths[["mgi"]][["phecode_dsb_file"]]))
setnames(mgi_phecode, "IID", "id")
mgi_hyper <- rbindlist(list(
    unique(mgi_phecode[phecode == "401", .(id)])[, .(id, hypertension = 1)],
    data.table(
        id = setdiff(mgi_phecode[, unique(id)],
        mgi_phecode[phecode == "401", unique(id)]),
        hypertension = 0
    )
), use.names = TRUE, fill = TRUE)

### weights
mgi_weights <- read_qs(glue("data/private/mgi/{opt$mgi_version}/weights_{opt$mgi_version}_comb.qs"))
mgi_weight_vars <- c("id", unlist(stringr::str_split(opt$mgi_weights, ",")))

### merge
mgi <- Reduce(
    \(x, y) merge.data.table(x, y, all = FALSE, by = "id"),
    list(
        mgi_demo[, .(id, female, age = round(Enrollment_DaysSinceBirth / 365.25, 1))],
        mgi_hyper[, .(id, hypertension)],
        mgi_weights[, ..mgi_weight_vars]
    )
)

## ukb
### demo
ukb_demo <- read_qs(file_paths[["ukb"]][["demo_file"]])[in_phenome == 1, ]

### pim
ukb_pim0 <- read_qs(file_paths[["ukb"]][["pim0_file"]])

### weights
ukb_weights <- fread(file_paths[["ukb"]][["weight_file"]], colClasses = "character")[
    , .(id = f.eid, weight = as.numeric(LassoWeight))
] |>
    na.omit()

### merge
ukb <- Reduce(
    \(x, y) merge.data.table(x, y, all = FALSE, by = "id"),
    list(
        ukb_demo[, .(id, age = age_at_consent, female = as.numeric(sex == "Female"))],
        ukb_pim0[, .(id, hypertension = ifelse(X401 > 0, 1, 0))],
        ukb_weights
    )
)

# analysis
## function
mod_fn <- function(out, ex, cov = NULL, data, w = NULL) {
    comb <- paste0(c(ex, cov), collapse = " + ")
    if (is.null(cov)) cov <- NA_character_
    mod <- glm(
        formula = as.formula(glue("{out} ~ {comb}")),
        family = quasibinomial(),
        data = data,
        weights = w
    ) |> summary()
    return(data.table(
        outcome    = out,
        exposure   = ex,
        covariates = cov,
        ex_coef    = coef(mod)[ex, 1],
        ex_se      = coef(mod)[ex, 2],
        ex_lower   = coef(mod)[ex, 1] - qnorm(0.975) * coef(mod)[ex, 2],
        ex_upper   = coef(mod)[ex, 1] + qnorm(0.975) * coef(mod)[ex, 2],
        ex_p       = coef(mod)[ex, 4],
        weighted   = !is.null(w)
    ))
}


## unweighted
mgi_um1 <- mod_fn(
    out = "hypertension",
    ex = "female",
    data = mgi
)[, data := "MGI"]
mgi_um2 <- mod_fn(
    out = "hypertension",
    ex = "female",
    cov = "age",
    data = mgi
)[, data := "MGI"]
ukb_um1 <- mod_fn(
    out = "hypertension",
    ex = "female",
    data = ukb
)[, data := "UKB"]
ukb_um2 <- mod_fn(
    out = "hypertension",
    ex = "female",
    cov = "age",
    data = ukb
)[, data := "UKB"]

## weighted
mgi_wm1 <- mod_fn(
    out = "hypertension",
    ex = "female",
    data = mgi,
    w = mgi[[mgi_weight_vars[2]]]
)[, `:=` (data = "MGI", weight = mgi_weight_vars[2])]
mgi_wm2 <- mod_fn(
    out = "hypertension",
    ex = "female",
    cov = "age",
    data = mgi,
    w = mgi[[mgi_weight_vars[2]]]
)[, `:=` (data = "MGI", weight = mgi_weight_vars[2])]
mgi_wm3 <- mod_fn(
    out = "hypertension",
    ex = "female",
    data = mgi,
    w = mgi[[mgi_weight_vars[3]]]
)[, `:=` (data = "MGI", weight = mgi_weight_vars[3])]
mgi_wm4 <- mod_fn(
    out = "hypertension",
    ex = "female",
    cov = "age",
    data = mgi,
    w = mgi[[mgi_weight_vars[3]]]
)[, `:=` (data = "MGI", weight = mgi_weight_vars[3])]
ukb_wm1 <- mod_fn(
    out = "hypertension",
    ex = "female",
    data = ukb,
    w = ukb[, weight]
)[, data := "UKB"]
ukb_wm2 <- mod_fn(
    out = "hypertension",
    ex = "female",
    cov = "age",
    data = ukb,
    w = ukb[, weight]
)[, data := "UKB"]

# combine
combined <- rbindlist(
    list(
        mgi_um1, mgi_um2, ukb_um1, ukb_um2,
        mgi_wm1, mgi_wm2, mgi_wm3, mgi_wm4, ukb_wm1, ukb_wm2
    ),
    use.names = TRUE, fill = TRUE
)

fwrite(
    combined,
    glue("results/targeted_analysis_mgi{opt$mgi_version}_ukb{opt$ukb_version}.csv")
)

cli_alert("Done! 🎉")
