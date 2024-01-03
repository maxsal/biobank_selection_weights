# Conduct targeted hypertension ~ sex analysis in MGI and UKB
# author:  max salvatore
# date:    20230718

ms::libri(
    ms, data.table, qs, glue, cli, optparse, survey
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

mgi_phecode <- read_qs(glue("data/private/mgi/{opt$mgi_version}/MGI_FULL_PHECODE_DSB_{opt$mgi_version}.qs"))
setnames(mgi_phecode, "IID", "id")

### weights
mgi_weights <- read_qs(glue("data/private/mgi/{opt$mgi_version}/weights_{opt$mgi_version}_comb.qs"))
mgi_weight_vars <- c("id", unlist(stringr::str_split(opt$mgi_weights, ",")))

### merge
mgi <- merge_list(
        list(
            mgi_demo[, .(id, female, age = round(Enrollment_DaysSinceBirth / 365.25, 1), smoker, cancer)],
            mgi_weights[, ..mgi_weight_vars]
        ),
        verbose = TRUE
)

## ukb
### demo
ukb_demo <- read_qs(file_paths[["ukb"]][["demo_file"]])[in_phenome == 1, ]
ukb_demo[, smoker := as.numeric(smk_ev == "Ever")]

### pim
ukb_pim0 <- read_qs(file_paths[["ukb"]][["pim0_file"]])

### weights
ukb_weights <- fread(file_paths[["ukb"]][["weight_file"]], colClasses = "character")[
    , .(id = f.eid, weight = as.numeric(LassoWeight))
] |>
    na.omit()

### merge
ukb <- merge_list(
    list(
        ukb_demo[, .(id, age = age_at_consent, female = as.numeric(sex == "Female"), smoker, cancer)],
        ukb_weights
    ),
    verbose = TRUE
)


# analysis
## function
mod_fn <- function(out, ex, cov = NULL, data, w = NULL, verbose = FALSE) {
    comb <- paste0(c(ex, cov), collapse = " + ")
    if (is.null(cov)) cov <- NA_character_
    if (sum(is.na(w)) > 0) {
        if (verbose) cli_alert("Missing weights! subsetting... 🚨")
        data <- data[!is.na(w)]
        w    <- w[!is.na(w)]
    }
    svydsn <- svydesign(
        id      = ~1,
        weights = w,
        data    = data
    )
    mod <- svyglm(
        formula = as.formula(glue("{out} ~ {comb}")),
        family = quasibinomial(),
        design = svydsn,
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

# ukb <- merge.data.table(
#     ukb, ukb_pim0[, .(id, hypertension = ifelse(X401 > 0, 1, 0), crc = as.numeric(X153 > 0))]
# )

# ukb_mods <- targeted_analysis(
#     out = c("hypertension", "crc"),
#     ex = "female",
#     cov = "age",
#     data = ukb
# )

results <- list()
outcomes <- c("X153")
cli_progress_bar(total = length(outcomes), format = "estimating {outcome} {pb_bar} {pb_percent} | ETA: {pb_eta}")
for (i in seq_along(outcomes)) {
        outcome <- outcomes[i]
        cli_progress_update()
        tmp_outcome <- c("id", outcome)
        mgi_out_ids <- mgi_phecode[phecode == gsub("X", "", outcome), unique(id)]
        mgi_no_out_ids <- setdiff(mgi_phecode[, unique(id)], mgi_out_ids)
        mgi_out_tmp <- data.table(
            id      = c(mgi_out_ids, mgi_no_out_ids),
            outcome = c(rep(1, length(mgi_out_ids)), rep(0, length(mgi_no_out_ids)))
        )
        setnames(mgi_out_tmp, "outcome", outcome)
        mgi_tmp <- merge.data.table(
            mgi,
            mgi_out_tmp,
            by = "id"
        )
        mgi_tmp <- mgi_tmp[!is.na(get(outcome))]

        ukb_tmp_out <- c("id", outcome)
        ukb_pim_tmp <- ukb_pim0[, ..ukb_tmp_out]
        ukb_pim_tmp[[outcomes[i]]] <- as.numeric(ukb_pim_tmp[[outcomes[i]]] > 0)
        ukb_tmp <- merge.data.table(
            ukb,
            ukb_pim_tmp,
            by = "id"
        )
        ukb_tmp <- ukb_tmp[!is.na(get(outcomes[i]))]

        mgi_um1 <- mod_fn(
            out = outcomes[i],
            ex = "female",
            data = mgi_tmp
        )[, data := "MGI"]
        mgi_um2 <- mod_fn(
            out = outcomes[i],
            ex = "female",
            cov = "age",
            data = mgi_tmp
        )[, data := "MGI"]
        ukb_um1 <- mod_fn(
            out = outcomes[i],
            ex = "female",
            data = ukb_tmp
        )[, data := "UKB"]
        ukb_um2 <- mod_fn(
            out = outcomes[i],
            ex = "female",
            cov = "age",
            data = ukb_tmp
        )[, data := "UKB"]

        ## weighted
        mgi_wm1 <- mod_fn(
            out = outcomes[i],
            ex = "female",
            data = mgi_tmp,
            w = mgi_tmp[[mgi_weight_vars[2]]]
        )[, `:=`(data = "MGI", weight = mgi_weight_vars[2])]
        mgi_wm2 <- mod_fn(
            out = outcomes[i],
            ex = "female",
            cov = "age",
            data = mgi_tmp,
            w = mgi_tmp[[mgi_weight_vars[2]]]
        )[, `:=`(data = "MGI", weight = mgi_weight_vars[2])]
        mgi_wm3 <- mod_fn(
            out = outcomes[i],
            ex = "female",
            data = mgi_tmp,
            w = mgi_tmp[[mgi_weight_vars[3]]]
        )[, `:=`(data = "MGI", weight = mgi_weight_vars[3])]
        mgi_wm4 <- mod_fn(
            out = outcomes[i],
            ex = "female",
            cov = "age",
            data = mgi_tmp,
            w = mgi_tmp[[mgi_weight_vars[3]]]
        )[, `:=`(data = "MGI", weight = mgi_weight_vars[3])]
        ukb_wm1 <- mod_fn(
            out = outcomes[i],
            ex = "female",
            data = ukb_tmp,
            w = ukb_tmp[, weight]
        )[, data := "UKB"]
        ukb_wm2 <- mod_fn(
            out = outcomes[i],
            ex = "female",
            cov = "age",
            data = ukb_tmp,
            w = ukb_tmp[, weight]
        )[, data := "UKB"]

    mgi_results <- rbindlist(list(
        mgi_um1, mgi_um2, mgi_wm1, mgi_wm2, mgi_wm3, mgi_wm4
    ), use.names = TRUE, fill = TRUE)
    ukb_results <- rbindlist(list(
        ukb_um1, ukb_um2, ukb_wm1, ukb_wm2
    ), use.names = TRUE, fill = TRUE)

    results[[i]] <- rbindlist(list(
        mgi_results, ukb_results
    ), use.names = TRUE, fill = TRUE)
}
cli_progress_done()

# results <- list()
# outcomes <- c("X165", "X165.1", "X411", "X150", "X157", "X250.2", "X362.2")
# cli_progress_bar("Running models...", total = length(outcomes))
# for (i in seq_along(outcomes)) {
#     cli_progress_update()
#     tmp_outcome <- c("id", outcomes[i])
#     mgi_out_ids <- mgi_phecode[phecode == gsub("X", "", outcomes[i]), unique(id)]
#     mgi_no_out_ids <- setdiff(mgi_phecode[, unique(id)], mgi_out_ids)
#     mgi_out_tmp <- data.table(
#         id      = c(mgi_out_ids, mgi_no_out_ids),
#         outcome = c(rep(1, length(mgi_out_ids)), rep(0, length(mgi_no_out_ids)))
#     )
#     setnames(mgi_out_tmp, "outcome", outcomes[i])
#     mgi_tmp <- merge.data.table(
#         mgi,
#         mgi_out_tmp,
#         by = "id"
#     )
#     mgi_tmp <- mgi_tmp[!is.na(get(outcomes[i]))]

#     ukb_tmp_out <- c("id", outcomes[i])
#     ukb_pim_tmp <- ukb_pim0[, ..ukb_tmp_out]
#     ukb_pim_tmp[[outcomes[i]]] <- as.numeric(ukb_pim_tmp[[outcomes[i]]] > 0)
#     ukb_tmp <- merge.data.table(
#         ukb,
#         ukb_pim_tmp,
#         by = "id"
#     )
#     ukb_tmp <- ukb_tmp[!is.na(get(outcomes[i]))]

#     mgi_um1 <- mod_fn(
#         out = outcomes[i],
#         ex = "smoker",
#         data = mgi_tmp
#     )[, data := "MGI"]
#     mgi_um2 <- mod_fn(
#         out = outcomes[i],
#         ex = "smoker",
#         cov = "age",
#         data = mgi_tmp
#     )[, data := "MGI"]
#     ukb_um1 <- mod_fn(
#         out = outcomes[i],
#         ex = "smoker",
#         data = ukb_tmp
#     )[, data := "UKB"]
#     ukb_um2 <- mod_fn(
#         out = outcomes[i],
#         ex = "smoker",
#         cov = "age",
#         data = ukb_tmp
#     )[, data := "UKB"]

#     ## weighted
#     mgi_wm1 <- mod_fn(
#         out = outcomes[i],
#         ex = "smoker",
#         data = mgi_tmp,
#         w = mgi_tmp[[mgi_weight_vars[2]]]
#     )[, `:=`(data = "MGI", weight = mgi_weight_vars[2])]
#     mgi_wm2 <- mod_fn(
#         out = outcomes[i],
#         ex = "smoker",
#         cov = "age",
#         data = mgi_tmp,
#         w = mgi_tmp[[mgi_weight_vars[2]]]
#     )[, `:=`(data = "MGI", weight = mgi_weight_vars[2])]
#     mgi_wm3 <- mod_fn(
#         out = outcomes[i],
#         ex = "smoker",
#         data = mgi_tmp,
#         w = mgi_tmp[[mgi_weight_vars[3]]]
#     )[, `:=`(data = "MGI", weight = mgi_weight_vars[3])]
#     mgi_wm4 <- mod_fn(
#         out = outcomes[i],
#         ex = "smoker",
#         cov = "age",
#         data = mgi_tmp,
#         w = mgi_tmp[[mgi_weight_vars[3]]]
#     )[, `:=`(data = "MGI", weight = mgi_weight_vars[3])]
#     ukb_wm1 <- mod_fn(
#         out = outcomes[i],
#         ex = "smoker",
#         data = ukb_tmp,
#         w = ukb_tmp[, weight]
#     )[, data := "UKB"]
#     ukb_wm2 <- mod_fn(
#         out = outcomes[i],
#         ex = "smoker",
#         cov = "age",
#         data = ukb_tmp,
#         w = ukb_tmp[, weight]
#     )[, data := "UKB"]

#     mgi_results <- rbindlist(list(
#         mgi_um1, mgi_um2, mgi_wm1, mgi_wm2, mgi_wm3, mgi_wm4
#     ), use.names = TRUE, fill = TRUE)
#     ukb_results <- rbindlist(list(
#         ukb_um1, ukb_um2, ukb_wm1, ukb_wm2
#     ), use.names = TRUE, fill = TRUE)

#     results[[i]] <- rbindlist(list(
#         mgi_results, ukb_results
#     ), use.names = TRUE, fill = TRUE)
# }
# cli_progress_done()

# combine
combined <- rbindlist(
    results,
    use.names = TRUE, fill = TRUE
)

fwrite(
    combined,
    glue("results/targeted_analysis_mgi{opt$mgi_version}_ukb{opt$ukb_version}.csv")
)

cli_alert("Done! 🎉")
