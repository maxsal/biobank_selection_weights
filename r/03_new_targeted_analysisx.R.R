# libraries -------------------------------------------------------------------
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
        type = "character", default = "ps_selection,ip_selection",
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

# data ------------------------------------------------------------------------
## mgi
### demo
mgi_demo <- read_qs(glue("data/private/mgi/{opt$mgi_version}/datax_{opt$mgi_version}_comb.qs"))
setnames(mgi_demo, "DeID_PatientID", "id")

mgi_pim <- read_qs(glue("data/private/mgi/{opt$mgi_version}/MGI_PIM0X_{opt$mgi_version}.qs"))

### weights
mgi_weights <- read_qs(glue("data/private/mgi/{opt$mgi_version}/weightsx_{opt$mgi_version}_comb.qs"))
mgi_weight_vars <- c("id", unlist(stringr::str_split(opt$mgi_weights, ",")))

### merge
mgi <- merge_list(
    list(
        mgi_demo[, .(id, female, age = age_at_last_diagnosisx, smoker, cancer = cancerx)],
        mgi_weights[, ..mgi_weight_vars],
        mgi_pim[, .(id, CA_101.41, CA_110.1, CA_108.5, CA_101.71, CA_101.6, CA_101.8, CA_102.1, CV_401, EM_202.2)]
    ),
    verbose = TRUE
)

## ukb
### demo
ukb_demo <- read_qs(glue("data/private/ukb/{opt$ukb_version}/datax_{opt$ukb_version}_comb.qs"))[in_phenome == 1, ]
ukb_demo[, smoker := as.numeric(smk_ev == "Ever")]

### pim
ukb_pim0 <- read_qs(glue("/net/junglebook/home/mmsalva/projects/dissertation/aim_one/data/private/ukb/{opt$ukb_version}/UKB_PIM0X_{opt$ukb_version}.qs"))

### weights
ukb_weights <- fread(file_paths[["ukb"]][["weight_file"]], colClasses = "character")[
    , .(id = f.eid, weight = as.numeric(LassoWeight))
] |>
    na.omit()

### merge
ukb <- merge_list(
    list(
        ukb_demo[, .(id, age = age_at_consent, female = as.numeric(sex == "Female"), smoker, cancer = cancerx)],
        ukb_weights,
        ukb_pim0[, .(id, CA_101.41, CA_110.1, CA_108.5, CA_101.71, CA_101.6, CA_101.8, CA_102.1, CV_401, EM_202.2)]
    ),
    verbose = TRUE
)

# fit models ------------------------------------------------------------------
outcome <- "CA_108.5"
## mgi
### unweighted
mgi_u1 <- glm(as.formula(paste0(outcome, " ~ female")), data = mgi, family = quasibinomial())
mgi_u2 <- glm(as.formula(paste0(outcome, " ~ female + age")), data = mgi, family = quasibinomial())

### ip-weighted
mgi_ip_dsn <- svydesign(
    id = ~id,
    weights = ~ip_selection,
    data = mgi[!is.na(ip_selection), ]
)
mgi_ip1 <- svyglm(as.formula(paste0(outcome, " ~ female")), design = mgi_ip_dsn, family = quasibinomial())
mgi_ip2 <- svyglm(as.formula(paste0(outcome, " ~ female + age")), design = mgi_ip_dsn, family = quasibinomial())

### ps-weighted
mgi_ps_dsn <- svydesign(
    id = ~id,
    weights = ~ps_selection,
    data = mgi[!is.na(ps_selection), ]
)
mgi_ps1 <- svyglm(as.formula(paste0(outcome, " ~ female")), design = mgi_ps_dsn, family = quasibinomial())
mgi_ps2 <- svyglm(as.formula(paste0(outcome, " ~ female + age")), design = mgi_ps_dsn, family = quasibinomial())

## ukb
### unweighted
ukb_u1 <- glm(as.formula(paste0(outcome, " ~ female")), data = ukb, family = quasibinomial())
ukb_u2 <- glm(as.formula(paste0(outcome, " ~ female + age")), data = ukb, family = quasibinomial())

### ip-weighted
ukb_ip_dsn <- svydesign(
    id = ~id,
    weights = ~weight,
    data = ukb[!is.na(weight), ]
)
ukb_ip1 <- svyglm(as.formula(paste0(outcome, " ~ female")), design = ukb_ip_dsn, family = quasibinomial())
ukb_ip2 <- svyglm(as.formula(paste0(outcome, " ~ female + age")), design = ukb_ip_dsn, family = quasibinomial())

# abstract model results ------------------------------------------------------
# helper function
model_results <- function(model_object, outcome = "CA_101.41", exposure = "female", weighted = NA_character_, data = "MGI", r = 3) {
    mod_sum <- summary(model_object)
    data.table(
        outcome    = outcome,
        exposure   = exposure,
        covariates = paste0(setdiff(names(model_object$coefficients), "(Intercept)"), collapse = ", "),
        ex_coef    = coef(mod_sum)[exposure, 1],
        ex_se      = coef(mod_sum)[exposure, 2],
        ex_lower   = coef(mod_sum)[exposure, 1] - qnorm(0.975) * coef(mod_sum)[exposure, 2],
        ex_upper   = coef(mod_sum)[exposure, 1] + qnorm(0.975) * coef(mod_sum)[exposure, 2],
        ex_p       = coef(mod_sum)[exposure, 4],
        weight     = weighted
    )[, `:=`(data = data, ex_print = paste0(
        format(round(ex_coef, r), nsmall = r), " (",
        format(round(ex_lower, r), nsmall = r), ", ",
        format(round(ex_upper, r), nsmall = r), ")"
    ))]
}

summary_table <- rbindlist(list(
    # unweighted
    model_results(mgi_u1, weighted = "unweighted", data = "MGI"),
    model_results(mgi_u2, weighted = "unweighted", data = "MGI"),
    model_results(ukb_u1, weighted = "unweighted", data = "UKB"),
    model_results(ukb_u2, weighted = "unweighted", data = "UKB"),

    # IP-weighted
    model_results(mgi_ip1, weighted = "IP-weighted", data = "MGI"),
    model_results(mgi_ip2, weighted = "IP-weighted", data = "MGI"),
    model_results(ukb_ip1, weighted = "IP-weighted", data = "UKB"),
    model_results(ukb_ip2, weighted = "IP-weighted", data = "UKB"),

    # PS-weighted
    model_results(mgi_ps1, weighted = "PS-weighted", data = "MGI"),
    model_results(mgi_ps2, weighted = "PS-weighted", data = "MGI")
))

# save
fwrite(
    summary_table,
    glue("results/targeted_analysisx_{outcome}_mgi{opt$mgi_version}_ukb{opt$ukb_version}.csv")
)

cli_alert_success("done! 🎉")
