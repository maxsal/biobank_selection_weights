suppressPackageStartupMessages({
  library(data.table)
  library(survey)
  library(cli)
  library(parallel)
  library(future)
  library(parallelly)
  library(doParallel)
  library(foreach)
})

calculate_weighted_prevalences <- function(
    pim_data,
    cov_data,
    pim_id_var   = "IID",
    cov_id_var   = "DeID_PatientID",
    weight_var   = "weights",
    sex_var      = "Sex",
    male_val     = "M",
    female_val   = "F",
    pheinfo      = NULL,
    pheinfo_path = "https://raw.githubusercontent.com/maxsal/public_data/main/phewas/Phecode_Definitions_FullTable_Modified.txt",
    verbose      = TRUE,
    n_cores      = parallelly::availableCores() / 4,
    parallelize  = NULL
) {
    
    # initialize
    if (!is.data.table(pim_data)) pim_data <- as.data.table(pim_data)
    if (!is.data.table(cov_data)) cov_data <- as.data.table(cov_data)
    if (!weight_var %in% names(pim_data)) {
        message(paste0(weight_var), " not found in pim_data. looking in cov_data")
        if (!weight_var %in% names(cov_data)) {
            stop(paste0(weight_var), " not found in cov_data. `weight_var` must be in one of the two data.tables")
        } else {
            keep_vars <- c(cov_id_var, weight_var)
            pim_data <- merge.data.table(pim_data, cov_data[, ..keep_vars], by = cov_id_var, all.x = TRUE)
        }
    }
    
    # load pheinfo
    phe_starts_with_num <- function(x) {
        suppressWarnings({
            !is.na(as.numeric(substr(x[1], 1, 1)))
        })
    }

    # load pheinfo
    if (is.null(pheinfo)) {
        pheinfo <- fread(pheinfo_path, colClasses = "character", showProgress = FALSE)
    }
    xcode <- phe_starts_with_num(pheinfo[, phecode])

    # identify sex specific phecodes
    both   <- pheinfo[sex == "Both", paste0(ifelse(xcode, "X", ""), phecode)]
    male   <- pheinfo[sex == "Male", paste0(ifelse(xcode, "X", ""), phecode)]
    female <- pheinfo[sex == "Female", paste0(ifelse(xcode, "X", ""), phecode)]

    # list ids by sex
    male_ids   <- cov_data[cov_data[[sex_var]] == male_val, ][[cov_id_var]]
    female_ids <- cov_data[cov_data[[sex_var]] == female_val, ][[cov_id_var]]

    # initialize output
    out <- data.table(
        phecode = pheinfo[, paste0("X", phecode)],
        n = NA_real_,
        N = NA_real_
    )

    # subset pim data by sex
    male_pim_data   <- pim_data[which(pim_data[[pim_id_var]] %in% male_ids), ]
    female_pim_data <- pim_data[which(pim_data[[pim_id_var]] %in% female_ids), ]

    # fill out count and length
    if (verbose) cli_alert("phecodes for both sexes...")
    both_design <- svydesign(ids = ~1, data = pim_data, weights = ~get(weight_var))
    if (is.null(parallelize)) {
        both <- both[both %in% names(pim_data)]
        both_out <- list()
        cli_progress_bar(total = length(both), name = "both sexes")
        for (i in seq_along(both)) {
            tmp <- svymean(as.formula(paste0("~`", both[i], "`")), design = both_design, na.rm = TRUE)
            both_out[[i]] <- data.table(
                phecode = both[i],
                prev    = tmp[1],
                se      = SE(tmp)[1]
            )
            cli_progress_update()
        }
        cli_progress_done()
    } else if (parallelize %in% c("parallel", "mclapply")) {
        both_out <- mclapply(
            both,
            \(x) {
                if (!x %in% names(pim_data)) {
                    data.table()
                } else {
                tmp <- svymean(as.formula(paste0("~`", x, "`")), design = both_design, na.rm = TRUE)
                data.table(
                    phecode = x,
                    prev    = tmp[1],
                    se      = SE(tmp)[1]
                )
                }
            },
            mc.cores = n_cores
        )
    } else if (parallelize %in% c("foreach", "doParallel")) {
        cl <- parallel::makeCluster(n_cores, type = "PSOCK")
        doParallel::registerDoParallel(cl)
        both <- both[both %in% names(pim_data)]
        both_out <- foreach(x = both) %dopar% {
            library(survey)
            tmp <- svymean(as.formula(paste0("~`", x, "`")), design = both_design, na.rm = TRUE)
            data.table(
                phecode = x,
                prev    = tmp[1],
                se      = SE(tmp)[1]
            )
        }
    }

    if (verbose) cli_alert("phecodes for males...")
    male_design <- svydesign(ids = ~1, data = male_pim_data, weights = ~get(weight_var))
    if (is.null(parallelize)) {
        male <- male[male %in% names(pim_data)]
        male_out <- list()
        cli_progress_bar(total = length(male), name = "male sex")
        for (i in seq_along(male)) {
            tmp <- svymean(as.formula(paste0("~`", male[i], "`")), design = male_design, na.rm = TRUE)
            male_out[[i]] <- data.table(
                phecode = male[i],
                prev    = tmp[1],
                se      = SE(tmp)[1]
            )
            cli_progress_update()
        }
        cli_progress_done()
    } else if (parallelize %in% c("parallel", "mclapply")) {
        male_out <- mclapply(
            male,
            \(x) {
                if (!x %in% names(pim_data)) {
                    data.table()
                } else {
                    tmp <- svymean(as.formula(paste0("~`", x, "`")), design = male_design, na.rm = TRUE)
                    data.table(
                        phecode = x,
                        prev    = tmp[1],
                        se      = SE(tmp)[1]
                    )
                }
            },
            mc.cores = n_cores
        )
    } else if (parallelize %in% c("foreach", "doParallel")) {
        male <- male[male %in% names(pim_data)]
        male_out <- foreach(x = male) %dopar% {
            library(survey)
            tmp <- svymean(as.formula(paste0("~`", x, "`")), design = male_design, na.rm = TRUE)
            data.table(
                phecode = x,
                prev    = tmp[1],
                se      = SE(tmp)[1]
            )
        }
    }

    if (verbose) cli_alert("phecodes for females...")
    female_design <- svydesign(ids = ~1, data = female_pim_data, weights = ~get(weight_var))
    if (is.null(parallelize)) {
        female <- female[female %in% names(pim_data)]
        female_out <- list()
        cli_progress_bar(total = length(female), name = "female sex")
        for (i in seq_along(female)) {
            tmp <- svymean(as.formula(paste0("~`", female[i], "`")), design = female_design, na.rm = TRUE)
            female_out[[i]] <- data.table(
                phecode = female[i],
                prev    = tmp[1],
                se      = SE(tmp)[1]
            )
            cli_progress_update()
        }
        cli_progress_done()
    } else if (parallelize %in% c("parallel", "mclapply")) {
        female_out <- mclapply(
            female,
            \(x) {
                if (!x %in% names(pim_data)) {
                    data.table()
                } else {
                    tmp <- svymean(as.formula(paste0("~`", x, "`")), design = female_design, na.rm = TRUE)
                    data.table(
                        phecode = x,
                        prev    = tmp[1],
                        se      = SE(tmp)[1]
                    )
                }
            },
            mc.cores = n_cores
        )
    } else if (parallelize %in% c("foreach", "doParallel")) {
        female <- female[female %in% names(pim_data)]
        female_out <- foreach(x = female) %dopar% {
            library(survey)
            tmp <- svymean(as.formula(paste0("~`", x, "`")), design = female_design, na.rm = TRUE)
            data.table(
                phecode = x,
                prev    = tmp[1],
                se      = SE(tmp)[1]
            )
        }
        parallel::stopCluster(cl)
    }

    # calculate prevalence
    out <- rbindlist(list(
        rbindlist(both_out),
        rbindlist(male_out),
        rbindlist(female_out)
    ))
    out[phecode %in% names(pim_data)][]

}
