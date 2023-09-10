suppressPackageStartupMessages({
  library(data.table)
  library(cli)
})

calculate_prevalences <- function(
    pim_data,
    cov_data,
    pim_id_var   = "IID",
    cov_id_var   = "DeID_PatientID",
    sex_var      = "Sex",
    male_val     = "M",
    female_val   = "F",
    pheinfo      = NULL,
    pheinfo_path = "https://raw.githubusercontent.com/maxsal/public_data/main/phewas/Phecode_Definitions_FullTable_Modified.txt",
    progress     = TRUE
) {
    if (!is.data.table(pim_data)) { pim_data <- as.data.table(pim_data) }
    if (!is.data.table(cov_data)) { cov_data <- as.data.table(cov_data) }

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
        phecode = pheinfo[, paste0(ifelse(xcode, "X", ""), phecode)],
        n = NA_real_,
        N = NA_real_
    )

    # subset pim data by sex
    male_pim_data   <- pim_data[which(pim_data[[pim_id_var]] %in% male_ids), ]
    female_pim_data <- pim_data[which(pim_data[[pim_id_var]] %in% female_ids), ]

    # fill out count and length
    if (progress) cli_progress_bar("Both sex phecodes", total = length(both))
    for (i in both) {
        out[phecode == i, `:=` (
            n = sum(pim_data[[i]], na.rm = TRUE),
            N = length(na.omit(pim_data[[i]]))
        )]
        if (progress) cli_progress_update()
    }
    if (progress) cli_progress_done()

    if (progress) cli_progress_bar("Male phecodes", total = length(male))
    for (i in male) {
        out[phecode == i, `:=` (
            n = sum(male_pim_data[[i]], na.rm = TRUE),
            N = length(na.omit(male_pim_data[[i]]))
        )]
        if (progress) cli_progress_update()
    }
    if (progress) cli_progress_done()

    if (progress) cli_progress_bar("Female phecodes", total = length(female))
    for (i in female) {
        out[phecode == i, `:=` (
            n = sum(female_pim_data[[i]], na.rm = TRUE),
            N = length(na.omit(female_pim_data[[i]]))
        )]
        if (progress) cli_progress_update()
    }
    if (progress) cli_progress_done()

    # calculate prevalence
    out <- out[phecode %in% names(pim_data)][]
    out[, prev := n/N][]
}
