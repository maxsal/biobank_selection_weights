require(data.table)
require(cli)

calculate_prevalences <- function(
    pim_data,
    cov_data,
    pim_id_var   = "IID",
    cov_id_var   = "DeID_PatientID",
    sex_var      = "Sex",
    male_val     = "M",
    female_val   = "F",
    pheinfo_path = "https://gitlab.com/maxsal/public_data/-/raw/main/phewas/Phecode_Definitions_FullTable_Modified.txt",
    progress     = TRUE
) {
    if (!is.data.table(pim_data)) { pim_data <- as.data.table(pim_data) }
    if (!is.data.table(cov_data)) { cov_data <- as.data.table(cov_data) }

    # load pheinfo
    pheinfo <- fread(pheinfo_path, colClasses = "character", showProgress = FALSE)
    # identify sex specific phecodes
    both   <- pheinfo[sex == "Both", paste0("X", phecode)]
    male   <- pheinfo[sex == "Male", paste0("X", phecode)]
    female <- pheinfo[sex == "Female", paste0("X", phecode)]

    # list ids by sex
    male_ids <- cov_data[cov_data[[sex_var]] == male_val, ][[cov_id_var]]
    female_ids <- cov_data[cov_data[[sex_var]] == female_val, ][[cov_id_var]]

    # initialize output
    out <- data.table(
        phecode = pheinfo[, paste0("X", phecode)],
        n = NA_real_,
        N = NA_real_
    )

    # fill out count and length
    if (progress) cli_progress_bar(total = length(both), name = "Both")
    for (i in both) {
        out[phecode == i, `:=` (
            n = sum(pim_data[[i]], na.rm = TRUE),
            N = length(pim_data[[i]])
        )]
        if (progress) cli_progress_update()
    }

    if (progress) cli_progress_bar(total = length(male), name = "Male")
    for (i in male) {
        out[phecode == i, `:=` (
            n = sum(pim_data[pim_data[[pim_id_var]] %in% male_ids, ][[i]], na.rm = TRUE),
            N = length(pim_data[pim_data[[pim_id_var]] %in% male_ids, ][[i]])
        )]
        if (progress) cli_progress_update()
    }

    if (progress) cli_progress_bar(total = length(female), name = "Female")
    for (i in female) {
        out[phecode == i, `:=` (
            n = sum(pim_data[pim_data[[pim_id_var]] %in% female_ids, ][[i]], na.rm = TRUE),
            N = length(pim_data[pim_data[[pim_id_var]] %in% female_ids, ][[i]])
        )]
        if (progress) cli_progress_update()
    }

    # calculate prevalence
    out <- out[phecode %in% names(pim_data)][]
    out[, prev := n/N][]
}
