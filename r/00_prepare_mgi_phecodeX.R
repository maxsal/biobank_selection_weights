# libraries --------------------------------------------------------------------
ms::libri(
    ms, qs, data.table, tidyverse, PheWAS/PheWAS, glue, optparse
)

# optparse list ---
option_list <- list(
    make_option("--mgi_version",
        type = "character", default = "20220822",
        help = "MGI cohort version [default = '20220822']"
    )
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args   <- parse_args(parser, positional_arguments = 0)
opt    <- args$options
print(opt)

# load data --------------------------------------------------------------------
## internal
mgi_icd  <- qread(glue("data/private/mgi/{opt$mgi_version}/MGI_ICD_{opt$mgi_version}.qs"))
mgi_demo <- qread(glue("data/private/mgi/{opt$mgi_version}/data_{opt$mgi_version}_comb.qs"))
setnames(mgi_icd, c("IID", "DiagnosisCode", "lexicon", "DaysSinceBirth"), c("id", "code", "vocabulary_id", "dsb"))

## external
### mapping table
phecodex_map <- fread("https://raw.githubusercontent.com/PheWAS/PhecodeX/main/phecodeX_R_map.csv",
                      colClasses = "character")
### rollup map
phecodex_rollup <- fread("https://raw.githubusercontent.com/PheWAS/PhecodeX/main/phecodeX_R_rollup_map.csv",
                         colClasses = "character")
### sex restriction map
phecodex_sex <- fread("https://raw.githubusercontent.com/PheWAS/PhecodeX/main/phecodeX_R_sex.csv")
### phecodex info
phecodex_info <- fread("https://raw.githubusercontent.com/PheWAS/PhecodeX/main/phecodeX_R_labels.csv")

# map phenotypes ---------------------------------------------------------------
mapped <- mapCodesToPhecodes(
    input = mgi_icd, vocabulary.map = phecodex_map,
    rollup.map = phecodex_rollup, make.distinct = TRUE
)

id_sex <- mgi_demo[, .(id = DeID_PatientID, sex = Sex)]

# remove sex-specific phecodes that are discordant with individual's reported sex at birth
restrictPhecodesBySex_mod <- function(phenotypes, id.sex, by_var = "person_id", sex_var = "sex") {
    data <- merge.data.table(
        as.data.table(phenotypes),
        as.data.table(id.sex),
        by = by_var, all.x = TRUE
    )
    # Get the restrictions found in the phenotypes data frame
    current_sex_restriction <- phecodex_sex[phecodex_sex$phecode %in% unique(data[, phecode]), ] |>
        as.data.table()
    # Get male and female-only phenotypes
    male_only <- current_sex_restriction[current_sex_restriction$male_only, phecode]
    female_only <- current_sex_restriction[current_sex_restriction$female_only, phecode]
    # Set row column matches to NA where inds of a sex meet restricted phenotypes
    data[phecode %in% male_only & sex == "F", phecode := NA]
    data[phecode %in% female_only & sex == "M", phecode := NA]

    na.omit(data)[, (sex_var) := NULL][]
}

## remove sex person-phenotype discordant pairs
restricted <- restrictPhecodesBySex_mod(mapped, id.sex = id_sex, by_var = "id")

## save
qsave(
    restricted,
    file = glue("data/private/mgi/{opt$mgi_version}/MGI_FULL_PHECODEX_DSB_{opt$mgi_version}.qs")
)

# get first phecode per person -------------------------------------------------
first_restricted <- restricted[ restricted[, .I[which.min(dsb)], by = c("id", "phecode")][["V1"]], ]

## save
qsave(
    first_restricted,
    file = glue("data/private/mgi/{opt$mgi_version}/MGI_FIRST_PHECODEX_DSB_{opt$mgi_version}.qs")
)

# get phecode info -------------------------------------------------------------
phecodex_n <- first_restricted[, .N, phecode][order(-N)]
phecodex_n <- merge.data.table(
    phecodex_n, phecodex_info, by.x = "phecode", by.y = "phenotype",
    all.x = TRUE
)[order(-N), ]
phecodex_n[, prev := N / length(unique(first_restricted$id))]

## save
qsave(
    phecodex_n,
    file = glue("data/private/mgi/{opt$mgi_version}/MGI_PHECODEX_N_{opt$mgi_version}.qs")
)

cli_alert_success("Finished! 🎉")
