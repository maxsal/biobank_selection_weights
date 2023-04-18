# read and prepare UKB data
# based almost exclusively on lars' scripts from: https://github.com/umich-cphds/createUKBphenome
# author: max salvatore
# date:   20230418

# libraries --------------------------------------------------------------------
library(data.table)

if (getwd() != "/net/junglebook/home/mmsalva/createUKBphenome") {
  message("Changing directory to '/net/junglebook/home/mmsalva/createUKBphenome'")
  setwd("/net/junglebook/home/mmsalva/createUKBphenome")
}
message(paste0(
  "Make sure your working directory is set to '/net/junglebook/",
  "home/mmsalve/createUKBphenome/', the './data/baskets.txt' and ",
  "withdrawn ('./data/w***_***.txt') files are up-to-date, and ",
  "that you've recently sourced './scripts/function.",
  "summarizeAvailableData.r' and './scripts/function.",
  "createUKBphenome.r'!"
))

source("./scripts/function.reformatUKB.r")
source("./scripts/function.harmonizeICD9.r")
source("/net/junglebook/home/mmsalva/projects/dissertation/aim_one/fn/files-utils.R")

# pull data --------------------------------------------------------------------
## icd9 data
message("Pulling ICD9 data...")
icd9 <- reformatUKB(fields = c(41271, 41281), dataCoding = FALSE, onlyInfo = FALSE)
setnames(icd9, old = c("41271", "41281"), new = c("diagnosis_code", "date"))
icd9[, date := as.Date(date)]

## icd10 data
message("Pulling ICD10 data...")
icd10 <- reformatUKB(fields = c(41270, 41280), dataCoding = FALSE, onlyInfo = FALSE)
setnames(icd10, old = c("41270", "41280"), new = c("diagnosis_code", "date"))
icd10[, date := as.Date(date)]

# function - pull all fields including subfields
# reformatUKB drops fields with multiple subfields :(
pull_ukb_subfields <- function(fields, fieldnames) {
  named_fields        <- fields
  names(named_fields) <- fieldnames
  named_fields        <- sort(named_fields)
  fieldinfo <- all.tables[which(all.tables$"Field ID" %in% named_fields & !is.na(all.tables$basket)), ]
  keep <- which(fieldinfo[, subfield_id] != "n/a")
  if (nrow(fieldinfo) != length(keep)) {
    message(paste0(nrow(fieldinfo) - length(keep), " field(s) are unavailable. continuing with available fields"))
    fieldinfo    <- fieldinfo[keep, ]
    named_fields <- named_fields[which(named_fields %in% fieldinfo[, `Field ID`])]
  }
  allfields <- fieldinfo$subfield_id[which(fieldinfo$subfield_id != "n/a")]

  field_subfields <- strsplit(allfields, ",")
  out <- list()
  for (i in seq_along(named_fields)) {
    out[[i]] <- data.table(
      field    = named_fields[[i]],
      subfield = field_subfields[[i]],
      name     = gsub("_0", "", paste0(names(named_fields)[[i]], "_", seq_len(length(field_subfields[[i]])) - 1))
    )
  }
  out <- do.call(rbind, out)

  selected_columns <- c("id", unique(unlist(field_subfields)))
  message(paste("Reading all entries across", length(selected_columns), "columns"))

  colnumbers <- which(ukb_data_columns %in% selected_columns)

  reformatted <- fread(
    cmd = paste("cut -f", getRanges(colnumbers), data.file), header = T,
    colClasses = "character"
  )

  fieldpattern <- "f\\.(.+)\\..+\\..+"
  entrypattern <- "f\\..+\\.(.+\\..+)"
  keep <- which(rowSums(is.na(reformatted[, -1]) | reformatted[, -1] == "", na.rm = T) != ncol(reformatted[, -1]))
  reformatted <- reformatted[keep, ]

  setnames(
    reformatted,
    out[["subfield"]],
    out[["name"]]
  )

  keep <- which(rowSums(is.na(reformatted[, -1]) | reformatted[, -1] == "", na.rm = T) != ncol(reformatted[, -1]))
  reformatted <- reformatted[keep, ]
  print(paste0(nrow(reformatted), " lines after filtering empty lines"))
  dim(reformatted[keep, ])

  reformatted
}

# pull UKB data by field
dob <- pull_ukb_subfields(
  fields = c(31, 33, 34, 52, 22001, 6141, 21000, 20116, 20117, 21001, 200),
  fieldnames = c("sex", "birth_date", "birth_year", "birth_month", "genetic_sex", "living_with", "ethnicity", "smoker", "drinker", "bmi", "consent_date")
)

dob[, id := as.character(id)]

dob[, `:=`(
  consent_date = as.Date(consent_date),
  dob = as.Date(paste0(birth_year, "-", ifelse(nchar(birth_month) == 2, birth_month, paste0("0", birth_month)), "-15")),
  sex = fcase(
    sex == "0", "Female",
    sex == "1", "Male"
  ),
  genetic_sex = fcase(
    sex == "0", "Female",
    sex == "1", "Male"
  ),
  smoker = fcase(
    smoker == "0", "Never",
    smoker == "1", "Past",
    smoker == "2", "Current"
  )
)]

bmi_cols <- grep("bmi", names(dob), value = TRUE)
dob[, (bmi_cols) := lapply(.SD, as.numeric), .SDcols = bmi_cols]

lw_cols <- grep("living_with", names(dob), value = TRUE)
dob[, (lw_cols) := lapply(.SD, \(x) fcase(
  x == "1", "Husband, wife, or partner",
  x == "2", "Son and/or daughter (include step-children)",
  x == "3", "Brother and/or sister",
  x == "4", "Mother and/or father",
  x == "5", "Grandparent",
  x == "6", "Grandchild",
  x == "7", "Other related",
  x == "8", "Other unrelated",
  x == "-3", "Prefer not to answer"
)), .SDcols = lw_cols]

dob[, `:=`(
  bmi_cat = fcase(
    bmi < 18.5, "Underweight (<18.5)",
    bmi < 25, "Healthy [18.5, 25)",
    bmi < 30, "Overweight [25, 30)",
    bmi >= 30, "Obese [30+)"
  ),
  ethnicity = fcase(
    ethnicity == "-1", "Do not know",
    ethnicity == "-3", "Prefer not to answer",
    ethnicity == "5", "Chinese",
    ethnicity == "6", "Other ethnic group",
    stringr::str_ends(ethnicity, "1"), "White",
    stringr::str_ends(ethnicity, "2"), "Mixed",
    stringr::str_ends(ethnicity, "3"), "Asian",
    stringr::str_ends(ethnicity, "4"), "Black",
    ethnicity == "", NA_character_,
    default = NA_character_
  ),
  age_at_consent = round(as.numeric(consent_date - dob) / 365.25, 1)
)][, `:=` (
  age_cat = cut(age_at_consent, seq(0, 85, 10), right = FALSE)
)]


# process data -----------------------------------------------------------------
message("Calculating DSB for ICD9 data")
icd9 <- merge(icd9, dob[, .(id, dob)], by = "id", all.x = TRUE)
icd9[, dsb := round(as.numeric(((date - dob))))]
icd9[, c("date", "dob") := NULL][, lexicon := "icd9"]
print(paste(nrow(icd9), "lines of ICD9 codes for", length(unique(icd9[, id])), "people"))

message("Calculating DSB for ICD10 data")
icd10 <- merge(icd10, dob[, .(id, dob)], by = "id", all.x = TRUE)
icd10[, dsb := round(as.numeric(((date - dob))))]
icd10[, c("date", "dob") := NULL][, lexicon := "icd10"]
print(paste(nrow(icd10), "lines of ICD10 codes for", length(unique(icd10[, id])), "people"))

## stack icd9 and icd10 data
ukb_icd_dsb <- rbindlist(list(
  icd9, icd10
))

## remove unncessary quotes
ukb_icd_dsb[, diagnosis_code := gsub("\"", "", diagnosis_code)]
print(paste(nrow(ukb_icd_dsb), "lines of ICD codes for", length(unique(ukb_icd_dsb[, id])), "people"))

## save
save_stamp <- format(as.Date(file.info("./data/merged_baskets.txt")$ctime), "%Y%m%d")
fwrite(
  x = ukb_icd_dsb,
  file = paste0("./results/UKB_PHENOME_ICD_DSB_", save_stamp, ".txt")
)
save_qs(
  x = ukb_icd_dsb,
  file = paste0("./results/UKB_PHENOME_ICD_DSB_", save_stamp, ".qs")
)
message(paste0("File saved as './results/UKB_PHENOME_ICD_DSB_", save_stamp, ".txt'!"))

# load necessary information for phecode mapping -------------------------------
# Phecode information
load(file = "./PheWAS/data/pheinfo.rda")
pheinfo <- data.table(pheinfo)
pheinfo <- pheinfo[, .(phecode, description, groupnum, group, color)]

# Gender rules
load(file = "./PheWAS/data/gender_restriction.rda")
gender_restriction <- data.table(gender_restriction)
gender_restriction <- rbind(
  data.table("phecode" = gender_restriction[male_only == FALSE & female_only == FALSE, phecode], "sex" = "Both"),
  data.table("phecode" = gender_restriction[male_only == TRUE, phecode], "sex" = "Male"),
  data.table("phecode" = gender_restriction[female_only == TRUE, phecode], "sex" = "Female")
)

# Exclusion / roll up rules
pheinfo2 <- unique(fread("./data/phecode_icd9_rolled.csv",
  colClasses = "character",
  select = c("PheCode", "Excl. Phecodes", "Excl. Phenotypes", "Rollup", "Leaf"),
  col.names = c("phecode", "phecode_exclude_range", "phecode_exclude_phenotypes", "rollup", "leaf")
))

# create one data.table with all criteria
pheinfo <- merge(pheinfo, pheinfo2, by = "phecode")
pheinfo <- merge(pheinfo, gender_restriction, by = "phecode")
pheinfo <- pheinfo[, c("phecode", "description", "sex", "rollup", "leaf", "groupnum", "group", "color", "phecode_exclude_range", "phecode_exclude_phenotypes")]

# add phecode information that's missing (collected form earlier versions)
pheinfoOLD <- fread("./data/Phecode_Definitions_FullTable_Modified.txt", colClasses = "character")
pheinfo <- rbind(pheinfo, pheinfoOLD[!phecode %in% pheinfo$phecode, ])

# Phecode that should not be rolled up
norollup <- pheinfo$phecode[which(pheinfo$rollup == 0)]
# add manual no rollup rules:
norollup <- c(norollup, fread("./data/no_rollup_extra.txt", colClasses = "character")$phecode)

# map ICD9 codes to phecodes ---------------------------------------------------
message("Mapping ICD9 codes to phecodes")

# read PheWAS map (downloaded from https://phewascatalog.org/phecodes; selected "export all" top right corner)
icd9map <- fread("./data/phecode_icd9_rolled.csv", colClasses = "character")

ICD9codes <- fread("./data/coding87.tsv")
ICD9codes[!grepl("Block", coding), ICD9 := sapply(coding, harmonizeICD9)]
codeICD9 <- ICD9codes[, sort(unique(ICD9))]

mappedICD9Codes <- NULL
icd9map_new <- list()
pb <- txtProgressBar(max = nrow(icd9map), width = 50, style = 3)
for (i in seq_len(nrow(icd9map))) {
  mapped9 <- grep(icd9map$ICD9[i], codeICD9)
  if (length(mapped9) > 0) {
    mappedICD9Codes <- unique(c(mappedICD9Codes, codeICD9[mapped9]))
    icd9map_new[[icd9map$ICD9[i]]] <- data.table(
      "phecode" = icd9map$PheCode[i],
      "ICD9"    = codeICD9[mapped9]
    )
  }
  setTxtProgressBar(pb, i)
}
close(pb)
icd9key <- rbindlist(icd9map_new)
icd9unmapped <- codeICD9[!codeICD9 %in% mappedICD9Codes]

# roll up PheWAS codes
icd9key$added <- "original"
pcodes <- unique(c(icd9key$phecode, gsub("\\..+", "", icd9key$phecode), gsub("(\\..).", "\\1", icd9key$phecode)))
pcodes <- sort(pcodes)

for (p in seq_along(pcodes)) {
  if (grepl("\\.", pcodes[p])) {
    iSub <- which(icd9key$phecode %in% pcodes[which(grepl(paste0("^", pcodes[p]), pcodes) &
      nchar(pcodes) > nchar(pcodes[p]))] &
      !icd9key$phecode %in% norollup)
  } else {
    iSub <- which(icd9key$phecode %in% pcodes[grep(paste0("^", pcodes[p], "\\."), pcodes)] &
      !icd9key$phecode %in% norollup)
  }
  if (length(iSub) == 0) next
  iTop <- icd9key$ICD9[which(icd9key$phecode == pcodes[p])]
  addTop <- which(icd9key$ICD9 %in% unique(icd9key$ICD9[iSub]) & !icd9key$ICD9 %in% iTop)
  if (length(addTop) == 0) next

  addKey <- icd9key[addTop, ]
  addKey$phecode <- pcodes[p]
  addKey$added <- "rolled up PheWAS code"
  icd9key <- rbind(icd9key, addKey)
}

## Add ICD code description and phecode description
icd9key <- unique(icd9key)
print("Rollup of phewas codes (ICD9 code)")
print(table(icd9key$added))

icd9key <- merge(icd9key, ICD9codes[, .(ICD9, meaning, node_id, parent_id, selectable)], by = "ICD9")
icd9key <- merge(icd9key, pheinfo, by = "phecode")
icd9key <- icd9key[, c(
  "ICD9", "meaning", "node_id", "parent_id", "selectable", "phecode", "description",
  "group", "groupnum", "added", "sex", "rollup", "leaf"
)]

# remove any issues you might have found
icd9map_remove <- fread("./data/remove_icd9_map.txt", colClasses = "character")
icd9key <- icd9key[!paste(ICD9, phecode, sep = ":") %in% icd9map_remove[, paste(ICD9, phecode, sep = ":")], ]

## map ICD9 data
icd9 <- ukb_icd_dsb[lexicon == "icd9"]
icd9[, diagnosis_code := sapply(diagnosis_code, harmonizeICD9)]
icd9 <- unique(icd9)

phecode1 <- merge(icd9, icd9key[, .(ICD9, phecode)],
  by.x = "diagnosis_code", by.y = "ICD9",
  allow.cartesian = TRUE
)

# map ICD10 codes to phecodes --------------------------------------------------
message("Mapping ICD10 codes to phecodes")

# read PheWAS map (downloaded from https://phewascatalog.org/phecodes_icd10; selected "export all" top right corner)
icd10map <- fread("./data/phecode_icd10.csv", colClasses = "character")

# read UKB coding
ICD10codes <- fread("./data/coding19.tsv")
ICD10codes <- ICD10codes[!grepl("Block", coding), ]
ICD10codes[, ICD10category := gsub("([A-Z][0-9]{2}).+", "\\1", coding)]
ICD10codes[, ICD10suffix := gsub("[A-Z].+", "", gsub("^[A-Z][0-9]{2}", "", coding))]
ICD10codes[, ICD10 := paste0(ICD10category, ifelse(ICD10suffix == "", "", "."), ICD10suffix)]
ICD10codes <- ICD10codes[, c("ICD10category", "ICD10suffix") := NULL]

codeICD10 <- ICD10codes[, sort(unique(ICD10))]
mappedICD10Codes <- NULL
icd10map_new <- list()
pb <- txtProgressBar(max = nrow(icd10map), width = 50, style = 3)
for (i in seq_len(nrow(icd10map))) {
  mapped10 <- grep(icd10map$ICD10[i], codeICD10)
  if (length(mapped10) > 0) {
    mappedICD10Codes <- unique(c(mappedICD10Codes, codeICD10[mapped10]))
    icd10map_new[[icd10map$ICD10[i]]] <- data.table(
      "phecode" = icd10map$PheCode[i],
      "ICD10" = codeICD10[mapped10]
    )
  }
  setTxtProgressBar(pb, i)
}
close(pb)
icd10key <- rbindlist(icd10map_new)
icd10unmapped <- codeICD10[!codeICD10 %in% mappedICD10Codes]

#### roll up phewas codes
icd10key$added <- "original"
pcodes <- unique(c(icd10key$phecode, gsub("\\..+", "", icd10key$phecode), gsub("(\\..).", "\\1", icd10key$phecode)))
pcodes <- sort(pcodes)
for (p in seq_along(pcodes)) {
  if (grepl("\\.", pcodes[p])) {
    iSub <- which(icd10key$phecode %in% pcodes[which(grepl(paste0("^", pcodes[p]), pcodes) &
      nchar(pcodes) > nchar(pcodes[p]))] &
      !icd10key$phecode %in% norollup)
  } else {
    iSub <- which(icd10key$phecode %in% pcodes[grep(paste0("^", pcodes[p], "\\."), pcodes)] &
      !icd10key$phecode %in% norollup)
  }
  if (length(iSub) == 0) next
  iTop <- icd10key$ICD10[which(icd10key$phecode == pcodes[p])]
  addTop <- which(icd10key$ICD10 %in% unique(icd10key$ICD10[iSub]) & !icd10key$ICD10 %in% iTop)
  if (length(addTop) == 0) next

  addKey <- icd10key[addTop, ]
  addKey$phecode <- pcodes[p]
  addKey$added <- "rolled up PheWAS code"
  icd10key <- rbind(icd10key, addKey)
}

icd10key <- unique(icd10key)

print("Rollup of phewas codes (ICD10 code)")
print(table(icd10key$added))

icd10key <- merge(icd10key, ICD10codes[, .(ICD10, meaning, node_id, parent_id, selectable)], by = "ICD10")
icd10key <- merge(icd10key, pheinfo, by = "phecode")
icd10key <- icd10key[, c(
  "ICD10", "meaning", "node_id", "parent_id", "selectable", "phecode", "description",
  "group", "groupnum", "added", "sex", "rollup", "leaf"
)]


## remove any issues you might have found
icd10map_remove <- fread("./data/remove_icd10_map.txt", colClasses = "character")
icd10key <- icd10key[!paste(ICD10, phecode, sep = ":") %in% icd10map_remove[, paste(ICD10, phecode, sep = ":")], ]

## map ICD10 data
icd10 <- ukb_icd_dsb[lexicon == "icd10", ]
icd10[, ICD10category := gsub("([A-Z][0-9]{2}).+", "\\1", diagnosis_code)]
icd10[, ICD10suffix := gsub("[A-Z].+", "", gsub("^[A-Z][0-9]{2}", "", diagnosis_code))]
icd10[, diagnosis_code := paste0(ICD10category, ifelse(ICD10suffix == "", "", "."), ICD10suffix)]
icd10[, c("ICD10category", "ICD10suffix") := NULL]
icd10 <- unique(icd10)

phecode2 <- merge(icd10, icd10key[, .(ICD10, phecode)],
  by.x = "diagnosis_code", by.y = "ICD10",
  allow.cartesian = TRUE
)

## stack icd-phecode mapped dsb data
ukb_phecode <- unique(rbindlist(list(
  phecode1,
  phecode2
)))

# save mapped data -------------------------------------------------------------
fwrite(
  x = ukb_phecode,
  file = paste0("./results/UKB_PHECODE_DSB_MAPPED_", save_stamp, ".txt")
)
save_qs(
  x = ukb_phecode,
  file = paste0("./results/UKB_PHECODE_DSB_MAPPED_", save_stamp, ".qs")
)
message(paste0("File saved as './results/UKB_PHECODE_DSB_MAPPED_", save_stamp, ".txt'!"))

# save sex data ----------------------------------------------------------------
first_phe <- ukb_phecode[ukb_phecode[, .I[which.min(dsb)], id]$V1]
last_phe <- ukb_phecode[ukb_phecode[, .I[which.max(dsb)], id]$V1]
ehr_followup <- merge.data.table(
  first_phe[, .(id, first_dsb = dsb, age_at_first_diagnosis = round(dsb / 365.25, 1))],
  last_phe[, .(id, last_dsb = dsb, age_at_last_diagnosis = round(dsb / 365.25, 1))],
  by = "id"
)[, ehr_days := as.numeric(last_dsb - first_dsb)][, ehr_years := round(ehr_days / 365.25, 1)]
dob <- merge.data.table(
  dob,
  ehr_followup,
  by = "id",
  all.x = TRUE
)
drop <- c("birth_year", "birth_month", lw_cols)

## add comorbidity indicators
cancer_phecodes <- fread("/net/junglebook/home/mmsalva/projects/dissertation/aim_one/data/public/cancer_phecodes.txt",
                         colClasses = "character"
)[[1]]

comorbid <- list(
  "cad"                = list("phecodes" = c("411.4")),
  "diabetes"           = list("phecodes" = c("250")),
  "hypertension"       = list("phecodes" = c("272.12")),
  "mixed_hypertension" = list("phecodes" = c("272.13")),
  "vitamin_d"          = list("phecodes" = c("261.4")),
  "depression"         = list("phecodes" = c("296.2")),
  "anxiety"            = list("phecodes" = c("300")),
  "bipolar"            = list("phecodes" = c("296.1")),
  "cancer"             = list("phecodes" = cancer_phecodes)
)

message("identifying cases and creating indicator variables...")
for (i in names(comorbid)) {
  comorbid[[i]][["ids"]] <- ukb_phecode[phecode %in% comorbid[[i]][["phecodes"]], id] |>
    unique()
}

for (i in names(comorbid)) {
  set(dob, j = i, value = fifelse(dob[["id"]] %in% comorbid[[i]][["ids"]], 1, 0))
}

dob[, triglycerides := fifelse(hypertension == 0 & mixed_hypertension == 0, 0, 1)]


## save
fwrite(
  x    = dob[, !..drop],
  file = paste0("./results/UKB_SEX_", save_stamp, ".txt")
)
message(paste0("File save as './results/UKB_SEX_", save_stamp, ".txt'!"))
save_qs(
  x    = dob[, !..drop],
  file = paste0("./results/UKB_SEX_", save_stamp, ".qs")
)
