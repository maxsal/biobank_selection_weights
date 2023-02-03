# estimating weighted and unweighted female-cancer log odds ratio
# author: max salvatore
# date:   20230127

# 1. libraries, functions, and options (outcome agnostic) ----------------------
options(stringsAsFactors = FALSE)

library(data.table)
library(glue)
library(logistf)
library(cli)
library(optparse)
library(survey)

set.seed(61787)

for (i in list.files("fn/", full.names = TRUE)) source(i) # load functions

# optparse list ----------------------------------------------------------------
option_list <- list(
  make_option("--mgi_version", type = "character", default = "20210318",
              help = "Version of MGI data [default = 20210318]")
)

parser <- OptionParser(usage="%prog [options]", option_list = option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

# 2. specifications (specifies outcome) ----------------------------------------
## extract file paths
file_paths <- get_files(mgi_version = opt$mgi_version)

# 3. read data -----------------------------------------------------------------
cli_alert_info("loading data...")

cancer_phecodes <- fread("data/public/cancer_phecodes.txt", colClasses = "character")[[1]]
cancer_vars <- c(glue("X{cancer_phecodes}"))

## mgi
mgi_pim0    <- fread(file_paths[["mgi"]]$pim0_file, select = c("IID", cancer_vars))
mgi_pim0[, cancers := as.numeric(rowSums(.SD, na.rm = TRUE) > 0), .SDcols = cancer_vars]

mgi_demo    <- fread(file_paths[["mgi"]]$demo_file)
mgi_weights <- fread(glue("results/mgi/{opt$mgi_version}/cancer_weights.txt"))

merged <- Reduce(merge, list(
  mgi_pim0[, .(id = IID, cancers)],
  mgi_demo[, .(id = Deid_ID, female = as.numeric(Sex == "F"), age = Age)],
  mgi_weights))

extract_coef <- function(var_name, mod, mod_name = "Unweighted") {
  tmp <- suppressWarnings(suppressMessages(confint(mod)))
  data.table(
    est = coef(mod)[[var_name]],
    lo = tmp[var_name, 1],
    hi = tmp[var_name, 2],
    name = mod_name
  )
}

m1 <- glm(cancers ~ female, data = merged, family = binomial)
m2 <- glm(cancers ~ female, data = merged[!is.na(weights_no_cancer)], family = binomial, weights = weights_no_cancer)
m3 <- glm(cancers ~ female + age, data = merged[!is.na(weights_cancer_indirect)], family = binomial, weights = weights_cancer_indirect)
m4 <- glm(cancers ~ female, data = merged, family = binomial, weights = weights_cancer_direct)

out <- rbindlist(list(
  extract_coef(var_name = "female", mod = m1),
  extract_coef(var_name = "female", mod = m2, mod_name = "No cancer weights"),
  extract_coef(var_name = "female", mod = m3, mod_name = "Cancer weights - indirect"),
  extract_coef(var_name = "female", mod = m4, mod_name = "Cancer weights - direct")
))

# fwrite(out, "bin/female_cancer_logodds_est.csv")
