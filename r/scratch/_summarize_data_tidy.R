# libraries, functions, and options --------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(MatchIt)
  library(glue)
  library(qs)
  library(optparse)
  library(tidyverse)
  library(survey)
})

set.seed(61787)

for (i in list.files("fn/", full.names = TRUE)) source(i)

# optparse list ----------------------------------------------------------------
option_list <- list(
  make_option("--mgi_version",
    type = "character", default = "20220822",
    help = "Version of MGI data [default = %default]"
  ),
  make_option("--mgi_cohort",
    type = "character", default = "comb",
    help = "Cohort of MGI used in weighting (comb, bb, mend, mhb) [default = %default]"
  ),
  make_option("--ukb_version",
    type = "character", default = "20221117",
    help = "Version of UKB data [default = %default]"
  )
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

## extract file paths
file_paths <- get_files(
  mgi_version = opt$mgi_version,
  ukb_version = opt$ukb_version
)

# read data --------------------------------------------------------------------
## ukb
message("loading ukb data...")
### demographics
ukb_demo <- fread(file_paths[["ukb"]][["demo_file"]],
  na.strings = c("", "NA", "."),
  colClass = "character"
)
ukb_demo <- ukb_demo[, .(
  id   = as.character(id),
  dob  = as.Date(dob),
  age  = as.numeric(age_at_consent),
  ethn = ethnicity,
  sex,
  smoker
)]
ukb_demo <- ukb_demo[complete.cases(ukb_demo), ]

weights <- read_delim(
  "/net/junglebook/home/mmsalva/createUKBphenome/data/UKBSelectionWeights.tab",
  col_types = "c",
  show_col_types = FALSE
) |>
  rename(
    id     = "f.eid",
    weight = "LassoWeight"
  )

ukb_demo <- as_tibble(ukb_demo)

merged <- inner_join(
  ukb_demo,
  weights,
  by = "id"
)

dsn <- svydesign(id = ~1, weights = ~weight, data = merged)

svymean(~age, dsn)
svymean(~ as.numeric(sex == "Female"), dsn)
tmp <- svymean(~ethn, dsn)
svymean(~smoker, dsn)

tmp |> as_tibble(rownames = "rowname")

# summarize character and factor variables
ch_sum <- function(d, v, dsn = NULL) {
  if (!is.character(d[[v]]) & !is.factor(d[[v]])) {
    stop(paste0("ch_sum function only works with character or factor variables. offending variable: ", v))
  }
  out <- d |>
    count(.data[[v]]) |>
    rename(sub_var = 1) |>
    mutate(
      raw_mean = n / sum(n),
      var      = v,
      rowname  = paste0(v, sub_var)
    ) |>
    relocate(rowname)
  if (!is.null(dsn)) {
    wgtd <- svymean(
      as.formula(paste0("~", v)),
      design = dsn
    ) |>
      as_tibble(rownames = "rowname") |>
      rename(
        wmean = mean,
        wse   = SE
      )
    out <- left_join(
      out,
      wgtd,
      by = "rowname"
    )
  }
  return(out)
}

# summarize continuous variables
bin_sum <- function(d, v, dsn) {
  out <- d |>
    summarize(
      raw_mean = mean(.data[[v]]),
      sd = sd(.data[[v]])
    ) |>
    mutate(
      var       = v,
      sub_var   = v,
      n         = length(d[[v]]),
      "rowname" = v
    ) |>
    relocate(rowname)
  if (!is.null(dsn)) {
    wgtd <- svymean(
      as.formula(paste0("~", v)),
      design = dsn
    ) |>
      as_tibble(rownames = "rowname") |>
      rename(
        wmean = mean,
        wse   = all_of(v)
      )
    out <- left_join(
      out,
      wgtd,
      by = "rowname"
    )
  }
  return(out)
}

# single function that summarizes variable based on its class
qsum <- function(d, v, dsn) {
  if (d |> pull({{ v }}) |> (\(x) is.character(x) | is.factor(x))()) {
    return(ch_sum(d, v, dsn))
  }
  if (d |> pull({{ v }}) |> is.numeric()) {
    return(bin_sum(d, v, dsn))
  }
}

ch_sum(merged, "age_verbose", dsn)
d <- merged
v <- "age_verbose"
qsum(d = merged, v = "age", dsn = dsn)
qsum(d = merged, v = "ethn", dsn = dsn)

# function that loops over a vector of variables and summarizes based on its class
lqsum <- function(d, vars, dsn) {
  out <- list()
  pb <- txtProgressBar(max = length(vars), width = 50, style = 3)
  for (i in seq_along(vars)) {
    out[[i]] <- qsum(merged, vars[i], dsn)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  bind_rows(out) |>
    relocate(rowname, var, sub_var)
}

# test
vars <- c(
  "age", "age_verbose", "sex", "race_eth",
  "cancer", "diabetes", "cad", "anxiety", "depression", "smoker"
)

ukb_sum <- lqsum(
  merged,
  vars,
  dsn
)

ukb_sum |>
  filter(rowname == "age") |>
  pull(raw_mean)
