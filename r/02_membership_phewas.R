# perform a membership phewas (m = 1 when in MGI)
# author:   max salvatore
# date:     20230418

# libraries, functions, and options --------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(MatchIt)
  library(logistf)
  library(glue)
  library(qs)
  library(optparse)
  library(purrr)
})

set.seed(61787)

lapply(list.files("fn/", full.names = TRUE), source) |> # load functions
  invisible()

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

time_thresholds <- as.numeric(strsplit(opt$time_thresholds, ",")[[1]])

## extract file paths
file_paths <- get_files(
  mgi_version = opt$mgi_version,
  ukb_version = opt$ukb_version
)

# read data --------------------------------------------------------------------
## mgi
### pim0
mgi_pim0 <- fread(file_paths[["mgi"]][["pim0_file"]])
setnames(mgi_pim0, "IID", "id")
mgi_pim0[, id := as.character(id)]
mgi_pim0_more_than_10 <- names(mgi_pim0[, !"id"])[colSums(mgi_pim0[, !"id"], na.rm = TRUE) > 10]

### cov
mgi_cov <- read_qs(file = glue("data/private/mgi/{opt$mgi_version}/data_{opt$mgi_version}_{opt$mgi_cohort}.qs"))[
  ,
  ehr_years := round((LastDaySinceBirth - FirstDaySinceBirth) / 365.25, 1)
]
setnames(mgi_cov, "DeID_PatientID", "id")
mgi_cov[, id := as.character(id)]

## ukb
### pim0
ukb_pim0 <- fread(file_paths[["ukb"]][["pim0_file"]])
setnames(ukb_pim0, "IID", "id")
replace_missing(ukb_pim0)
ukb_pim0[, id := as.character(id)]
ukb_pim0_more_than_10 <- names(ukb_pim0[, !"id"])[colSums(ukb_pim0[, !"id"], na.rm = TRUE) > 10]

### cov
ukb_cov <- fread(file_paths[["ukb"]][["demo_file"]])[, female := as.numeric(sex == "Female")][!is.na(dob)]
ukb_cov[, id := as.character(id)]

# select codes in both cohorts -------------------------------------------------
removr <- function(x, y) {
  x[!(x %in% y)]
}
both_phecodes <- removr(intersect(mgi_pim0_more_than_10, ukb_pim0_more_than_10), c("IID", "id"))

# membership model -------------------------------------------------------------
stacked_cov <- rbindlist(list(
  mgi_cov[, data := "MGI"],
  ukb_cov[, data := "UKB"]
), use.names = TRUE, fill = TRUE)[, case := as.numeric(data == "MGI")]
stacked_cov[case == 1, last_dsb := LastDaySinceBirth]

mem_mod <- glm(case ~ age_at_last_diagnosis + female + ehr_years, data = stacked_cov, family = "binomial")
summary(mem_mod)

# membership phewas ------------------------------------------------------------
keep_cols <- c("id", both_phecodes)
stacked_pim <- rbindlist(list(
  mgi_pim0[, ..keep_cols][, data := "MGI"],
  ukb_pim0[, ..keep_cols][, data := "UKB"]
), use.names = TRUE)

both_ids <- intersect(stacked_pim[, id], stacked_cov[, id])

pim_data <- stacked_pim[id %in% both_ids]
cov_data <- stacked_cov[id %in% both_ids]

# 1. identify analytic phecodes
possible_phecodes <- names(pim_data)[names(pim_data) %in% both_phecodes]
phecodes_to_consider <- melt(
  pim_data[
    ,
    ..possible_phecodes
  ][
    ,
    lapply(.SD, \(x) sum(x, na.rm = TRUE))
  ],
  variable.name = "phecode", value.name = "n", id.vars = character()
)[
  n >= 10,
  as.character(phecode)
]

# 2. merge covariates
merged <- merge.data.table(
  pim_data[, !c("case", "data")],
  cov_data,
  by = "id",
  all.x = TRUE
)

mem_phewas <- list()
pb <- txtProgressBar(max = length(phecodes_to_consider), style = 3, width = 50)
for (i in seq_along(phecodes_to_consider)) {
  mem_phewas[[i]] <- quick_cooccur_mod(
    dat        = merged,
    covs       = c("last_dsb", "female", "ehr_years"),
    ex_code    = phecodes_to_consider[i],
    mod_type   = "glm"
  )
  setTxtProgressBar(pb, i)
}
close(pb)
mem_phewas <- rbindlist(mem_phewas)

save_qs(
  x = mem_phewas,
  file = glue(
    "results/mgi/{opt$mgi_version}/",
    "mgi_ukb_membership_phewas_results.qs"
  )
)

# plot phewas results ----------------------------------------------------------
### libraries ----------
library(data.table)
library(ggplot2)
library(ggtext)

#### !!! ########
### need to fix the group names in plotting functions
#################

### specs ----------
out_phe <- "155"
t <- 5

### source plotting functions ----------
purrr::walk(
  list.files("~/Dropbox (University of Michigan)/projects/dissertation/aim_one/fn/"),
  ~ source(paste0("~/Dropbox (University of Michigan)/projects/dissertation/aim_one/fn/", .x))
)

### load data ----------
mem_phewas[, `:=`(
  phecode = gsub("X", "", phecode),
  shape = ifelse(beta > 0, "Up", "Down"),
  log10p = -log10(p_value),
  phe_num = as.numeric(gsub("X", "", phecode))
)]
p <- fread("data/public/Phecode_Definitions_FullTable_Modified.txt", colClasses = "character")[, `:=`(
  groupnum = as.numeric(groupnum)
)]

### prep data ----------
merged <- merge.data.table(
  mem_phewas,
  p[, .(phecode, description, group, groupnum, color)],
  all.x = TRUE,
  by = "phecode"
)[order(groupnum, phe_num), num := 1:.N]
merged[, slide_num := num + ((groupnum - 1) * 30)]
thresh <- -log10(0.05 / merged[, .N])
merged[, sig := ifelse(p_value < 0.05 / merged[, .N], "Significant", "Not significant")]
merged[, group := stringr::str_to_sentence(group)]

# gather data for x-axis labels
x_axis_dat <- unique(merged[, .(mean = mean(slide_num), color), group])

col_vals <- x_axis_dat[, color]
names(col_vals) <- x_axis_dat[, group]

data <- merged
phe_data <- p
sci_p_thresh <- formatC(1 / 10^(thresh), format = "e", digits = 2)

mem_phewas_plot <- data[!is.na(beta)] |>
  ggplot(aes(x = slide_num, y = log10p, fill = group, color = group)) +
  geom_hline(yintercept = thresh, color = "red") +
  annotate(geom = "text", x = x_axis_dat[, max(mean)], y = thresh - 0.5, label = sci_p_thresh, vjust = 1, color = "red", fontface = "italic", size = 3) +
  geom_point(aes(shape = shape)) +
  scale_color_manual(values = col_vals) +
  scale_fill_manual(values = col_vals) +
  scale_shape_manual(values = c("Up" = 24, "Down" = 25)) +
  scale_x_continuous(
    breaks = x_axis_dat[, mean],
    labels = paste0(
      "<span style=\"color: ", x_axis_dat[, color], "\">",
      x_axis_dat[, group],
      "</span>"
    )
  ) +
  labs(
    title = "MGI membership PheWAS plot",
    x = "",
    y = "-log10(p-value)",
    caption = paste0("N = ", format(data[, .N], big.mark = ","))
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1),
    plot.caption = element_text(hjust = 0)
  )

ggsave(
  plot = mem_phewas_plot,
  filename = "bin/membership_phewas_plot.pdf",
  width = 6, height = 6, device = cairo_pdf
)
