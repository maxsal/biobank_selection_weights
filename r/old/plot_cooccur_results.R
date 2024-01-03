### libraries ----------
library(data.table)
library(ggplot2)
library(plotly)
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
  ~source(paste0("~/Dropbox (University of Michigan)/projects/dissertation/aim_one/fn/", .x))
)

### load data ----------
d <- fread(paste0("~/Downloads/mgi_X", gsub("X", "", out_phe),"_t", t, "_20210318_results.txt"))[, `:=` (
  phecode = gsub("X", "", phecode),
  shape = ifelse(beta > 0, "Up", "Down"),
  log10p = -log10(p_value),
  phe_num = as.numeric(gsub("X", "", phecode))
  )][]
p <- fread("~/Downloads/Phecode_Definitions_FullTable_Modified.txt", colClasses = "character")[, `:=` (
  groupnum = as.numeric(groupnum)
)]

### prep data ----------
merged <- merge.data.table(
  d,
  p[, .(phecode, description, group, groupnum, color)],
  all.x = TRUE,
  by = "phecode"
)[order(groupnum, phe_num), num := 1:.N]
merged[, slide_num := num + ((groupnum - 1)*30)]
thresh <- -log10(0.05/merged[, .N])
merged[, sig := ifelse(p_value < 0.05/merged[, .N], "Significant", "Not significant")]
merged[, group := stringr::str_to_sentence(group)]

### ggplot-based plots ----------
ggmanplot_log10p(data = merged, phe_data = p, outcome = out_phe)
ggmanplot_beta(data = merged, phe_data = p, outcome = out_phe)


### plotly-based plots ----------
plotlyman_log10p(data = merged, phe_data = p, outcome = out_phe)
plotlyman_beta(data = merged, phe_data = p, outcome = out_phe)

### table
DT::datatable(merged[order(p_value), .(Phecode = phecode, Description = description, Beta = round(beta, 3), SEbeta = round(se_beta, 3), `P-value` = formatC(p_value, format = "e", digits = 2), Significant = sig, `Phecode category` = group)])
