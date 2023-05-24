# libraries --------------------------------------------------------------------
suppressPackageStartupMessages({
  library(qs)
  library(data.table)
  library(tidyverse)
  library(glue)
  library(cowplot)
  library(htmlwidgets)
  library(plotly)
  library(patchwork)
})

# load data --------------------------------------------------------------------
mgi <- qread("~/Downloads/mgi_prevs.qs")
setnames(mgi, c("n", "N", "prev_unweighted", "prev_weighted", "se"), paste0("mgi_", c("n", "N", "prev_unweighted", "prev_weighted", "se")))

ukb <- qread("~/Downloads/ukb_prevs.qs")
setnames(ukb, c("n", "N", "prev_unweighted", "prev_weighted", "se"), paste0("ukb_", c("n", "N", "prev_unweighted", "prev_weighted", "se")))

aou <- qread("~/Downloads/aou_prevs.qs")
setnames(aou, c("n", "N", "prev_unweighted", "prev_weighted", "se"), paste0("aou_", c("n", "N", "prev_unweighted", "prev_weighted", "se")))


# process data -----------------------------------------------------------------
merged <- Reduce(
  \(x, y) merge.data.table(x, y, by = "phecode", all = TRUE),
  list(mgi, ukb, aou)
)

merged[, `:=`(
  mgi_aou = mgi_prev_unweighted / aou_prev_unweighted,
  mgi_ukb = mgi_prev_unweighted / ukb_prev_unweighted,
  ukb_aou = ukb_prev_unweighted / aou_prev_unweighted,
  aou_ukb = aou_prev_unweighted / ukb_prev_unweighted,
  mgi_pr  = mgi_prev_weighted / mgi_prev_unweighted,
  ukb_pr  = ukb_prev_weighted / ukb_prev_unweighted,
  aou_pr  = aou_prev_weighted / aou_prev_unweighted
)]

plot_dat <- merged[
  aou_n >= 20 & ukb_n >= 20 & mgi_n >= 20
][!is.na(mgi_aou) & !is.na(aou_ukb) & !is.na(mgi_ukb), ]

pheinfo <- fread("https://raw.githubusercontent.com/maxsal/public_data/main/phewas/Phecode_Definitions_FullTable_Modified.txt",
  colClasses = "character", showProgress = FALSE
)
phe_groups <- fread("https://raw.githubusercontent.com/maxsal/public_data/main/phewas/phecat_alt_colors.txt")
phe_group_cols <- phe_groups[, color]
names(phe_group_cols) <- phe_groups[, group]

plot_dat <- merge.data.table(
  plot_dat,
  pheinfo[, .(phecode = paste0("X", phecode), description, group, color)]
)

# plots ------------------------------------------------------------------------
(full_plot <- plot_dat |>
  ggplot(aes(
    x = mgi_aou, y = ukb_aou, color = group, label = phecode,
    label1 = mgi_prev_unweighted, label2 = ukb_prev_unweighted, label3 = aou_prev_unweighted,
    label4 = description
  )) +
  annotate(
    geom = "rect",
    xmin = 0.75, xmax = 1.25, ymin = 0, ymax = Inf,
    fill = "#0072B2", alpha = 0.2
  ) +
  annotate(
    geom = "rect",
    xmin = 0.9, xmax = 1.1, ymin = 0, ymax = Inf,
    fill = "#0072B2", alpha = 0.5
  ) +
  annotate(
    geom = "rect",
    ymin = 0.75, ymax = 1.25, xmin = 0, xmax = Inf,
    fill = "#009E73", alpha = 0.2
  ) +
  annotate(
    geom = "rect",
    ymin = 0.9, ymax = 1.1, xmin = 0, xmax = Inf,
    fill = "#009E73", alpha = 0.5
  ) +
  geom_abline() +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 1) +
  scale_color_manual(values = phe_group_cols) +
  geom_point() +
  scale_x_continuous(
    trans = "log10",
    breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
    labels = as.character(c(0.001, 0.01, 0.1, 1, 10, 100))
  ) +
  scale_y_continuous(
    trans = "log10",
    breaks = c(0.001, 0.01, 0.1, 1, 10, 100),
    labels = as.character(c(0.001, 0.01, 0.1, 1, 10, 100))
  ) +
  labs(
    color = "Phecode category",
    x = "MGI / AOU prevelance ratio",
    y = "UKB / AOU prevelance ratio",
    title = "Comparison of phecode prevalence ratios in AOU, MGI, and UKB",
    caption = stringr::str_wrap(glue(
      "Shaded areas represent ranges where MGI (blue) and UKB (green) ",
      "are within 10% (darker) and 25% (lighter) than AOU. There are ",
      "{format(nrow(plot_dat), big.mark = ',')} phecodes with prevalences in ",
      "all three cohorts. Phecodes in AOU required a minimum count of 20."
    ),
    width = 200
    )
  ) +
  cowplot::theme_half_open() +
  theme(
    plot.caption = element_text(hjust = 0)
  )
)

(zoom_plot <- plot_dat |>
  ggplot(aes(x = mgi_aou, y = ukb_aou, color = group)) +
  annotate(
    geom = "rect",
    xmin = 0.75, xmax = 1.25, ymin = 0, ymax = Inf,
    fill = "#0072B2", alpha = 0.2
  ) +
  annotate(
    geom = "rect",
    xmin = 0.9, xmax = 1.1, ymin = 0, ymax = Inf,
    fill = "#0072B2", alpha = 0.5
  ) +
  annotate(
    geom = "rect",
    ymin = 0.75, ymax = 1.25, xmin = 0, xmax = Inf,
    fill = "#009E73", alpha = 0.2
  ) +
  annotate(
    geom = "rect",
    ymin = 0.9, ymax = 1.1, xmin = 0, xmax = Inf,
    fill = "#009E73", alpha = 0.5
  ) +
  geom_abline() +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 1) +
  scale_color_manual(values = phe_group_cols) +
  coord_cartesian(xlim = c(0.1, 2), ylim = c(0.1, 2)) +
  geom_point() +
  scale_x_continuous(
    trans = "log10",
    breaks = c(0.1, .5, 1, 1.5, 2),
    labels = as.character(c(0.1, .5, 1, 1.5, 2))
  ) +
  scale_y_continuous(
    trans = "log10",
    breaks = c(0.1, .5, 1, 1.5, 2),
    labels = as.character(c(0.1, .5, 1, 1.5, 2))
  ) +
  labs(
    color = "Phecode category",
    x = "MGI / AOU prevelance ratio",
    y = "UKB / AOU prevelance ratio",
    title = "Comparison of phecode prevalence ratios in AOU, MGI, and UKB",
    caption = str_wrap(glue(
      "Shaded areas represent ranges where MGI (blue) and UKB (green) ",
      "are within 10% (darker) and 25% (lighter) than AOU. There are ",
      "{format(nrow(plot_dat), big.mark = ',')} phecodes with prevalences in ",
      "all three cohorts. Phecodes in AOU required a minimum count of 20."
    ),
    width = 200
    )
  ) +
  cowplot::theme_half_open() +
  theme(
    plot.caption = element_text(hjust = 0)
  )
)

# save output ------------------------------------------------------------------
ggsave(
  plot = full_plot,
  filename = "~/Downloads/full_prev_plot.pdf",
  width = 10 * 1.68, height = 10,
  device = cairo_pdf
)

ggsave(
  plot = zoom_plot,
  filename = "~/Downloads/zoom_prev_plot.pdf",
  width = 10 * 1.68, height = 10,
  device = cairo_pdf
)

# full_plotly <- full_plot |>
#   ggplotly(tooltip = c(
#     "label", "x", "y",
#     "label1", "label2", "label3", "label4"
#   ))

# htmlwidgets::saveWidget(
#   widget = full_plotly, # the plotly object
#   file = "~/Downloads/full_prevalence_plot.html", # the path & file name
#   selfcontained = TRUE # creates a single html file
# )

# prevalence ratio plots -------------------------------------------------------
prevalence_ratio_plot <- function(prevalence_data,
                                  ratio_var,
                                  breaks = 1 * 10^c(-2, -1, 0, 1, 3),
                                  savefile = NULL,
                                  title_text,
                                  y_axis_label = "Prevalence ratio",
                                  out_height = 7) {

  out <- as_tibble(plot_dat) |>
    select(phecode, mgi_aou, ukb_aou, mgi_ukb, aou_ukb, mgi_pr, ukb_pr, aou_pr) |>
    pivot_longer(cols = -1) |>
    filter(name == ratio_var) |>
    left_join(
      pheinfo[, .(phecode = paste0("X", phecode), group, color)],
      by = "phecode"
    ) |>
    ggplot(aes(x = group, y = value, fill = group)) +
    geom_hline(yintercept = 1) +
    geom_boxplot() +
    scale_fill_manual(values = phe_group_cols) +
    scale_y_continuous(
      trans = "log10",
      breaks = breaks,
      labels = ifelse(breaks >= 1000,
        format(round(breaks),
          nsmall = 0, scientific = FALSE,
          big.mark = ","
        ),
        as.character(breaks)
      )
    ) +
    coord_cartesian(ylim = c(min(breaks), max(breaks))) +
    labs(
      x = "",
      y = y_axis_label,
      title = str_wrap(title_text, width = 75)
    ) +
    theme_half_open() +
    theme(
      axis.text.x = element_text(angle = -45, hjust = 0, vjust = 0.5),
      legend.position = "none",
      plot.margin = margin(t = 3, r = 10, b = 3, l = 3, unit = "mm")
    )

  if (!is.null(savefile)) {
    ggsave(
      filename = savefile,
      plot     = out,
      height   = out_height,
      width    = out_height * 1.618,
      device   = cairo_pdf
    )
  }
  return(out)
}

mgi_aou_plot <- prevalence_ratio_plot(
  prevalence_data = plot_dat,
  title_text = "Boxplots of phecode unweighted prevalence ratios in Michigan Genomics Initiative and All of Us",
  y_axis_label = "Prevelance ratio (MGI / AOU)",
  ratio_var = "mgi_aou",
  savefile = "~/Downloads/mgi_aou_prevalence_plot.pdf"
)

aou_ukb_plot <- prevalence_ratio_plot(
  prevalence_data = plot_dat,
  title_text = "Boxplots of phecode unweighted prevalence ratios in All of Us and UK Biobank",
  y_axis_label = "Prevelance ratio (AOU / UKB)",
  ratio_var = "aou_ukb",
  breaks = 1 * 10^c(-3, -2, -1, 0, 1, 3),
  savefile = "~/Downloads/aou_ukb_prevalence_plot.pdf"
)

mgi_ukb_plot <- prevalence_ratio_plot(
  prevalence_data = plot_dat,
  title_text = "Boxplots of phecode unweighted prevalence ratios in Michigan Genomics Initiative and UK Biobank",
  y_axis_label = "Prevelance ratio (MGI / UKB)",
  ratio_var = "mgi_ukb",
  breaks = 1 * 10^c(-1, 0, 1, 3),
  savefile = "~/Downloads/mgi_ukb_prevalence_plot.pdf"
)

patched <- (mgi_aou_plot + labs(title = "MGI / AOU") + theme(axis.text.x = element_blank())) /
  (aou_ukb_plot + labs(title = "AOU / UKB") + theme(axis.text.x = element_blank())) /
  (mgi_ukb_plot + labs(title = "MGI / UKB"))

patched2 <- patched +
  plot_annotation(
    # title = "Boxplots of unweighted phecode prevalence ratios in All of Us, the Michigan Genomics Initiative, and UK Biobank",
    tag_levels = "A"
  )

ggsave(
  plot = patched2,
  filename = "~/Downloads/patched_unweighted_pr_plot.pdf",
  width = 7, height = 7 * 1.68,
  device = cairo_pdf
)


mgi_pr_plot <- prevalence_ratio_plot(
  prevalence_data = plot_dat,
  title_text = "Boxplots of phecode weighted / unweighted prevalence ratios in Michigan Genomics Initiative",
  y_axis_label = "Prevelance ratio\n(Weighted / Unweighted)",
  ratio_var = "mgi_pr",
  breaks = 1 * 10^c(-1, 0, 1),
  savefile = "~/Downloads/mgi_w_prevalence_plot.pdf"
)

ukb_pr_plot <- prevalence_ratio_plot(
  prevalence_data = plot_dat,
  title_text = "Boxplots of phecode weighted / unweighted prevalence ratios in UK Biobank",
  y_axis_label = "Prevelance ratio\n(Weighted / Unweighted)",
  ratio_var = "ukb_pr",
  breaks = 1 * 10^c(-1, 0, 1),
  savefile = "~/Downloads/ukb_w_prevalence_plot.pdf"
)

aou_pr_plot <- prevalence_ratio_plot(
  prevalence_data = plot_dat,
  title_text = "Boxplots of phecode weighted / unweighted prevalence ratios in All of Us",
  y_axis_label = "Prevelance ratio\n(Weighted / Unweighted)",
  ratio_var = "aou_pr",
  breaks = 1 * 10^c(-1, 0, 1),
  savefile = "~/Downloads/aou_w_prevalence_plot.pdf"
)

patched3 <- (aou_pr_plot + labs(title = "All of Us") + theme(axis.text.x = element_blank())) /
  (mgi_pr_plot + labs(title = "Michigan Genomics Initiative") + theme(axis.text.x = element_blank())) /
  (ukb_pr_plot + labs(title = "UK Biobank"))

patched4 <- patched3 +
  plot_annotation(
    # title = "Boxplots of weighted vs unweighted phecode prevalence ratios in All of Us, the Michigan Genomics Initiative, and UK Biobank",
    tag_levels = "A"
  )

ggsave(
  plot = patched4,
  filename = "~/Downloads/patched_weighted_pr_plot.pdf",
  width = 7, height = 7 * 1.68,
  device = cairo_pdf
)
