if (!"ms" %in% as.data.frame(installed.packages())[["Package"]]) remotes::install_github("maxsal/ms")
suppressPackageStartupMessages({
  library(ms)
  library(data.table)
  library(ggplot2)
  library(scales)
})

stacked_pca_plot <- function(x, cohort = "mgi") {
  
  if (!is.data.table(x)) {
    x <- as.data.table(x)
  }
  
  x[, stat := c("sd", "prop", "cum_prop")]
  
  plot_data <- melt(
    x,
    id.vars = "stat"
  )[, pc := as.numeric(gsub("PC", "", variable))][]
  
  cumulative_plot <- plot_data[stat == "cum_prop"] |>
    ggplot(aes(x = pc, y = value)) +
    geom_vline(xintercept = plot_data[stat == "cum_prop"][order(pc)][ value > 0.95, ][1, pc], linetype = 2, linewidth = 1, color = "goldenrod") +
    geom_vline(xintercept = plot_data[stat == "cum_prop"][order(pc)][ value > 0.99, ][1, pc], linetype = 2, linewidth = 1, color = "darkred") +
    geom_label(data = plot_data[stat == "cum_prop"][order(pc)][value > 0.95][1, ],
               aes(x = pc, y = 0.3, label = paste0("95%: ", format(pc, big.mark = ","))), color = "goldenrod") +
    geom_label(data = plot_data[stat == "cum_prop"][order(pc)][value > 0.99][1, ],
               aes(x = pc, y = 0.3, label = paste0("99%: ", format(pc, big.mark = ","))), color = "darkred") +
    geom_line(linewidth = 1) +
    scale_x_continuous(labels = scales::comma) +
    labs(
      title = "Cumulative proportion of variance explained",
      x     = "Principal component",
      y     = "Variance explained",
      caption = paste0("N PCs: ", format(max(plot_data[, pc], na.rm = TRUE), big.mark = ","))
    ) +
    theme_ms(font_family = "sans")
    # theme_minimal()
  
  proportion_plot <- plot_data[stat == "prop"][value > 0.01] |>
    ggplot(aes(x = pc, y = value)) +
    geom_point(size = 2) +
    geom_line(linewidth = 1) +
    labs(
      title = "Proportion of variance explained by principal components",
      x = "Principal component",
      y = "Variance explained",
      caption = "Only showing PCs that explain at least 1% of variance"
    ) +
    theme_ms(font_family = "sans")
    # theme_minimal() +
    # theme(
    #   plot.caption = element_text(hjust = 0)
    # )
  
  patched <- proportion_plot / cumulative_plot
  
  patched + plot_annotation(
    title = glue::glue("PCA of {toupper(cohort)} phenome"),
    theme = theme(plot.title = element_text(size = 18, face = "bold")),
    tag_levels = 'A') &
    theme(
      plot.tag.position = c(0, 0.98),
      plot.tag = element_text(hjust = 0.5, vjust = 0.5, face = "bold")
    )
  
}
