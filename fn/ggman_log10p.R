ggmanplot_log10p <- function(data, phe_data, thresh = NULL, outcome) {
  
  # check that necessary columns exist in dataset
  if (!all(c("slide_num", "log10p", "group", "shape", "beta") %in% names(data))) {
    not_in <- c("slide_num", "log10p", "color", "group", "shape")[!(c("slide_num", "log10p", "group", "shape") %in% names(date))]
    stop(paste0("Check that the following variables are in the dataset: ", paste0(not_in, collapse = ", ")))
  }
  
  # check that threshold is present - if not, calculate
  if (is.null(thresh)) {
    thresh <- -log10(0.05/data[, .N])
  }
  
  sci_p_thresh <- formatC(1/10^(thresh), format = "e", digits = 2)
  
  # gather data for x-axis labels
  x_axis_dat <- unique(data[, .(mean = mean(slide_num), color), group])
  
  col_vals <- x_axis_dat[, color]
  names(col_vals) <- x_axis_dat[, group]
  
  data[!is.na(beta)] |>
    ggplot(aes(x = slide_num, y = log10p, fill = group, color = group)) +
    geom_hline(yintercept = thresh, color = "red") +
    annotate(geom = "text", x = x_axis_dat[, max(mean)], y = thresh - 0.5, label = sci_p_thresh, vjust = 1, color = "red", fontface = "italic", size = 3) +
    geom_point(aes(shape = shape)) +
    scale_color_manual(values = col_vals) +
    scale_fill_manual(values = col_vals) +
    scale_shape_manual(values = c("Up" = 24, "Down" = 25)) +
    scale_x_continuous(breaks = x_axis_dat[, mean],
                       labels = paste0("<span style=\"color: ", x_axis_dat[, color], "\">",
                                       x_axis_dat[, group],
                                       "</span>")) +
    labs(
      title = paste0("X", gsub("X", "", outcome), ": ", phe_data[phecode == gsub("X", "", outcome), description]),
      subtitle = paste0("Threshold = ", t, " years"),
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
}


