ggmanplot_beta <- function(data, phe_data, thresh = NULL, outcome) {
  
  # check that necessary columns exist in dataset
  if (!all(c("slide_num", "sig", "group", "shape", "beta") %in% names(data))) {
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
    ggplot(aes(x = slide_num, y = beta, fill = group, color = group)) +
    geom_hline(yintercept = 0, color = "red") +
    geom_point(aes(shape = sig)) +
    scale_color_manual(values = col_vals) +
    scale_fill_manual(values = col_vals) +
    scale_shape_manual(values = c("Significant" = 21, "Not significant" = 1)) +
    scale_x_continuous(breaks = x_axis_dat[, mean],
                       labels = paste0("<span style=\"color: ", x_axis_dat[, color], "\">",
                                       x_axis_dat[, group],
                                       "</span>")) +
    labs(
      title = paste0("X", gsub("X", "", outcome), ": ", phe_data[phecode == gsub("X", "", outcome), description]),
      subtitle = paste0("Threshold = ", t, " years"),
      x = "",
      y = "Beta coefficient",
      caption = paste0("N = ", format(data[, .N], big.mark = ","), ". Filled shapes indicate statistical significance at multiple testing threshold (p < ", sci_p_thresh, ").")
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1),
      plot.caption = element_text(hjust = 0)
    )
}


