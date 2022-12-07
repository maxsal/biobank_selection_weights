plotlyman_beta <- function(data, phe_data, outcome) {
  
  # check that necessary columns exist in dataset
  if (!all(c("slide_num", "log10p", "group", "shape", "beta") %in% names(data))) {
    not_in <- c("slide_num", "log10p", "color", "group", "shape")[!(c("slide_num", "log10p", "group", "shape") %in% names(date))]
    stop(paste0("Check that the following variables are in the dataset: ", paste0(not_in, collapse = ", ")))
  }
  
  symbols = c("Significant" = "circle", "Not significant" = "circle-open")
  
  # gather data for x-axis labels
  x_axis_dat <- unique(data[, .(mean = mean(slide_num), color), group])
  
  col_vals <- x_axis_dat[, color]
  names(col_vals) <- x_axis_dat[, group]
  
  plot_ly(data = data[!is.na(beta)]) |>
    layout(shapes = list(
      type = "line", 
      x0 = 0, 
      x1 = 1, 
      xref = "paper",
      y0 = 0, 
      y1 = 0, 
      line = list(color = "red")
    )) |>
    add_trace(
      x = ~slide_num, y = ~beta,
      type = "scatter", mode = "markers",
      color = ~group, colors = col_vals,
      symbol = ~sig, symbols = symbols,
      marker = list(size = 10),
      text = ~paste0("Phecode: ", phecode, "<br>",
                     "Description: ", description, "<br>",
                     "Group: ", snakecase::to_title_case(group), "<br>",
                     "Beta: ", round(beta, 2), "<br>",
                     "P-value: ", signif(p_value, 3)),
      hovertemplate = "%{text}",
      name = ""
    ) |>
    layout(
      title = list(
        text = paste0("X", gsub("X", "", outcome), ": ", phe_data[phecode == gsub("X", "", outcome), description]),
        x = 0.1
      ),
      yaxis = list(
        title = "Beta coefficient"
      ),
      xaxis = list(
        title = "",
        tickangle = -45,
        ticktext = as.list(paste0("<span style='color:", x_axis_dat[, color], "'>", x_axis_dat[, group], "</span>")),
        tickvals = as.list(x_axis_dat[, mean])
      ), 
      showlegend = FALSE,
      margin = list(l = 50, r = 50, b = 75, t = 50),
      annotations = list(x = 1, y = -0.1, text = paste0("N = ", data[, .N]),
                         xref='paper', yref='paper', showarrow = F, 
                         xanchor='right', yanchor='auto', xshift=0, yshift=0,
                         font = list(size = 10))
    )
}


