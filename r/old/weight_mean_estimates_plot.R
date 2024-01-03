library(data.table)
library(ggplot2)
library(colorblindr)
library(snakecase)
library(glue)
library(patchwork)

d <- fread("~/Downloads/cancer_weighted_summary_means.txt")

d

d[weight == "unweighted", weight := "Unweighted"]
d[weight == "mgi_can_dsn", weight := "Weight with cancer (direct)"]
d[weight == "mgi_ncan_cor_dsn", weight := "Weight with cancer (indirect)"]
d[weight == "mgi_dsn", weight := "Weight without cancer"]
d[, weight := factor(weight,
                     levels = c("Unweighted", "Weight without cancer", "Weight with cancer (indirect)", "Weight with cancer (direct)", "nhanes_design"))]

weighting_plot <- function(
    data,
    var,
    var_units = NULL,
    nhanes_est = "nhanes_design",
    round_n = 1,
    dodge_n = 0.2,
    displace_n = 0.5,
    rect_fill = "black",
    rect_alpha = 0.1,
    titl = NULL,
    y_axis_title = NULL
) {
  
  plot_data <- data[variable == var]
  
  plot_tab <- data.table(
    mean = plot_data[weight == "nhanes_design", mean],
    lower = plot_data[weight == "nhanes_design", mean] - (qnorm(0.975) * plot_data[weight == "nhanes_design", SE]),
    upper = plot_data[weight == "nhanes_design", mean] + (qnorm(0.975) * plot_data[weight == "nhanes_design", SE])
  )
  
  plot_data <- plot_data[weight != "nhanes_design"]
  
  plot_data[, `:=` (
    lower = mean - (qnorm(0.975) * SE),
    upper = mean + (qnorm(0.975) * SE)
  )]
  
  plot_data |>
    ggplot(aes(x = var, y = mean)) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = plot_tab[, lower], ymax = plot_tab[, upper]),
              fill = rect_fill, alpha = rect_alpha) +
    geom_pointrange(aes(color = weight, shape = weight,
                        ymin = lower, ymax = upper), position = position_dodge2(width = dodge_n)) +
    geom_text(aes(x = var, y = plot_data[, mean] - displace_n, label = format(round(plot_data[, mean], round_n), nsmall = round_n), color = weight),
              position = position_dodge2(width = dodge_n), show.legend = FALSE) +
    scale_color_OkabeIto(order = c(6, 2, 3, 7)) +
    scale_shape_manual(values = c(15:18)) +
    labs(
      title = ifelse(is.null(titl),
                     glue("{to_title_case(var)} in MGI cohort before and after NHANES weighting"),
                     titl),
      x = "",
      y = ifelse(
        is.null(y_axis_title),
        glue("{to_title_case(var)}{ifelse(!is.null(var_units), paste0(' (', var_units, ')'), '')}"),
        y_axis_title),
      caption = glue("NHANES est: {format(round(plot_tab[, mean], round_n), nsmall = round_n)} ({format(round(plot_tab[, lower], round_n), nsmall = round_n)}, {format(round(plot_tab[, upper], round_n), nsmall = round_n)})")
    ) +
    theme_minimal() +
    theme(
      legend.title = element_blank(),
      plot.caption = element_text(hjust = 0),
      legend.position = "top"
    )
}

age_plot <- weighting_plot(
  data = d,
  var = "age",
  var_units = "years",
  displace_n = 1,
  dodge_n = 0.9,
  titl = "Age"
) + theme(legend.position = "none")


cancer_plot <- weighting_plot(
  data = d,
  var = "cancer",
  var_units = "%",
  displace_n = 0.02,
  round_n = 2,
  dodge_n = 0.9,
  titl = "Cancer"
) + theme(legend.position = "none")


bmi_plot <- weighting_plot(
  data = d,
  var = "bmi",
  displace_n = 0.15,
  dodge_n = 0.9,
  y_axis_title = "BMI",
  titl = "BMI"
) + theme(legend.position = "none")


female_plot <- weighting_plot(
  data = d,
  var = "female",
  var_units = "%",
  round_n = 2,
  displace_n = 0.0075,
  dodge_n = 0.9,
  titl = "Female"
) + theme(legend.position = "none")


nhw_plot <- weighting_plot(
  data = d,
  var = "nhw",
  var_units = "%",
  round_n = 2,
  displace_n = 0.015,
  dodge_n = 0.9,
  titl = "Non-Hispanic White",
  y_axis_title = "Non-Hispanic White (%)"
) + theme(legend.position = "right")

patched <- age_plot + female_plot + nhw_plot + bmi_plot + cancer_plot + guide_area()+
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Unweighted and weighted mean values for selected variables",
    caption = glue(
      "Shaded area represents NHANES 95% confidence interval estimate.\n",
      "'Weights without cancer' weights were calculated using a Beta regression model adjusted for age, coronary heart diease, diabetes, BMI, smoking,\nand non-Hispanic White.\n",
      "'Weights with cancer (indirect)' are the 'weights without cancer' weights multiplied by a ratio of predicted cancer in MGI and NHANES.\n",
      "'Weights with cancer (direct)' uses the same model as 'weights without cancer' but additionally adjusts for cancer."
      ),
    theme =
      theme(
        plot.title = element_text(size = 18, face = "bold"),
        plot.caption = element_text(hjust = 0)))

ggsave(
  plot = patched,
  filename = "plots/weighted_variable_means_patched_plot.pdf",
  width = 9, height = 9, device = cairo_pdf
)
