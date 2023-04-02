# functions for obtaining quick summaries from datasets
library(data.table)

summarizer_func <- function(x, var_name = NULL) {
  if (is.numeric(x) & all(x %in% c(0:1, NA_real_))) {
    # dichotomous
    out <- data.table(
      val = as.character(c(0, 1)),
      N   = c(sum(x == 0, na.rm = TRUE), sum(x == 1, na.rm = TRUE)),
      nas = length(is.na(x))
    )[,
      prop := round(N * 100 / length(na.omit(x)), 1)
      ][, `:=` (
            print      = paste0(format(round(prop, 1), nsmall = 1), " (", trimws(format(N, big.mark = ",")), ")"),
            print_unit = "% (n)"
        )
      ][]
  } else if (is.numeric(x)) {
    # numeric
    out <- data.table(
      mean   = mean(x, na.rm = TRUE),
      sd     = sd(x, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      N      = length(na.omit(x)),
      nas    = length(is.na(x))
    )[, `:=` (
      print = paste0(
        format(round(mean, 2), big.mark = ",", nsmall = 2),
        " (",
        format(round(sd, 2), big.mark = ",", nsmall = 2),
        ")"),
        print_unit = "mean (sd)"
        )
      ][]
  } else if (is.factor(x)) {
    out <- data.table(
      val = levels(x),
      N   = map_dbl(levels(x), \(level) sum(x == level, na.rm = TRUE)),
      nas = length(is.na(x))
    )[,
      prop := round(N * 100 / length(na.omit(x)), 1)
      ][, `:=` (
        print      = paste0(format(round(prop, 1), nsmall = 1), " (", trimws(format(N, big.mark = ",")), ")"),
        print_unit = "% (n)"
        )
        ][]
  } else if (is.character(x)) {
    out <- data.table(
      val = unique(x),
      N   = map_dbl(unique(x), \(level) sum(x == level, na.rm = TRUE)),
      nas = length(is.na(x))
    )[,
      prop := round(N * 100 / length(na.omit(x)), 1)
    ][, `:=` (
      print      = paste0(format(round(prop, 1), nsmall = 1), " (", trimws(format(N, big.mark = ",")), ")"),
      print_unit = "% (n)"
      )
    ][]
  } else {
    out <- data.table()
  }
  if (!is.null(var_name)) {
    out[, variable := var_name]
  }
  return(out)
}

summarizer <- function(x, col_names = NULL, as_tibble = FALSE) {
  if (is.null(col_names)) {
    out <- rbindlist(mapply(summarizer_func, x, names(x)), use.names = TRUE, fill = TRUE)
  } else {
    out <- rbindlist(mapply(summarizer_func, x[, ..col_names], col_names), use.names = TRUE, fill = TRUE)
  }
  if (as_tibble) {
    tibble::as_tibble(out)
  } else {
    out
  }
}
