quick_mgi_partial_correlation <- function(x) {
  cca <- complete.cases(sub_pim[, .SD, .SDcols = c(x[1], x[2])])

  if (sum(cca) == dim(sub_pim)[1]) {
    a <- ppcor::pcor.test(sub_pim[[x[1]]],
                          sub_pim[[x[2]]],
                          x1_mgi,
                          method = "pearson")
  } else if (sum(cca) < dim(sub_pim)[1] && sum(cca != 0)) {
    a <- ppcor::pcor.test(sub_pim[cca, .SD, .SDcols = x[1]],
                          sub_pim[cca, .SD, .SDcols = x[2]],
                          x2_mgi[cca, ],
                          method = "pearson")
  }
  cbind(
    as.data.table::data.table(to = x[1], from = x[2]),
    data.table::as.data.table(a)
  )
}
