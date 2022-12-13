library(data.table)
library(ppcor)

quick_mgi_partial_correlation <- function(x) {
  cca <- complete.cases(sub_pim[, .SD, .SDcols = c(x[1], x[2])])

  if (sum(cca) == dim(sub_pim)[1]) {
    a <- pcor.test(sub_pim[[x[1]]],
                   sub_pim[[x[2]]],
                   x1_mgi,
                   method = "pearson")
  } else if (sum(cca) < dim(sub_pim)[1] && sum(cca != 0)) {
    a <- pcor.test(sub_pim[cca, .SD, .SDcols = x[1]],
                   sub_pim[cca, .SD, .SDcols = x[2]],
                   x2_mgi[cca, ],
                   method = "pearson")
  }
  
  cbind(
    data.table(to = x[1], from = x[2]),
    as.data.table(a)
  )
  
}
