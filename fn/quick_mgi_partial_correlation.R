quick_mgi_partial_correlation <- function(x) {
  CCA = complete.cases(sub_pim[, .SD, .SDcols = c(x[1], x[2])])
  
  if (sum(CCA) == dim(sub_pim)[1]) {
    a <- pcor.test(sub_pim[[x[1]]], sub_pim[[x[2]]], X1_MGI, method = "pearson")
  } else if (sum(CCA) < dim(sub_pim)[1] & sum(CCA != 0)) {
    a <- pcor.test(sub_pim[CCA, .SD, .SDcols = x[1]], sub_pim[CCA, .SD, .SDcols = x[2]], X2_MGI[CCA, ], method = "pearson")
  }
  cbind(
    data.table(to = x[1], from = x[2]),
    as.data.table(a)
  )
}