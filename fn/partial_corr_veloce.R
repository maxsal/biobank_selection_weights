# super fast parallelize partial correlation
# ADAPTED FROM https://gist.github.com/melissamlwong/b489e7765a620c5fb1f3af0ece3024b1#file-bigcorpar-r-L12
# author: max salvatore
# date:   20230109

require(data.table)

partial_corr_veloce <- function(pim, ncore = detectCores()/2, covs1, covs2 = NULL) {
  require(data.table)
  require(doMC)
  registerDoMC(cores = ncore)
  column <- colnames(pim)
  cols   <- 1:ncol(pim)
  output <- foreach(i = cols) %dopar% {
    out <- list()
    for (j in cols) {
      if (j >= i) next
      if (j >= i) next
      cca <- complete.cases(a[, .SD, .SDcols = c(column[c(i, j)])])
      if (sum(cca) == nrow(a)) {
        cor_out <- pcor.test(a[, ..i],
                             a[, ..j],
                             covs1,
                             method = "pearson")
      } else if (sum(cca) < nrow(a) && sum(cca != 0)) {
        cor_out <- pcor.test(a[cca, ..i],
                             a[cca, ..j],
                             fifelse(!is.null(covs2),
                                     covs2[cca, ],
                                     covs1[cca, ]),
                             method = "pearson")
      }
      
      if (exists("cor_out")) {
        out[[j]] <- data.table(
          "phe1" = column[i],
          "phe2" = column[j],
          as.data.table(cor_out)
        )
      } else {
        out[[j]] <-  data.table(
          "phe1" = column[i],
          "phe2" = column[j]
        )
      }
    }
    rbindlist(out, use.names = TRUE, fill = TRUE)
  }
  gc(verbose = FALSE)
  return(output)
}