##############################
### Calculate odds ratios  ###
##############################
## based on code by Lauren  ##
## Beesley ###################
##############################

###### !!! NOT READY FOR USE !!! ########

# libraries and options --------------------------------------------------------
options(stringsAsFactors=F)

library("data.table",quietly=T)
library(tableone)
library(knitr)
library(xtable)
library(grDevices)
library(ggplot2)
library(gridExtra)
library(GGally)
library(corrplot)
library(network)
library(sna)
library(dplyr)
library(ppcor)
library(progress)
mytheme <- gridExtra::ttheme_default(
  core = list(fg_params=list(cex = 0.5)),
  colhead = list(fg_params=list(cex = 0.5)),
  rowhead = list(fg_params=list(cex = 0.5)))

# load functions ---------------------------------------------------------------
invisible(lapply(list.files("fn/", full.names = TRUE), source))
use_geno <- TRUE

# specifications ---------------------------------------------------------------
mgi_version <- "20210318" # mgi data version

file_paths <- get_files(mgi_version = mgi_version)

# read in data -----------------------------------------------------------------

# phecode indicator matrix (PEDMASTER_0)
pim0 <- data.table::fread(file_paths[["mgi"]]$pim0_file)
data.table::setnames(pim0, old = "IID", new = "id")

# demographics data
demo <- data.table::fread(file_paths[["mgi"]]$demo_file)[, .(id = Deid_ID, age = Age, alive = as.numeric(AliveYN == "Y"), dead_dsb = Deceased_DaysSinceBirth, ethn = EthnicityName, marital = MaritalStatusCode, sex = Sex, race = RaceName)]

# pc data 
if (use_geno == TRUE) {
  pcs <- data.table::fread("/net/junglebook/magic_data/MGI_GenotypeData_Freeze4/MGI_Freeze4_Sample_Info_60215samples.txt")
  data.table::setnames(pcs, old = "Deid_ID", new = "id")
  
  # merge pcs in demo data
  merged <- merge(demo, pcs)
} else {
  merged <- demo
}

sub_pim <- pim0[id %in% merged[, id]]    # subset pim
merged  <- merged[id %in% sub_pim[, id]] # subset merged

# covariates -------------------------------------------------------------------
if (use_geno == TRUE) {
  X1_MGI = merged[, .(
    Age = age,
    FEMALE = as.numeric(sex == "F"),
    PC1, PC2, PC3, PC4)]
  X2_MGI = merged[, .(Age = age, PC1, PC2, PC3, PC4)]
} else {
  X1_MGI = merged[, .(Age = age, FEMALE = as.numeric(sex == "F"))]
  X2_MGI = merged[, .(Age = age)]
}

# MGI Partial Correlations -----------------------------------------------------
phecodes <- names(sub_pim[, !c("id")])
sub_pim <- sub_pim[, 2:ncol(sub_pim)]

SAVED = c()
pb <- progress_bar$new(total = length(phecodes) - 1)
for (i in 1:(length(phecodes) - 1)) {
  for (j in (i+1):length(phecodes)) {
    
    CCA = complete.cases(sub_pim[, .SD, .SDcols = c(i, j)])
    
    if (sum(CCA) == dim(sub_pim)[1]) {
      a <- pcor.test(sub_pim[[i]], sub_pim[[j]], X1_MGI, method = "pearson")
      SAVED = rbind(SAVED, c(i, j, a[1, 1]))
    } else if (sum(CCA) < dim(sub_pim)[1] & sum(CCA != 0)) {
      a <- pcor.test(sub_pim[CCA, .SD, .SDcols = i], sub_pim[CCA, .SD, .SDcols = j], X2_MGI[CCA, ], method = "pearson")
      SAVED = rbind(SAVED, c(i, j, a[1, 1]))
    }
  }
  pb$tick()
}

SAVED2 = data.frame(to = phecodes[SAVED[,1]], from = phecodes[SAVED[,2]], pcor = SAVED[,3])

output_file_name <- paste0("mgi_phenome_partial_correlations_",
                           ifelse(use_geno == TRUE, "w_geno_pcs_", ""),
                           version, ".txt")

data.table::fwrite(x = SAVED2,
                   file = paste0("data/", version, "/", output_file_name), sep = "\t")
