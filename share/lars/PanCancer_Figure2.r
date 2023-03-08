options(stringsAsFactors=F)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)

version <- "20200211"
source("/net/junglebook/home/larsf/Rfunctions/function.num2phecode.r")

load(file=paste0("/net/junglebook/home/larsf/Projects/panCancer/results/manuscript/All_PRS_Evaluations_",version,".Rsav"))
# There is one lassosum PRS with less than 5 SNPs
prseval <- prseval[SNPs >= 5,]


# Figure 2 / four PRS distributions
source("/net/junglebook/home/larsf/Rfunctions/function.fig_label.r")

# load phenome data
PEDMASTER_UKB <- fread("/net/junglebook/home/larsf/Projects/PRS_PIPELINE/data/Phenomes/UKB_20181102/UKB_20181102_PEDMASTER_MATCHED.txt")
setnames(PEDMASTER_UKB,c("FID","IID"),c("#FAM_ID","IND_ID"))
PEDMASTER_MGI <- fread("/net/junglebook/home/larsf/Projects/PheWAS_ICD/MGI_20190429/MGI_20190429_PEDMASTER_MATCHED.txt")
phenomes <- list(
	'UKB'=list(
		'covars'="Sex,birthYear,genotyping.array,PC1,PC2,PC3,PC4",
		'keepSamples'=readLines("/net/junglebook/home/larsf/Projects/PRS_PIPELINE/data/Phenomes/UKB_20181102/UKB_20181102_PEDMASTER_MATCHED_KeepSamples.txt"),
		'PEDMASTER'=PEDMASTER_UKB
	),
	'MGI'=list(
		'covars'="SEX,AGE,BATCH,PC1,PC2,PC3,PC4",
		'keepSamples'=PEDMASTER_MGI[[1]],
		'PEDMASTER'=PEDMASTER_MGI
	))

library(RColorBrewer)
library(tools) #toTitleCase
pcols <- brewer.pal(6, "Set1")
tcols <- grDevices::rgb(t(grDevices::col2rgb(pcols)), alpha = 50, maxColorValue = 255)
prsdistribution <- function(i){
	if(prseval$Top_Underpowered[i]) next
	phecode <- prseval$phecode[i]
	covars <- strsplit(phenomes[[prseval$genomePRS[i]]]$covars,",")[[1]]	
	if(prseval$sex[i] != "both") covars <- covars[-1]
	# matched results
	cohort <- phenomes[[prseval$genomePRS[i]]]$PEDMASTER[,c("IND_ID",phecode,gsub("X","S",phecode),covars),with=F]
	cohort <- cohort[IND_ID %in% phenomes[[prseval$genomePRS[i]]]$keepSamples,]
	cohort <- cohort[which(!is.na(cohort[,gsub("X","S",phecode),with=F])),]

	# load PRS
	file.prs <- paste0("/net/junglebook/home/larsf/Projects/panCancer/results/prsweb_",version,"/",prseval$prswebprefix[i],"_PRS.txt")
	prs <- fread(file.prs)
	cohort <- merge(cohort,prs,by.x="IND_ID",by.y="IID")
	
	#Normalized across the full analytical dataset
	cohort[,PRS:=scale(PRS)]

	samplesize <- formatC(table(cohort[[phecode]]),big.mark=",")
	plottitle <- paste0(toTitleCase(prseval$phecode_description[i]),"\n",
		prseval$genomePRS[i],": ",samplesize[2]," Cases & ",samplesize[1]," Matched Controls")

	# TOP 1, 2, 5, 10 and 25%
	probs <- c(0.25,0.1,0.05,0.02,0.01)

	topbreaks <- c(as.numeric(quantile(cohort$PRS, probs = 1-probs)),max(cohort$PRS))
	print("Mean difference")
	print(mean(cohort$PRS[which(cohort[[2]] == 1)]) - mean(cohort$PRS[which(cohort[[2]] == 0)]))
	
	d_controls <- density(cohort$PRS[which(cohort[[2]] == 0)],n=2^15)
	d_cases <- density(cohort$PRS[which(cohort[[2]] == 1)],n=2^15)
	xlim <- range(d_cases$x,d_controls$x)
	ylim <- range(0,0.55)#d_cases$y,d_controls$y)

	par(las=1)
	plot(d_controls,lwd=2,col="black",xlab=paste0("PRS"),main=plottitle,xlim=xlim,ylim=ylim*1.2)
	lines(d_cases,lwd=2,col="red")
	liney <- par("usr")[3] / length(probs)
	linex <- xlim[2]
	for(j in 2:length(topbreaks)){
		x1 <- topbreaks[j-1]
		x2 <- topbreaks[j]
		polygon(c(x1,d_controls$x[d_controls$x>=x1 & d_controls$x<=x2],x2),  c(0,d_controls$y[d_controls$x>=x1 & d_controls$x<=x2],0), col=tcols[j-1], border = NA)
		polygon(c(x1,d_cases$x[d_cases$x>=x1 & d_cases$x<=x2],x2),  c(0,d_cases$y[d_cases$x>=x1 & d_cases$x<=x2],0), col=tcols[j-1], border = NA)
		abline(v=x1,col=pcols[j-1])		
		lines(c(x1,linex),rep(liney*(j-2),2),col=pcols[j-1],lwd=2)	
	}
	tops <- probs*100
	or <- signif(prseval[i,paste0("Top_",probs,"_OR"),with=F],2)
	ci1 <- signif(prseval[i,paste0("Top_",probs,"_CI1"),with=F],2)
	ci2 <- signif(prseval[i,paste0("Top_",probs,"_CI2"),with=F],2)
	legend("topleft",c("Cases","Controls"),col=c("red","black"),lwd=2,bg="white")
	legend("topright",paste0("Top ",tops,"%: OR ",or," (",ci1,", ",ci2,")"),
		col=pcols,bg="white",pch=15,title="Top versus Rest")
}




pdf(paste0("/net/junglebook/home/larsf/Projects/panCancer/results/manuscript/Figures/Figure2_PRSdistributions.pdf"),
	width = 17.4 / cm(1), height = 17.4 / cm(1), pointsize = 8)
	layout(matrix(1:4,nrow=2,byrow=T))
	selected_prefix <- "Onco_iCOGS_Overall_BRCA"
	selected_method <- "Lassosum"
	index <- prseval[,.I[prefix == selected_prefix & method == selected_method],]

	prsdistribution(index[1])
	fig_label("a",cex=1.5)
	
	prsdistribution(index[2])
	fig_label("b",cex=1.5)

	selected_prefix <- "GWAS_Catalog_r2019-05-03_X204.12"
	selected_method <- "P&T"
	index <- prseval[,.I[prefix == selected_prefix & method == selected_method],]

	prsdistribution(index[2])
	fig_label("c",cex=1.5)

	prsdistribution(index[1])
	fig_label("d",cex=1.5)

dev.off()


