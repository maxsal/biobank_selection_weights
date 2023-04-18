# weighting variable subfields from https://github.com/sjoerdvanalten/UKBWeights
# getRanges, reformatUKB and other code borrowed from https://github.com/umich-cphds/createUKBphenome/blob/master/scripts/function.reformatUKB.r

options(stringsAsFactors=F)
library("data.table")
library("optparse")
library("intervals")

# data.table with the fields and their descriptions
all.tables <- fread("./data/Fields_in_Available_Data.txt",header=T)

# data.table with the merged baskets
data.file <- "./data/merged_baskets.txt"

ukb_data_columns <- scan(data.file,what=character(0),sep="\t",nlines=1,quiet=T)

# merge consecutive numbers to intervals and collapse
getRanges <- function(colnumbers,collapse=T){
	cx <- clusters(colnumbers,1)
	ux <- colnumbers[which(!colnumbers %in% unlist(cx))]
	ix <- sapply(cx,function(x) paste(range(x),collapse="-"))
	out <- c(ux,ix)
	out <- out[order(as.numeric(gsub("\\-.+","",out)))]
	paste(out,collapse=",")
}

weight_var_subfields <- c(
    "f.eid", "f.31.0.0", "f.34.0.0", "f.54.0.0", "f.20074.0.0", 
    "f.20075.0.0", "f.709.0.0", "f.6138.0.0", "f.6138.0.1", "f.6138.0.2", 
    "f.6138.0.3", "f.6138.0.4", "f.6138.0.5", "f.2178.0.0", "f.6142.0.0", 
    "f.680.0.0", "f.728.0.0", "f.21000.0.0", "f.40000.0.0", "f.845.0.0"
)

weight_var_subfields2 <- c(
    "id", "f.31.0.0", "f.34.0.0", "f.54.0.0", "20074-0.0", 
    "20075-0.0", "709-0.0", "6138-0.0", "6138-0.1", "6138-0.2", 
    "6138-0.3", "6138-0.4", "6138-0.5", "2178-0.0", "6142-0.0", 
    "680-0.0", "728-0.0", "f.21000.0.0", "f.40000.0.0", "845-0.0"
)

weight_var_fields <- sapply(strsplit(weight_var_subfields, "\\."),
        "[[", 2)

#reformatUKB <- function(fields,dataCoding=F,onlyInfo=F){
#	# print info about fields
#	fieldinfo <- all.tables[which(all.tables$'Field ID' %in% fields & !is.na(all.tables$basket)),]
#
#	if(dataCoding){
#		require("RCurl")
#		require("htmltidy")
#		require("XML")
#		fieldinfo$URL <- URL <- paste0("https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=",fieldinfo$'Field ID')
#		dataCodes <- NULL
#		for(i in 1:nrow(fieldinfo)){
#			u <- URL[i]
#			doc.raw <- getURL(u)
#			doc <- tidy_html(doc.raw)
#			html <- htmlTreeParse(doc, useInternal = TRUE)
#			txt <- xpathApply(html, "//body//text()[not(ancestor::script)][not(ancestor::style)][not(ancestor::noscript)]", xmlValue)
#			txt <- unlist(unlist(txt))
#			codedData <- which(grepl("Data-Coding",txt))
#			if(length(codedData)==1){
#				dataCodes <- c(dataCodes,txt[which(grepl("Data-Coding",txt))+1])
#			} else {
#				dataCodes <- c(dataCodes,NA)
#			}
#		}
#		fieldinfo$'Data Coding' <- paste0("https://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=",dataCodes)
#		print(data.table(fieldinfo))
#		if(onlyInfo) return(return(list('info'=data.table(fieldinfo,dataCodes))))
#	} else {
#		print(fieldinfo[,1:4])
#	}
#	
#	
#	# Extract all fields 
#	allfields <- fieldinfo$subfield_id[which(fieldinfo$subfield_id != "n/a")]
#	
#	selected_columns <- c("id",unique(unlist(strsplit(allfields,","))))	
# 	print(paste("Reading all entries across",length(selected_columns),"columns"))
#
#	colnumbers <- which(ukb_data_columns %in% selected_columns)
#	
#	reformatted <- fread(cmd=paste("cut -f",getRanges(colnumbers),data.file),header=T,
#		colClasses="character")
#	
#	fieldpattern <- "f\\.(.+)\\..+\\..+"
#	entrypattern <- "f\\..+\\.(.+\\..+)"
#
#	keep <- which(rowSums(is.na(reformatted[,-1]) | reformatted[,-1] == "",na.rm = T) != ncol(reformatted[,-1]))
#	
#	reformatted <- reformatted[keep,]
#
#	entries <- names(reformatted)[-1]
#	entryNumbers <-	gsub(entrypattern,"\\1",entries)
#	fieldCols <- unique(gsub(fieldpattern,"\\1",entries))
#	uentries <- unique(entryNumbers)
#
#	# reformat cancer registry (allow multiple entries by person)
#	reformatted2 <- list()
#	for(e in uentries){
#		newEntries <- entries[which(entryNumbers == e)]
#		newRows <- reformatted[,c("id",newEntries),with=F]
#		setnames(newRows,newEntries,gsub(fieldpattern,"\\1",newEntries))
#		missingFields <- fields[which(!fields %in% names(newRows))]
#		for(missingField in missingFields){
#			newRows[[as.character(missingField)]] <- NA
#		}
#		reformatted2[[e]] <- newRows[,c("id",fields),with=F]
#	}
#
#	reformatted2 <- rbindlist(reformatted2)
#
#	# remove empty entries	
#	keep <- which(rowSums(is.na(reformatted2[,-1]) | reformatted2[,-1] == "",na.rm = T) != ncol(reformatted2[,-1]))
#	reformatted2	<- reformatted2[keep,]
#	print(paste(nrow(reformatted2),"lines after filtering empty lines"))
#	if(!dataCoding){
#		return(reformatted2)
#	} else {
#		return(list('data'=reformatted2,'info'=data.table(fieldinfo,dataCodes)))
#	}
#}

fieldinfo <- all.tables[which(all.tables$'Field ID' %in% fields & !is.na(all.tables$basket)),]

	# Extract all fields 
allfields <- fieldinfo$subfield_id[which(fieldinfo$subfield_id != "n/a")]

selected_columns <- c("id",unique(unlist(strsplit(allfields,","))))	
print(paste("Reading all entries across",length(selected_columns),"columns"))
colnumbers <- which(ukb_data_columns %in% selected_columns)

reformatted <- fread(cmd=paste("cut -f",getRanges(colnumbers),data.file),header=T,
	colClasses="character")

fieldpattern <- "f\\.(.+)\\..+\\..+"
entrypattern <- "f\\..+\\.(.+\\..+)"
keep <- which(rowSums(is.na(reformatted[,-1]) | reformatted[,-1] == "",na.rm = T) != ncol(reformatted[,-1]))
	
reformatted <- reformatted[keep,]

reformatted2 <- reformatted[, ..weight_var_subfields2]

keep <- which(rowSums(is.na(reformatted2[,-1]) | reformatted2[,-1] == "",na.rm = T) != ncol(reformatted2[,-1]))
reformatted2
reformatted2 <- reformatted2[keep,]

reformatted3 <- copy(reformatted2)
setnames(reformatted3, names(reformatted3), weight_var_subfields)

fwrite(x = reformatted3, "data/UKB-weight-variables.tab", sep = "\t")
