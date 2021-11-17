library(coop)

countTableLarge <- read.csv("~/FlatheadMicrobes/SubsampledData/rrBothLarge.count_table", 
                            stringsAsFactors = FALSE, row.names = 1)
countTableSmall <- read.csv("~/FlatheadMicrobes/SubsampledData/rrBothSmall.count_table", 
                            stringsAsFactors = FALSE, row.names = 1)

taxTable <- read.csv("RawData/prok.98.taxonomy")

for (i in 1){
  countTableLarge$Representative_Sequence <- NULL
  countTableSmall$Representative_Sequence <- NULL
  
  colnames(countTableLarge) <- sub("[^_]*_(.*)", "\\1", colnames(countTableLarge))
  colnames(countTableLarge) <- gsub("_[^_]+$", "\\1", colnames(countTableLarge))
  colnames(countTableLarge) <- gsub("Sept", "Sep", colnames(countTableLarge))
  colnames(countTableLarge) <- gsub("June", "Jun", colnames(countTableLarge))
  
  colnames(countTableSmall) <- sub("[^_]*_(.*)", "\\1", colnames(countTableSmall))
  colnames(countTableSmall) <- gsub("_[^_]+$", "\\1", colnames(countTableSmall))
  colnames(countTableSmall) <- gsub("Sept", "Sep", colnames(countTableSmall))
  colnames(countTableSmall) <- gsub("June", "Jun", colnames(countTableSmall))
}
rm(i)

trimmedLarge <- countTableLarge[rowSums(countTableLarge)>1,]
trimmedSmall <- countTableSmall[rowSums(countTableSmall)>1,]

taxTable[,"tclass"] <- gsub("\\s*\\([^\\)]+\\)","",as.character(taxTable$class))

largeClass <- list()
smallClass <- list()
for (i in unique(taxTable$tclass)){
  targs <- taxTable[taxTable$tclass == i, "seqID"]
  largeDF <- trimmedLarge[rownames(trimmedLarge) %in% targs,]
  smallDF <- trimmedSmall[rownames(trimmedSmall) %in% targs,]
  
  largeClass <- append(largeClass, list(largeDF))
  smallClass <- append(smallClass, list(smallDF))
}

names(largeClass) <- unique(taxTable$tclass)
names(smallClass) <- unique(taxTable$tclass)

largeClass <- largeClass[sapply(largeClass, nrow)>1]
smallClass <- smallClass[sapply(smallClass, nrow)>1]








###Adding euks
###
eukTax <- read.table("~/FlatheadMicrobes/RawData/euks.taxonomy", sep = "\t", 
                     header = FALSE, stringsAsFactors = FALSE)
processTax <- function(taxTable){
  taxDat <- data.frame(do.call("rbind", strsplit(as.character(taxTable$V2), ";", 
                                                 fixed = TRUE)))
  newTax <- cbind(taxTable$V1, taxDat)
  colnames(newTax) <- c("seqID", "domain", "supergroup", "phylum", "class", "subclass", 
                        "order", "suborder", "family", "genus", "species")
  return(newTax)
}
eukTax <- processTax(eukTax)

eukTax[,"tclass"] <- gsub("\\s*\\([^\\)]+\\)","",as.character(eukTax$class))

largeClassEuks <- list()
smallClassEuks <- list()
for (i in unique(eukTax$tclass)){
  targs <- eukTax[eukTax$tclass == i, "seqID"]
  largeDF <- trimmedLarge[rownames(trimmedLarge) %in% targs,]
  smallDF <- trimmedSmall[rownames(trimmedSmall) %in% targs,]
  
  largeClassEuks <- append(largeClassEuks, list(largeDF))
  smallClassEuks <- append(smallClassEuks, list(smallDF))
}

names(largeClassEuks) <- unique(eukTax$tclass)
names(smallClassEuks) <- unique(eukTax$tclass)

numEuksLarge <- length(largeClassEuks)
numEuksSmall <- length(smallClassEuks)

largeClassEuks <- largeClassEuks[sapply(largeClassEuks, nrow)>1]
smallClassEuks <- smallClassEuks[sapply(smallClassEuks, nrow)>1]

largeClass <- append(largeClass, largeClassEuks)
smallClass <- append(smallClass, smallClassEuks)


#Based largely on equation A2
lehrmenStab2000 <- function(df){
  bioVal <- mean(colSums(df))
  miniMatt <- t(df)
  
  covarMatt <- sum(covar(miniMatt))
  sumCV <- sum(covarMatt)
  
  tempStab <- bioVal/(sumCV**(1/2))
  
  return(tempStab)
}

largeStabs <- c()
largeAbund <- c()
smallStabs <- c()
smallAbund <- c()
for (i in 1:length(largeClass)){
  df <- largeClass[[i]]
  stabilityVal <- lehrmenStab2000(df)
  
  largeStabs <- c(largeStabs, stabilityVal)
  largeAbund <- c(largeAbund, mean(colSums(df)))
}
for (i in 1:length(smallClass)){
  df <- smallClass[[i]]
  stabilityVal <- lehrmenStab2000(df)
  
  smallStabs <- c(smallStabs, stabilityVal)
  smallAbund <- c(smallAbund, mean(colSums(df)))
}



names(largeStabs) <- names(largeClass)
names(largeAbund) <- names(largeClass)
names(smallStabs) <- names(smallClass)
names(smallAbund) <- names(smallClass)