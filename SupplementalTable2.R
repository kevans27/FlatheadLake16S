countTableSmall <- 
  read.csv("~/FlatheadMicrobes/SubsampledData/small_updatedNames.count_table", 
                            stringsAsFactors = FALSE, row.names = 1)

smallSeqs <- rownames(countTableSmall)
nameKey <- data.frame("LongNames" = smallSeqs, 
                      "ShortNames" = countTableSmall$Representative_Sequence)

countTableSmall$blastname <- NULL
countTableSmall$Representative_Sequence <- NULL

for (i in 1){
  fullset <- read.csv('~/FlatheadPublic/hydrolab.csv', na.strings="")
  fullset <- fullset[fullset$Site == 'Flathead Lake, Midlake Deep',]
  units <- fullset[1,]
  fullset <- fullset[-1,]
  tempStart <- as.Date("2016-01-01")
  tempEnd  <- as.Date("2019-01-31")
  
  main <- "Temp"
  val <- main
  
  df <- fullset[complete.cases(fullset[,main]),]
  df <- df[seq(13969, nrow(df)),]
  
  df$Date <- as.Date(df$Date, format = "%m/%d/%y")
  suppressWarnings(df[, main] <- as.numeric(as.character(df[, main])))
  df <- df[(!is.na(df[, main])),]
  df <- df[(!is.na(df$Depth)),]
  df$Depth <- round(as.numeric(as.character(df$Depth)))
  df$Month <- as.numeric(format(df$Date, "%m"))
  
  df <- df[df$Depth %in% c(1, 5, 10, 15, 20, 50, 90),]
  
  tempDF <- df[df$Date >= tempStart & df$Date <= tempEnd,]
  tempDF <- tempDF[!is.na(tempDF$Time),]
  surfTemp <- tempDF[tempDF$Depth == 1,]
  
  DCMDat <- read.csv("~/FlatheadPublic/DCMData.csv", header = TRUE, stringsAsFactors = FALSE)
  DCMDat$Date <- as.Date(DCMDat$Date)
  DCMDat$Date.1 <- as.Date(DCMDat$Date.1)
  DCMDat$Date.1.1 <- as.Date(DCMDat$Date.1.1)
  DCMDat <- DCMDat[order(DCMDat$Date),]
  DCMDat <- aggregate(DCM~Date, DCMDat, FUN = mean)
  
  
  mldData <- 
    read.csv("~/FlatheadPublic/MLDData.csv", stringsAsFactors = FALSE, row.names = 1)
  mldData$Date <- as.Date(mldData$Date)
  mldData <- mldData[order(mldData$Date),]
}

taxTable <- read.csv("RawData/prok.98.taxonomy")

saveSeqNamesSmall <- smallSeqs[order(rowSums(countTableSmall), decreasing = TRUE)]
countTableSmall <- countTableSmall[order(rowSums(countTableSmall), 
                                         decreasing = TRUE),]
trimmedSNS <- saveSeqNamesSmall[1:500]
trimmedCTS <- countTableSmall[1:500,]
trimmedSS <- rowSums(countTableSmall)[1:500]
trimmedSTax <- taxTable[match(saveSeqNamesSmall, taxTable$seqID)[1:500],]

relAbund <- data.frame("SeqID" = trimmedSNS, "RelAbund" = trimmedSS/157/100)

colDataSmall <- data.frame("Depth" = rep(NA, ncol(trimmedCTS)), "SurfTemp" = NA, 
                           "Date" = as.Date(NA), "DOY" = NA)

for (i in 1){
  colnames(trimmedCTS) <- gsub("Sept", "Sep", colnames(trimmedCTS))
  colnames(trimmedCTS) <- gsub("June", "Jun", colnames(trimmedCTS))
  colnames(trimmedCTS) <- gsub("July", "Jul", colnames(trimmedCTS))
  
  twoPartsS <- strsplit(colnames(trimmedCTS), "_")
  subDatesS <- sapply(twoPartsS, `[[`, 2)
  subDatesS <- as.Date(subDatesS, format = "%d%b%y")
  depthsS <- sapply(twoPartsS, `[[`, 3)
  
  colDataSmall$Depth <- depthsS
  colDataSmall$SurfTemp <- tempDF[match(subDatesS, tempDF$Date), "Temp"]
  colDataSmall$Date <- subDatesS
  colDataSmall$DOY <- as.integer(format(subDatesS, "%j"))
  
  #colDataSmall <- merge(colDataSmall, mldData, by = "Date", all.x = TRUE)
  
  colDataSmall[colDataSmall$Depth == "DCM", "Depth"] <- 
    DCMDat[match(colDataSmall[colDataSmall$Depth == "DCM", "Date"], DCMDat$Date), "DCM"]
  
  colDataSmall$Depth <- gsub("m", "", colDataSmall$Depth)
  colDataSmall$Depth <- as.numeric(as.character(colDataSmall$Depth))
  
  removeSmall <- which(rowSums(is.na(colDataSmall)) > 0)
  colDataSmall <- colDataSmall[-removeSmall,]
  trimmedCTS <- trimmedCTS[, -removeSmall]
}
colDataSmall$Layer <- "Epi"
colDataSmall[colDataSmall$Depth > 30, "Layer"] <- "Hypo"

colDataSmall$MixingStat <- "Strat"
colDataSmall[(colDataSmall$Date > as.Date("2016-12-01") & colDataSmall$Date < as.Date("2017-05-01")) | (colDataSmall$Date > as.Date("2017-12-01") & colDataSmall$Date < as.Date("2018-05-15")), "MixingStat"] <- "Mixed"

noFullyMixed <- colDataSmall[colDataSmall$MixingStat == "Strat",]
noMixedAbund <- trimmedCTS[, colDataSmall$MixingStat == "Strat"]

onlyMixed <- colDataSmall[colDataSmall$MixingStat == "Mixed",]
mixedAbund <- trimmedCTS[, colDataSmall$MixingStat == "Mixed"]

library(DESeq2)


###Strat

deseqStrat <- DESeqDataSetFromMatrix(noMixedAbund, noFullyMixed, 
                                             ~ Layer)
deStrat <- DESeq(deseqStrat)
otuStratResults <- results(deStrat)

deseqTableStrat <- data.frame("baseMean" = otuStratResults$baseMean, 
                              "log2FoldChange" = otuStratResults$log2FoldChange, 
                              "lfcSE" = otuStratResults$lfcSE, 
                              "stat" = otuStratResults$stat, 
                              "pvalue" = otuStratResults$pvalue, 
                              "padj" = otuStratResults$padj)
rownames(deseqTableStrat) <- rownames(otuStratResults[1])

stratDE <- merge(deseqTableStrat, nameKey, by.x = "row.names", by.y = "LongNames")
stratDE$Rank <- as.numeric(as.character(gsub("Otu", "", stratDE$ShortNames)))

plotMA(otuStratResults)


###Mixed

deseqStrat <- DESeqDataSetFromMatrix(mixedAbund, onlyMixed, 
                                     ~ Layer)
deStrat <- DESeq(deseqStrat)
otuStratResults <- results(deStrat)

deseqTableStrat <- data.frame("baseMean" = otuStratResults$baseMean, 
                              "log2FoldChange" = otuStratResults$log2FoldChange, 
                              "lfcSE" = otuStratResults$lfcSE, 
                              "stat" = otuStratResults$stat, 
                              "pvalue" = otuStratResults$pvalue, 
                              "padj" = otuStratResults$padj)
rownames(deseqTableStrat) <- rownames(otuStratResults[1])

mixedDE <- merge(deseqTableStrat, nameKey, by.x = "row.names", by.y = "LongNames")
mixedDE$Rank <- as.numeric(as.character(gsub("Otu", "", mixedDE$ShortNames)))

plotMA(otuStratResults)


###Mixed vs strat

deseqStrat <- DESeqDataSetFromMatrix(trimmedCTS, colDataSmall, 
                                     ~ MixingStat)
deStrat <- DESeq(deseqStrat)
otuStratResults <- results(deStrat)

deseqTableStrat <- data.frame("baseMean" = otuStratResults$baseMean, 
                              "log2FoldChange" = otuStratResults$log2FoldChange, 
                              "lfcSE" = otuStratResults$lfcSE, 
                              "stat" = otuStratResults$stat, 
                              "pvalue" = otuStratResults$pvalue, 
                              "padj" = otuStratResults$padj)
rownames(deseqTableStrat) <- rownames(otuStratResults[1])

compDE <- merge(deseqTableStrat, nameKey, by.x = "row.names", by.y = "LongNames")
compDE$Rank <- as.numeric(as.character(gsub("Otu", "", compDE$ShortNames)))

plotMA(otuStratResults)




##In stacked bar
targOTUs <- c(1, 5, 2, 7, 9, 8, 17, 14, 46, 15, 23, 114, 3, 36, 28, 84, 6, 22, 10, 12, 
              25, 4, 105, 11, 20, 24, 26, 42, 18, 29, 51, 118, 120)

diffOTUsStrat <- stratDE[stratDE$Rank %in% targOTUs, c("ShortNames", "log2FoldChange",
                                                       "padj", "Rank")]
colnames(diffOTUsStrat) <- c("Names", "log2_Strat", "padj_Strat", "Rank")
diffOTUsMix <- mixedDE[mixedDE$Rank %in% targOTUs, c("ShortNames", "log2FoldChange",
                                                     "padj")]
colnames(diffOTUsMix) <- c("Names", "log2_Mixed", "padj_Mixed")
diffOTUsSeas <- compDE[compDE$Rank %in% targOTUs, c("ShortNames", "log2FoldChange",
                                                    "padj")]
colnames(diffOTUsSeas) <- c("Names", "log2_Seas", "padj_Seas")

merge0 <- merge(diffOTUsStrat, diffOTUsMix, by = "Names")
merge1 <- merge(merge0, diffOTUsSeas, by = "Names")

sortedDESeqTable <- merge1[match(targOTUs, merge1$Rank),]
sortedDESeqTable$Rank <- NULL

write.csv(sortedDESeqTable, "deseqTriTable.csv", quote = FALSE)
