library(vegan)
`%notin%` <- Negate(`%in%`)

#Ready tables
for (i in 1){
  eukCount <- read.table("~/FlatheadMicrobes/RawData/euks.count_table", 
                         sep = "\t", header = TRUE, 
                         stringsAsFactors = FALSE)
  rownames(eukCount) <- eukCount$Representative_Sequence
  eukCount$Representative_Sequence <- NULL
  relEukCount <- eukCount[, grep("05m|10m|DCM|50m|90m", colnames(eukCount))]
  relEukCount$Representative_Sequence <- paste0("Euk", sprintf("%07d", 1:nrow(relEukCount)))
  
  smallEukAbs <- relEukCount[, grep("3u", names(relEukCount), invert = TRUE)]
  largeEukAbs <- relEukCount[, grep("02u", names(relEukCount), invert = TRUE)]
  
  prokCount <- read.table("~/FlatheadMicrobes/RawData/proks.count_table", sep = "\t",
                          header = TRUE, stringsAsFactors = FALSE)
  prokTax <- read.csv("~/FlatheadMicrobes/RawData/prok.98.taxonomy", 
                      stringsAsFactors = FALSE)
  prokKey <- read.csv("~/FlatheadMicrobes/RawData/prok.name.key", stringsAsFactors = FALSE,
                      header = TRUE)

  relProkCount <- prokCount[, grep("05m|10m|DCM|50m|90m|Representative_Sequence", 
                                 colnames(prokCount))]
  
  relProkCount$blastname <- 
    prokKey[relProkCount$Representative_Sequence == prokKey$qseqids, "blastname"]
  
  rownames(relProkCount) <- relProkCount$blastname
  relProkCount$blastname <- NULL

  absCounts <- colSums(relProkCount[, grep("05m|DCM", colnames(relProkCount))])
  smallProkAbs <- relProkCount[, grep("3u", colnames(relProkCount), invert = TRUE)]
  largeProkAbs <- relProkCount[, grep("02u", colnames(relProkCount), invert = TRUE)]
  
  smallAbs <- rbind(smallEukAbs, smallProkAbs)
  largeAbs <- rbind(largeEukAbs, largeProkAbs)
}

#Subsample
for (i in 1){
  largeAbsCounts <- colSums(largeAbs[,grep("Repres", colnames(largeAbs), invert = TRUE)])
  #Min 12894, FMP188_5Dec17_10m_3u; Max 64813, FMP112_7June17_90m_3u
  #Subsample to 12000
  smallAbsCounts <- colSums(smallAbs[,grep("Repres", colnames(smallAbs), invert = TRUE)])
  #Min 2169, FMP137_7Aug17_DCM_02u; Max 96569, FMP10_9Sept16_90m_02u 
  #Remove FMP137_7Aug17_DCM_02u, subsample to 12000
  
  trimLAC <- largeAbs[, grep("3u", colnames(largeAbs))]
  trimSAC <- smallAbs[, colnames(smallAbs) %notin% c("FMP137_7Aug17_DCM_02u")]
  trimSAC <- trimSAC[, grep("02u", colnames(trimSAC))]
  
  
  tSmallAbs <- t(as.matrix(trimSAC))
  tLargeAbs <- t(as.matrix(trimLAC))
  
  rSmallAbs <- rrarefy(tSmallAbs, 12000)
  colnames(rSmallAbs) <- rownames(smallAbs)
  rrSmall <- as.data.frame(t(rSmallAbs))
  rrSmall$Representative_Sequence <- smallAbs$Representative_Sequence
  rrSmall <- rrSmall[rowSums(rrSmall[, grep("Represent", colnames(rrSmall), invert = TRUE)]) > 0,]
  
  
  rLargeAbs <- rrarefy(tLargeAbs, 12000)
  colnames(rLargeAbs) <- rownames(largeAbs)
  rrLarge <- as.data.frame(t(rLargeAbs))
  rrLarge$Representative_Sequence <- largeAbs$Representative_Sequence
  rrLarge <- rrLarge[rowSums(rrLarge[, grep("Represent", colnames(rrLarge), invert = TRUE)]) > 0,]
}


write.csv(rrLarge, file = "~/FlatheadMicrobes/SubsampledData/rrBothLarge.count_table", 
          quote = FALSE)
write.csv(rrSmall, file = "~/FlatheadMicrobes/SubsampledData/rrBothSmall.count_table", 
          quote = FALSE)

