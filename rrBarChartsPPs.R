largeAbund <- read.csv("~/FlatheadMicrobes/SubsampledData/rrBothLarge.count_table", 
                       stringsAsFactors = FALSE, row.names = 1)
smallAbund <- read.csv("~/FlatheadMicrobes/SubsampledData/rrBothSmall.count_table", 
                       stringsAsFactors = FALSE, row.names = 1)
#Rename stuff
for (i in 1){
  colnames(largeAbund) <- gsub("Sept", "Sep", colnames(largeAbund))
  colnames(largeAbund) <- gsub("June", "Jun", colnames(largeAbund))
  colnames(largeAbund) <- gsub("July", "Jul", colnames(largeAbund))
  colnames(largeAbund) <- sub("[^_]*_(.*)", "\\1", colnames(largeAbund))
  colnames(smallAbund) <- gsub("Sept", "Sep", colnames(smallAbund))
  colnames(smallAbund) <- gsub("June", "Jun", colnames(smallAbund))
  colnames(smallAbund) <- gsub("July", "Jul", colnames(smallAbund))
  colnames(smallAbund) <- sub("[^_]*_(.*)", "\\1", colnames(smallAbund))
  
  largeAbund$Sequence <- NULL
  smallAbund$Sequence <- NULL
}

eukTax <- read.table("~/FlatheadMicrobes/RawData/euks.taxonomy", sep = "\t", header = FALSE, 
                     stringsAsFactors = FALSE)
processTax <- function(taxTable){
  taxDat <- data.frame(do.call("rbind", strsplit(as.character(taxTable$V2), ";", 
                                                 fixed = TRUE)))
  newTax <- cbind(taxTable$V1, taxDat)
  colnames(newTax) <- c("seqID", "domain", "supergroup", "phylum", "class", "subclass", 
                        "order", "suborder", "family", "genus", "species")
  return(newTax)
}
eukTax <- processTax(eukTax)
prokTax <- read.csv("~/FlatheadMicrobes/RawData/prok.98.taxonomy", 
                     stringsAsFactors = FALSE)
#Just save PPs
for(i in 1){
  targSeqs <- prokTax[grep("Oxyphoto", prokTax$class), "seqID"]
  cyanoTax <- prokTax[grep("Oxyphoto", prokTax$class), ]
  
  targSeqs <- c(targSeqs, as.character(eukTax$seqID))
  
  smallPPs <- smallAbund[rownames(smallAbund) %in% targSeqs, 
                         grep("02u", colnames(smallAbund))]
  largePPs <- largeAbund[rownames(largeAbund) %in% targSeqs, 
                         grep("3u", colnames(largeAbund))]
  
  #Make into relabund
  smallSums <- colSums(smallPPs)
  largeSums <- colSums(largePPs)
  
  largeRels <- sweep(largePPs, 2, largeSums, FUN = '/')
  smallRels <- sweep(smallPPs, 2, smallSums, FUN = '/')
}

`%notin%` <- Negate(`%in%`)

largeRels$total <- rowSums(largeRels)
smallRels$total <- rowSums(smallRels)

largeRels <- largeRels[order(largeRels$total, decreasing = TRUE), ]
smallRels <- smallRels[order(smallRels$total, decreasing = TRUE), ]

ppFrame <- merge(largeRels, smallRels, by = 0, all = TRUE)
ppFrame[is.na(ppFrame)] <- 0
rownames(ppFrame) <- ppFrame$Row.names
ppFrame$Row.names <- NULL

ppFrame$total <- rowSums(ppFrame)

ppAbund <- ppFrame[order(ppFrame$total, decreasing = TRUE), ]
ppAbund$total <- NULL

relabundTax <- function(countTable, cyaLevel, eukLevel){
  cyaSubTax <- cyanoTax[,c("seqID", cyaLevel)]
  eukSubTax <- eukTax[,c("seqID", eukLevel)]
  cyaSubTax[,"replaced"] <- gsub("\\s*\\([^\\)]+\\)","",as.character(cyaSubTax[,cyaLevel]))
  eukSubTax[,"replaced"] <- gsub("\\s*\\([^\\)]+\\)","",as.character(eukSubTax[,eukLevel]))
  cyaSubTax$Type <- "Cyano"
  eukSubTax$Type <- "Euk"
  cyaSubTax[cyaSubTax$replaced == "unclassified", "replaced"] <- 'Cyano_Unclassified'
  eukSubTax[grep("unclassified", eukSubTax$replaced), "replaced"] <- 'Euk_Unclassified'
  
  subTax <- rbind(cyaSubTax[, c("seqID", "replaced", "Type")], 
                  eukSubTax[, c("seqID", "replaced", "Type")])
  
  countTable$level <- subTax$replaced[match(rownames(countTable), subTax$seqID)]
  countTable$type <- subTax$Type[match(rownames(countTable), subTax$seqID)]
  
  #Get rid of the sequence names when aggregating, then sum up every sample's thing
  targs <- colnames(countTable)[colnames(countTable) %notin% 
                                  c("level", "total.x", "total.y", "total", "type")]
  levelCounts <- aggregate(countTable[targs], countTable["level"], FUN=sum)
  
  levelCounts$Total <- rowSums(levelCounts[,-grep("level", colnames(levelCounts))])
  levelCounts <- levelCounts[order(-levelCounts$Total),]
  
  rownames(levelCounts) <- levelCounts[,"level"]
  levelCounts$level <- NULL
  
  return(levelCounts)
}
lineClass <- relabundTax(ppAbund, "lineage", "class")
lineClass$Total <- NULL

cyaSubTax <- cyanoTax[,c("seqID", "lineage")]
eukSubTax <- eukTax[,c("seqID", "class")]
cyaSubTax[,"replaced"] <- gsub("\\s*\\([^\\)]+\\)","",as.character(cyaSubTax[,"lineage"]))
eukSubTax[,"replaced"] <- gsub("\\s*\\([^\\)]+\\)","",as.character(eukSubTax[,"class"]))
cyaSubTax$Type <- "Cyano"
eukSubTax$Type <- "Euk"
cyaSubTax[cyaSubTax$replaced == "unclassified", "replaced"] <- 'Cyano_Unclassified'
eukSubTax[grep("unclassified", eukSubTax$replaced), "replaced"] <- 'Euk_Unclassified'

subTax <- rbind(cyaSubTax[, c("seqID", "replaced", "Type")], 
                eukSubTax[, c("seqID", "replaced", "Type")])


targFrame <- lineClass
ibmColors <- c("#DC267F", "#785EF0", "#FE6100", "#648FFF", "#FFB000")
tolColors <- c("#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", 
               "#AA4499", "#882255", "#332288")
wongColors <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
                "#D55E00", "#CC79A7")
tolColorsPlus <- c("#117733", "#44AA99", "#88CCEE", "#DDCC77", "#CC6677", 
                   "#AA4499", "#882255", "#332288", "#6A3D9A")

depthTarg <- c("05m", "10m", "DCM", "50m", "90m")
sizeTargs <- c("02u", "3u")
allPoss <- paste(rep(depthTarg, each = length(sizeTargs)), sizeTargs, sep = "_")

twoParts1 <- strsplit(colnames(targFrame), "_")
subDates1 <- sapply(twoParts1, `[[`, 1)
subDates1 <- as.Date(subDates1, format = "%d%b%y")
targFrame <- targFrame[,order(subDates1)]
subDates1 <- sort(subDates1)

if(length(depthTarg)>1){
  depthStr <- paste0(depthTarg, sep = "|", collapse = "")
  depthStr <- substr(depthStr, 1, nchar(depthStr)-1)
  miniFrame <- 
    targFrame[,grep(depthStr, colnames(targFrame))]
} else{
  miniFrame <- targFrame[,grep(depthTarg, colnames(targFrame))]
}
miniFrame <- miniFrame[rowSums(miniFrame) > 0,]
miniMat <- as.matrix(miniFrame)

subOrPhyls <- rownames(miniFrame)
assignColors <- function(colors, phylDat){
  targetN <- length(colors)
  
  otherTab <- phylDat[grep("Other", rownames(phylDat)),]
  remDat <- phylDat[grep("Other", rownames(phylDat), invert = TRUE),]
  excessTab <- remDat[seq((targetN-1), nrow(remDat)-1),]
  types <-  subTax$Type[match(rownames(excessTab), subTax$replaced)]
  goodDat <- remDat[seq(1, targetN-2),]
  emptyOpt <- phylDat[grep("empty", rownames(phylDat)),]
  
  otherCyas <- colSums(excessTab[types == "Cyano",])+otherTab[rownames(otherTab)=="OtherCyano",]
  otherEuks <- colSums(excessTab[types == "Euk",])+otherTab[rownames(otherTab)=="OtherEuk",]
  
  miniPhyl <- rbind(goodDat, otherCyas, otherEuks, emptyOpt)
  
  sumsPerPhyl <- rowSums(miniPhyl)
  sumsPerPhyl[length(sumsPerPhyl)] <- 0
  miniPhyl <- miniPhyl[order(sumsPerPhyl, decreasing = TRUE),]
  
}
reorderMat <- function(matFram, datAll){
  rowSumFram <- rowSums(matFram)
  
  types <-  subTax$Type[match(rownames(matFram), subTax$replaced)]
  
  tempDF <- matFram[grep("_", rownames(matFram)),]
  tempTypes <- types[grep("_", rownames(matFram))]
  guilty <- which(is.na(tempTypes))
  if(length(guilty) > 0){
    for (i in 1:length(guilty)){
      subJect <- rownames(tempDF)[guilty[i]]
      subSplit <- strsplit(subJect, "_")
      typeOpt <- sapply(subSplit, `[[`, 1)
      
      tempTypes[guilty[i]] <- typeOpt
    }
    
    rownames(tempDF) <- tempTypes
  }
  
  
  matFram <- matFram[grep("_", rownames(matFram), invert = TRUE),]
  if(is.vector(tempDF[tempTypes == "Euk",])){
    matFram <- rbind(matFram, tempDF[tempTypes == "Euk",])
  }else{
    matFram <- rbind(matFram, colSums(tempDF[tempTypes == "Euk",]))
  }
  rownames(matFram)[nrow(matFram)] <- "OtherEuk"
  if(is.vector(tempDF[tempTypes == "Cyano",])){
    matFram <- rbind(matFram, tempDF[tempTypes == "Cyano",])
  }else{
    matFram <- rbind(matFram, colSums(tempDF[tempTypes == "Cyano",]))
  }
  rownames(matFram)[nrow(matFram)] <- "OtherCyano"
  
  sumsPerPhyl <- rowSums(matFram)
  matFram <- matFram[order(sumsPerPhyl, decreasing = TRUE),]
  
  matFram <- rbind(matFram, rep(0, ncol(matFram)))  
  rownames(matFram)[nrow(matFram)] <- "emptyCol"
  
  allCols <- paste(rep(format(datAll, "%d%b%y"), each = length(allPoss)), allPoss, sep = "_")
  allCols <- gsub("^0", "", allCols)
  
  for (i in seq(1, length(allCols))){
    targCol <- allCols[i]
    colPick <- which(colnames(matFram) == targCol)
    if (length(colPick) != 1){
      print(i)
      if (i == 1){
        matFram <- cbind("trash" = rep(0, nrow(matFram)), matFram)
        colnames(matFram)[grep("trash", colnames(matFram))] <- targCol
      }else if (i > ncol(matFram)){
        matFram <- cbind(matFram, "trash" = rep(0, nrow(matFram)))
        colnames(matFram)[grep("trash", colnames(matFram))] <- targCol
      }else{
        matFram <- cbind(matFram[,1:i-1], "trash" = rep(0, nrow(matFram)), 
                         matFram[,i:ncol(matFram)])
        colnames(matFram)[grep("trash", colnames(matFram))] <- targCol
      }
      matFram["emptyCol", targCol] <- 1
    }
  }
  return(matFram)
}

sortedMat <- reorderMat(miniMat, unique(subDates1))
colorSorting <- assignColors(tolColorsPlus, sortedMat)
sortedColors <- c(tolColorsPlus, "white")

idCs <- which(rownames(colorSorting) == "otherCyas")
cyanoBiums <- which(rownames(colorSorting) == "Cyanobiaceae")
colorSorting <- 
  colorSorting[c(cyanoBiums, idCs, setdiff(seq(1, nrow(colorSorting)), c(idCs, cyanoBiums))),]

rownames(colorSorting)[rownames(colorSorting) == "otherCyas"] <- "Other Cyanos"
rownames(colorSorting)[rownames(colorSorting) == "otherEuks"] <- "Other Ph. Euks"

colorSorting <- colorSorting[c(1,2,8,9,3,4,5,6,7,10),]

listOfPhyls <- rownames(colorSorting)





png("FigBin/allDepthsPPPlusBar.png", width = 2400, height = 1200)
plot.new()

par(new = "TRUE",plt = c(0.06,0.48,0.815,0.965),las = 1, cex.axis = 2.5)
x <- barplot(colorSorting[,grep("5m*_02u", colnames(colorSorting))], 
             col = sortedColors, yaxt = "n", xaxt = "n", space = 0.1)
axis(2, cex.axis=2, tck = -0.03, lwd.ticks = 3, hadj = 1.2, font = 2, lwd = 2.5, 
     labels = c(0, 20, 40, 60, 80, 100), at = c(0, 0.2, 0.4, 0.6, 0.8, 1))
mtext("Small (0.2-3 μm)", 3, cex = 2.5, font = 2, line = 1)
mtext("A", 3, cex = 2.5, font = 2, line = 1, adj = 0.05)
abline(v = 4.4, lty = 2, col = "black", lwd = 5)
abline(v = 21, lty = 2, col = "black", lwd = 5)
mtext("5 m", 2, cex = 2.5, font = 2, las = 0, line = 5)

par(new = "TRUE",plt = c(0.06,0.48,0.635,0.785),las = 1, cex.axis = 2.5)
x <- barplot(colorSorting[,grep("10m*_02u", colnames(colorSorting))], 
             col = sortedColors, yaxt = "n", xaxt = "n", space = 0.1)
axis(2, cex.axis=2, tck = -0.03, lwd.ticks = 3, hadj = 1.2, font = 2, lwd = 2.5, 
     labels = c(0, 20, 40, 60, 80, 100), at = c(0, 0.2, 0.4, 0.6, 0.8, 1))
abline(v = 4.4, lty = 2, col = "black", lwd = 5)
abline(v = 21, lty = 2, col = "black", lwd = 5)
mtext("10 m", 2, cex = 2.5, font = 2, las = 0, line = 5)

par(new = "TRUE",plt = c(0.06,0.48,0.455,0.605),las = 1, cex.axis = 2.5)
x <- barplot(colorSorting[,grep("DCM*_02u", colnames(colorSorting))], 
             col = sortedColors, yaxt = "n", xaxt = "n", space = 0.1)
axis(2, cex.axis=2, tck = -0.03, lwd.ticks = 3, hadj = 1.2, font = 2, lwd = 2.5, 
     labels = c(0, 20, 40, 60, 80, 100), at = c(0, 0.2, 0.4, 0.6, 0.8, 1))
abline(v = 4.4, lty = 2, col = "black", lwd = 5)
abline(v = 21, lty = 2, col = "black", lwd = 5)
mtext("% of Photosynthetic Organism 16S rRNA Genes in Sample", 2, cex = 2.5, font = 2, 
      las = 0, 
      line = 8, adj = 0.5)
mtext("Chl. Max", 2, cex = 2.5, font = 2, las = 0, line = 5)

par(new = "TRUE",plt = c(0.06,0.48,0.275,0.425),las = 1, cex.axis = 2.5)
x <- barplot(colorSorting[,grep("50m*_02u", colnames(colorSorting))], 
             col = sortedColors, yaxt = "n", xaxt = "n", space = 0.1)
axis(2, cex.axis=2, tck = -0.03, lwd.ticks = 3, hadj = 1.2, font = 2, lwd = 2.5, 
     labels = c(0, 20, 40, 60, 80, 100), at = c(0, 0.2, 0.4, 0.6, 0.8, 1))
abline(v = 4.4, lty = 2, col = "black", lwd = 5)
abline(v = 21, lty = 2, col = "black", lwd = 5)
mtext("50 m", 2, cex = 2.5, font = 2, las = 0, line = 5)

par(new = "TRUE",plt = c(0.06,0.48,0.095,0.245),las = 1, cex.axis = 2.5)
x <- barplot(colorSorting[,grep("90m*_02u", colnames(colorSorting))], 
             col = sortedColors, yaxt = "n", xaxt = "n", space = 0.1)
axis(2, cex.axis=2, tck = -0.03, lwd.ticks = 3, hadj = 1.2, font = 2, lwd = 2.5, 
     labels = c(0, 20, 40, 60, 80, 100), at = c(0, 0.2, 0.4, 0.6, 0.8, 1))
text(cex=1.5, x=x, y=-0.2, unique(format(subDates1, "%b %d")), xpd=TRUE, srt=90, font = 2)
mtext("2016", side = 1, line = 6, cex = 1.5, font = 2, adj = 0.075)
mtext("2017", side = 1, line = 6, cex = 1.5, font = 2, adj = 0.375)
mtext("2018", side = 1, line = 6, cex = 1.5, font = 2, adj = 0.775)
abline(v = 4.4, lty = 2, col = "black", lwd = 5)
abline(v = 21, lty = 2, col = "black", lwd = 5)
mtext("90 m", 2, cex = 2.5, font = 2, las = 0, line = 5)


par(new = "TRUE",plt = c(0.48,0.87,0.815,0.965),las = 1, cex.axis = 2.5)
barplot(colorSorting[,grep("05m*_3u", colnames(colorSorting))], 
        col = sortedColors, yaxt = "n" , xaxt = "n", space = 0.1)
axis(2, cex.axis=2, tck = -0.03, lwd.ticks = 3, hadj = 1.2, font = 2, lwd = 2.5, 
     labels = NA, at = c(0, 0.2, 0.4, 0.6, 0.8, 1))
mtext("Large (> 3μm)", 3, cex = 2.5, font = 2, line = 1)
mtext("B", 3, cex = 2.5, font = 2, line = 1, adj = 0.05)
abline(v = 4.4, lty = 2, col = "black", lwd = 5)
abline(v = 21, lty = 2, col = "black", lwd = 5)

legend("topright", inset=c(-0.3,0), legend = listOfPhyls[1:4], 
       fill = unique(sortedColors[1:4]), xpd = NA, cex = 1.8, 
       title = "Cyanobacteria Lineage", box.lwd = 2)

legend("bottomright", inset=c(-0.3,-0.9), legend = listOfPhyls[5:9], 
       fill = unique(sortedColors[5:9]), xpd = NA, cex = 1.8, title = "Ph. Plastid Class",
       box.lwd = 2)

par(new = "TRUE",plt = c(0.48,0.87,0.635,0.785),las = 1, cex.axis = 2.5)
x <- barplot(colorSorting[,grep("10m*_3u", colnames(colorSorting))], 
             col = sortedColors, yaxt = "n", xaxt = "n", space = 0.1)
axis(2, cex.axis=2, tck = -0.03, lwd.ticks = 3, hadj = 1.2, font = 2, lwd = 2.5, 
     labels = NA, at = c(0, 0.2, 0.4, 0.6, 0.8, 1))
abline(v = 4.4, lty = 2, col = "black", lwd = 5)
abline(v = 21, lty = 2, col = "black", lwd = 5)

par(new = "TRUE",plt = c(0.48,0.87,0.455,0.605),las = 1, cex.axis = 2.5)
x <- barplot(colorSorting[,grep("DCM*_3u", colnames(colorSorting))], 
             col = sortedColors, yaxt = "n", xaxt = "n", space = 0.1)
axis(2, cex.axis=2, tck = -0.03, lwd.ticks = 3, hadj = 1.2, font = 2, lwd = 2.5, 
     labels = NA, at = c(0, 0.2, 0.4, 0.6, 0.8, 1))
abline(v = 4.4, lty = 2, col = "black", lwd = 5)
abline(v = 21, lty = 2, col = "black", lwd = 5)

par(new = "TRUE",plt = c(0.48,0.87,0.275,0.425),las = 1, cex.axis = 2.5)
x <- barplot(colorSorting[,grep("50m*_3u", colnames(colorSorting))], 
             col = sortedColors, yaxt = "n", xaxt = "n", space = 0.1)
axis(2, cex.axis=2, tck = -0.03, lwd.ticks = 3, hadj = 1.2, font = 2, lwd = 2.5, 
     labels = NA, at = c(0, 0.2, 0.4, 0.6, 0.8, 1))
abline(v = 4.4, lty = 2, col = "black", lwd = 5)
abline(v = 21, lty = 2, col = "black", lwd = 5)

par(new = "TRUE",plt = c(0.48,0.87,0.095,0.245),las = 1, cex.axis = 2.5)
x <- barplot(colorSorting[,grep("90m*_3u", colnames(colorSorting))], 
             col = sortedColors, yaxt = "n", xaxt = "n", space = 0.1)
axis(2, cex.axis=2, tck = -0.03, lwd.ticks = 3, hadj = 1.2, font = 2, lwd = 2.5, 
     labels = NA, at = c(0, 0.2, 0.4, 0.6, 0.8, 1))
text(cex=1.5, x=x, y=-0.2, unique(format(subDates1, "%b %d")), xpd=TRUE, srt=90, font = 2)
mtext("2016", side = 1, line = 6, cex = 1.5, font = 2, adj = 0.075)
mtext("2017", side = 1, line = 6, cex = 1.5, font = 2, adj = 0.375)
mtext("2018", side = 1, line = 6, cex = 1.5, font = 2, adj = 0.775)
abline(v = 4.4, lty = 2, col = "black", lwd = 5)
abline(v = 21, lty = 2, col = "black", lwd = 5)



dev.off()
