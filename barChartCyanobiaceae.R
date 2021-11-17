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
  
  largeAbund$Representative_Sequence <- NULL
  smallAbund$Representative_Sequence <- NULL
}

proktax <- read.csv("RawData/prok.98.taxonomy")
applyTaxonomy <- function(df){
  firstStepTax <- as.character(cyanoTax[match(rownames(df), cyanoTax$seqID), "clade"])
  problemChildren <- which(firstStepTax == "unclassified")
  for (i in problemChildren){
    firstStepTax[i] <- 
      as.character(cyanoTax[match(df$tip.label[i], cyanoTax$seqID),
                            as.numeric(apply(cyanoTax[match(df$tip.label[i], 
                                                            cyanoTax$seqID),],1, 
                                             match,x="unclassified")-1)])
  }
  firstStepTax <- gsub("\\s*\\([^\\)]+\\)","",firstStepTax)
  firstStepTax[problemChildren] <- paste0("unclassified (", firstStepTax[problemChildren], ")")
  
  df$Name <- firstStepTax
  return(firstStepTax)
}
`%notin%` <- Negate(`%in%`)

targetLine <- "Cyanobiaceae"
seqTargs <- as.character(proktax[grep(targetLine, proktax$lineage), "seqID"])

targetFrameLarge <- largeAbund[rownames(largeAbund) %in% seqTargs,]
targetFrameSmall <- smallAbund[rownames(smallAbund) %in% seqTargs,]

cyanobiumFrame <- merge(targetFrameLarge, targetFrameSmall, by = 0, all = TRUE)
cyanobiumFrame[is.na(cyanobiumFrame)] <- 0
rownames(cyanobiumFrame) <- cyanobiumFrame$Row.names
cyanobiumFrame$Row.names <- NULL
cyanobiumFrame$Sequence.x <- NULL
cyanobiumFrame$Sequence.y <- NULL

mainVal <- colSums(cyanobiumFrame[, grep("total", colnames(cyanobiumFrame), invert = TRUE)])
trimmedFram <- cyanobiumFrame[, grep("total", colnames(cyanobiumFrame), invert = TRUE)]
miniMat <- sweep(trimmedFram, 2, mainVal, FUN = '/')
miniMat$total <- rowSums(miniMat)

cyanobiumAbund <- miniMat[order(miniMat$total, decreasing = TRUE), ]
cyanobiumAbund$total <- NULL
prokKey <- read.csv("RawData/prok.name.key")

rownames(cyanobiumAbund) <- prokKey[match(rownames(cyanobiumAbund), prokKey$blastname),
                                    "qseqids"]

targFrame <- cyanobiumAbund

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
  
  unclassifiedDat <- phylDat[grep("Unclassified", rownames(phylDat)),]
  phylDat <- phylDat[grep("Unclassified", rownames(phylDat), invert = TRUE),]
  
  excessTab <- phylDat[seq((targetN), nrow(phylDat)-1),]
  goodDat <- phylDat[seq(1, targetN-1),]
  emptyOpt <- phylDat[grep("empty", rownames(phylDat)),]
  
  extraextra <- rbind(unclassifiedDat, excessTab)
  
  otherProks <- colSums(extraextra)
  
  miniPhyl <- rbind(goodDat, otherProks, emptyOpt)
  
  sumsPerPhyl <- rowSums(miniPhyl)
  sumsPerPhyl[length(sumsPerPhyl)] <- 0
  miniPhyl <- miniPhyl[order(sumsPerPhyl, decreasing = TRUE),]
  
}
reorderMat <- function(matFram, datAll){
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
colorSorting <- assignColors(tolColors, sortedMat)
sortedColors <- c(tolColors, "white")

rownames(colorSorting)[rownames(colorSorting) == "otherProks"] <- "Other Cyanobiaceae"

listOfPhyls <- rownames(colorSorting)





png("FigBin/allDepthsBarChartCyanobiaceaeTribePlus.png", width = 2400, height = 1200)
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
mtext("% of Cyanobiaceae 16S rRNA Genes in Sample", 2, cex = 2.5, font = 2, las = 0, 
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


legend("topright", inset=c(-0.3,0), legend = listOfPhyls[1:8], 
       fill = unique(sortedColors[1:8]), xpd = NA, cex = 1.8, 
       title = "Cyanobiaceae OTU",
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
