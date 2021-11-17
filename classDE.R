library(DESeq2)

countTableLarge <- read.csv("~/FlatheadMicrobes/SubsampledData/rrBothLarge.count_table", 
                            stringsAsFactors = FALSE, row.names = 1)
countTableLarge$Representative_Sequence <- NULL
countTableSmall <- read.csv("~/FlatheadMicrobes/SubsampledData/rrBothSmall.count_table", 
                            stringsAsFactors = FALSE, row.names = 1)
countTableSmall$Representative_Sequence <- NULL



#countTable <- merge(countTableLarge, countTableSmall, by = "row.names", all = TRUE)

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
}

colDataLarge <- data.frame("Depth" = rep(NA, ncol(countTableLarge)), "SurfTemp" = NA,
                           "Date" = as.Date(NA))
colDataSmall <- data.frame("Depth" = rep(NA, ncol(countTableSmall)), "SurfTemp" = NA, 
                           "Date" = as.Date(NA))

for (i in 1){
  colnames(countTableLarge) <- gsub("Sept", "Sep", colnames(countTableLarge))
  colnames(countTableLarge) <- gsub("June", "Jun", colnames(countTableLarge))
  colnames(countTableLarge) <- gsub("July", "Jul", colnames(countTableLarge))
  colnames(countTableSmall) <- gsub("Sept", "Sep", colnames(countTableSmall))
  colnames(countTableSmall) <- gsub("June", "Jun", colnames(countTableSmall))
  colnames(countTableSmall) <- gsub("July", "Jul", colnames(countTableSmall))
  
  twoPartsL <- strsplit(colnames(countTableLarge), "_")
  subDatesL <- sapply(twoPartsL, `[[`, 2)
  subDatesL <- as.Date(subDatesL, format = "%d%b%y")
  depthsL <- sapply(twoPartsL, `[[`, 3)
  
  twoPartsS <- strsplit(colnames(countTableSmall), "_")
  subDatesS <- sapply(twoPartsS, `[[`, 2)
  subDatesS <- as.Date(subDatesS, format = "%d%b%y")
  depthsS <- sapply(twoPartsS, `[[`, 3)
  
  colDataLarge$Depth <- depthsL
  colDataSmall$Depth <- depthsS
  
  colDataLarge$SurfTemp <- tempDF[match(subDatesL, tempDF$Date), "Temp"]
  colDataSmall$SurfTemp <- tempDF[match(subDatesS, tempDF$Date), "Temp"]
  
  colDataLarge$Date <- subDatesL
  colDataSmall$Date <- subDatesS
  
  colDataLarge[colDataLarge$Depth == "DCM", "Depth"] <- 
    DCMDat[match(colDataLarge[colDataLarge$Depth == "DCM", "Date"], DCMDat$Date), "DCM"]
  colDataSmall[colDataSmall$Depth == "DCM", "Depth"] <- 
    DCMDat[match(colDataSmall[colDataSmall$Depth == "DCM", "Date"], DCMDat$Date), "DCM"]
  
  colDataLarge$Depth <- gsub("m", "", colDataLarge$Depth)
  colDataLarge$Depth <- as.numeric(as.character(colDataLarge$Depth))
  colDataSmall$Depth <- gsub("m", "", colDataSmall$Depth)
  colDataSmall$Depth <- as.numeric(as.character(colDataSmall$Depth))
  
  removeLarge <- which(rowSums(is.na(colDataLarge)) > 0)
  colDataLarge <- colDataLarge[-removeLarge,]
  countTableLarge <- countTableLarge[, -removeLarge]
  removeSmall <- which(rowSums(is.na(colDataSmall)) > 0)
  colDataSmall <- colDataSmall[-removeSmall,]
  countTableSmall <- countTableSmall[, -removeSmall]
}

taxTable <- read.csv("RawData/prok.98.taxonomy")
taxTable[,"tclass"] <- gsub("\\s*\\([^\\)]+\\)","",as.character(taxTable$class))

for (i in unique(taxTable$tclass)){
  targs <- taxTable[taxTable$tclass == i, "seqID"]
  largeDF <- countTableLarge[rownames(countTableLarge) %in% targs,]
  smallDF <- countTableSmall[rownames(countTableSmall) %in% targs,]
  
  classVecLarge <- colSums(largeDF)
  classVecSmall <- colSums(smallDF)
  
  if(exists("largeClass")){
    largeClass <- rbind(largeClass, classVecLarge)
  } else{
    largeClass <- classVecLarge
  }
  
  if(exists("smallClass")){
    smallClass <- rbind(smallClass, classVecSmall)
  } else{
    smallClass <- classVecSmall
  }
}

rownames(largeClass) <- unique(taxTable$tclass)
rownames(smallClass) <- unique(taxTable$tclass)

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

for (i in unique(eukTax$tclass)){
  targs <- eukTax[eukTax$tclass == i, "seqID"]
  largeDF <- countTableLarge[rownames(countTableLarge) %in% targs,]
  smallDF <- countTableSmall[rownames(countTableSmall) %in% targs,]
  
  classVecLarge <- colSums(largeDF)
  classVecSmall <- colSums(smallDF)
  
  if(exists("largeEuk")){
    largeEuk <- rbind(largeEuk, classVecLarge)
  } else{
    largeEuk <- classVecLarge
  }
  
  if(exists("smallEuk")){
    smallEuk <- rbind(smallEuk, classVecSmall)
  } else{
    smallEuk <- classVecSmall
  }
}

rownames(largeEuk) <- unique(eukTax$tclass)
rownames(smallEuk) <- unique(eukTax$tclass)


largeClass <- rbind(largeClass, largeEuk)
smallClass <- rbind(smallClass, smallEuk)

largeSums <- rowSums(largeClass)
largeClass <- largeClass[order(largeSums, decreasing = TRUE),]
smallSums <- rowSums(smallClass)
smallClass <- smallClass[order(smallSums, decreasing = TRUE),]

deseqDSLargeClassST <- DESeqDataSetFromMatrix(largeClass, colDataLarge, 
                                              ~ Depth + SurfTemp)
deseqDSSmallClassST <- DESeqDataSetFromMatrix(smallClass, colDataSmall, 
                                              ~ Depth + SurfTemp)

deLargeOutputClassST <- DESeq(deseqDSLargeClassST)
deSmallOutputClassST <- DESeq(deseqDSSmallClassST)

deseqDSLargeClassD <- DESeqDataSetFromMatrix(largeClass, colDataLarge, 
                                             ~ SurfTemp + Depth)
deseqDSSmallClassD <- DESeqDataSetFromMatrix(smallClass, colDataSmall, 
                                             ~ SurfTemp + Depth)

deLargeOutputClassD <- DESeq(deseqDSLargeClassD)
deSmallOutputClassD <- DESeq(deseqDSSmallClassD)

classLargeTop10Depth <- results(deLargeOutputClassD)[grep("unclassified", 
                                                         rownames(deLargeOutputClassD),
                                                         invert = TRUE)[1:10],]
classLargeTop10ST <- results(deLargeOutputClassST)[grep("unclassified", 
                                                       rownames(deLargeOutputClassST),
                                                       invert = TRUE)[1:10],]
classSmallTop10Depth <- results(deSmallOutputClassD)[grep("unclassified", 
                                                         rownames(deSmallOutputClassD),
                                                         invert = TRUE)[1:10],]
classSmallTop10ST <- results(deSmallOutputClassST)[grep("unclassified", 
                                                       rownames(deSmallOutputClassST),
                                                       invert = TRUE)[1:10],]

plotReadyLarge <- data.frame("Organism" = rownames(classLargeTop10ST), 
                             "STDE" = classLargeTop10ST$log2FoldChange, 
                             "DepthDE" = classLargeTop10Depth$log2FoldChange,
                             "AvgAbundance" = classLargeTop10Depth$baseMean,
                             "OrgColor" = "#1A85FF")
plotReadyLarge$OrgColor <- as.character(plotReadyLarge$OrgColor)
plotReadyLarge$Organism <- as.character(plotReadyLarge$Organism)
plotReadyLarge[plotReadyLarge$Organism %in% eukTax$tclass, "OrgColor"] = "#D41159"
plotReadySmall <- data.frame("Organism" = rownames(classSmallTop10ST), 
                             "STDE" = classSmallTop10ST$log2FoldChange, 
                             "DepthDE" = classSmallTop10Depth$log2FoldChange,
                             "AvgAbundance" = classSmallTop10Depth$baseMean,
                             "OrgColor" = "#1A85FF")
plotReadySmall$OrgColor <- as.character(plotReadySmall$OrgColor)
plotReadySmall$Organism <- as.character(plotReadySmall$Organism)
plotReadySmall[plotReadySmall$Organism %in% eukTax$tclass, "OrgColor"] <- "#D41159"

#Temp range is 1 to 21.53 (20.53)
#Depth range is 5 to 90 (85)
#Multply axes by those ranges
#
#Small y range is about -2.3 to 1, large is about -2 to 0.7
#Small x range is about -1 to 1, large is -1 to 2.5

png("FigBinExtras/ClassDS/depthSTDE.png", width = 1800, height = 945)
plot.new()
par(new = "TRUE",plt = c(0.06,0.44,0.1,0.9),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(plotReadySmall$STDE*20.53, plotReadySmall$DepthDE*85, xlab = "", ylab = "",
     cex = plotReadySmall$AvgAbundance/100, lwd = 2, ylim = c(2.5, -2.5), 
     xlim = c(-2.5, 2.5), col = plotReadySmall$OrgColor,
     xaxt = "n", yaxt = "n")
abline(h = 0, lty = 2, lwd = 2.5, col = "gray")
abline(v = 0, lty = 2, lwd = 2.5, col = "gray")
x_vec <- plotReadySmall$STDE*20.53
x_vec[plotReadySmall$Organism %in% c("Alphaproteobacteria")] <-
  x_vec[plotReadySmall$Organism %in% c("Alphaproteobacteria")]-2
x_vec[plotReadySmall$Organism %in% c("Cryptophyceae")] <-
  x_vec[plotReadySmall$Organism %in% c("Cryptophyceae")]-1.8
x_vec[plotReadySmall$Organism %in% c("Planctomycetacia", "Oxyphotobacteria")] <-
  x_vec[plotReadySmall$Organism %in% c("Planctomycetacia", "Oxyphotobacteria")]-0.75
x_vec[plotReadySmall$Organism %in% c("Bacteroidia")] <-
  x_vec[plotReadySmall$Organism %in% c("Bacteroidia")]-0.55
x_vec[plotReadySmall$Organism %in% c("Verrucomicrobiae")] <-
  x_vec[plotReadySmall$Organism %in% c("Verrucomicrobiae")]+0.7
y_vec <- plotReadySmall$DepthDE*85
y_vec[plotReadySmall$Organism %in% c("Gammaproteobacteria")] <-
  y_vec[plotReadySmall$Organism %in% c("Gammaproteobacteria")]+0.8
y_vec[plotReadySmall$Organism %in% c("Bacteroidia")] <-
  y_vec[plotReadySmall$Organism %in% c("Bacteroidia")]+0.1
y_vec[plotReadySmall$Organism %in% c("Phycisphaerae")] <-
  y_vec[plotReadySmall$Organism %in% c("Phycisphaerae")]-0.1
y_vec[plotReadySmall$Organism %in% c("Acidimicrobiia")] <-
  y_vec[plotReadySmall$Organism %in% c("Acidimicrobiia")]-0.25
y_vec[plotReadySmall$Organism %in% c("Actinobacteria")] <-
  y_vec[plotReadySmall$Organism %in% c("Actinobacteria")]-1.05
text(x_vec, y_vec, labels=plotReadySmall$Organism,
     cex= 1.5)
arrows(0.68, 0.2, 0.68, -0.3, lwd = 2, length = 0.15)
arrows(-0.55, -1.15, 0.62, -1.15, lwd = 2, length = 0.15)
arrows(-0.5, -0.78, 0.75, -0.78, lwd = 2, length = 0.15)
box(lwd = 2.5)
mtext(bquote(bold("log2-Fold Change (Surface Temperature)")), 1, font = 2, line = 5, 
      cex = 2.5, las = 0, adj = 4.5)
mtext(bquote(bold("log2-Fold Change (Depth)")), 2, font = 2, line = 5, 
      cex = 2.5, las = 0)
mtext("Small (0.2-3 μm)", 3, font = 2, line = 2, cex = 2.5)
mtext("A", 3, font = 2, line = 2, cex = 2.5, adj = 0.05)
axis(1, cex.axis=2.5, tck = -0.02, lwd.ticks = 2.5, padj = 1, font = 2, 
     at = seq(-2, 2, by = 2), labels = seq(-2, 2, by = 2))
axis(1, cex.axis=2.5, tck = -0.015, lwd.ticks = 2.5, padj = 1.2, font = 2, 
     at = seq(-1, 3, by = 2), labels = NA)
axis(2, cex.axis=2.5, tck = -0.02, lwd.ticks = 2.5, hadj = 1.2, font = 2, 
     at = seq(-2, 2, by = 2), labels = seq(2, -2, by = -2))
axis(2, cex.axis=2.5, tck = -0.015, lwd.ticks = 2.5, hadj = 1.2, font = 2, 
     at = seq(-1, 3, by = 2), labels = NA)
# arrows(1.1, -0.38, 0.8, -0.45, lwd = 2, length = 0.15)
# arrows(1.18, -1.2, 0.88, -1.12, lwd = 2, length = 0.15)
text(0, -2.5, "Shallower", cex= 3, col = "gray30", font = 2)
text(0, 2.5, "Deeper", cex= 3, col = "gray30", font = 2)
text(-2.2, 0, "Cooler", cex= 3, col = "gray30", font = 2)
text(2.1, 0, "Warmer", cex= 3, col = "gray30", font = 2)


par(new = "TRUE",plt = c(0.49,0.87,0.1,0.9),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(plotReadyLarge$STDE*20.53, plotReadyLarge$DepthDE*85, xlab = "", ylab = "",
     cex = plotReadyLarge$AvgAbundance/100, lwd = 2, ylim = c(2.5, -2.5), 
     xlim = c(-2.5, 2.5), col = plotReadyLarge$OrgColor,
     xaxt = "n", yaxt = "n")
abline(h = 0, lty = 2, lwd = 2.5, col = "gray")
abline(v = 0, lty = 2, lwd = 2.5, col = "gray")
x_vec <- plotReadyLarge$STDE*20.53
x_vec[plotReadyLarge$Organism %in% c("Oxyphotobacteria")] <-
  x_vec[plotReadyLarge$Organism %in% c("Oxyphotobacteria")]-1.1
x_vec[plotReadyLarge$Organism %in% c("Planctomycetacia")] <-
  x_vec[plotReadyLarge$Organism %in% c("Planctomycetacia")]-0.85
x_vec[plotReadyLarge$Organism %in% c("Bacteroidia")] <-
  x_vec[plotReadyLarge$Organism %in% c("Bacteroidia")]-0.65
x_vec[plotReadyLarge$Organism %in% c("Actinobacteria")] <-
  x_vec[plotReadyLarge$Organism %in% c("Actinobacteria")]-0.18
x_vec[plotReadyLarge$Organism %in% c("Phycisphaerae")] <-
  x_vec[plotReadyLarge$Organism %in% c("Phycisphaerae")]+0.05
x_vec[plotReadyLarge$Organism %in% c("Verrucomicrobiae")] <-
  x_vec[plotReadyLarge$Organism %in% c("Verrucomicrobiae")]+0.7
x_vec[plotReadyLarge$Organism %in% c("Alphaproteobacteria")] <-
  x_vec[plotReadyLarge$Organism %in% c("Alphaproteobacteria")]+0.72
x_vec[plotReadyLarge$Organism %in% c("Gammaproteobacteria")] <-
  x_vec[plotReadyLarge$Organism %in% c("Gammaproteobacteria")]+0.85
y_vec <- plotReadyLarge$DepthDE*85
y_vec[plotReadyLarge$Organism %in% c("Alphaproteobacteria")] <-
  y_vec[plotReadyLarge$Organism %in% c("Alphaproteobacteria")]-0.12
y_vec[plotReadyLarge$Organism %in% c("Bacillariophyta")] <-
  y_vec[plotReadyLarge$Organism %in% c("Bacillariophyta")]-0.38
y_vec[plotReadyLarge$Organism %in% c("Phycisphaerae", "Chrysophyceae", "Actinobacteria")] <-
  y_vec[plotReadyLarge$Organism %in% 
          c("Phycisphaerae", "Chrysophyceae", "Actinobacteria")]-0.2
text(x_vec, y_vec, labels=plotReadyLarge$Organism,
     cex= 1.5)
mtext("Large (>3 μm)", 3, font = 2, line = 2, cex = 2.5)
mtext("B", 3, font = 2, line = 2, cex = 2.5, adj = 0.05)
box(lwd = 2.5)
axis(1, cex.axis=2.5, tck = -0.02, lwd.ticks = 2.5, padj = 1, font = 2, 
     at = seq(-2, 2, by = 2), labels = seq(-2, 2, by = 2))
axis(1, cex.axis=2.5, tck = -0.015, lwd.ticks = 2.5, padj = 1.2, font = 2, 
     at = seq(-1, 3, by = 2), labels = NA)
axis(2, cex.axis=2.5, tck = -0.02, lwd.ticks = 2.5, hadj = 1.2, font = 2, 
     at = seq(-2, 2, by = 2), labels = seq(2, -2, by = -2))
axis(2, cex.axis=2.5, tck = -0.015, lwd.ticks = 2.5, hadj = 1.2, font = 2, 
     at = seq(-1, 3, by = 2), labels = NA)
text(0, -2.5, "Shallower", cex= 3, col = "gray30", font = 2)
text(0, 2.5, "Deeper", cex= 3, col = "gray30", font = 2)
text(-2.2, 0, "Cooler", cex= 3, col = "gray30", font = 2)
text(2.1, 0, "Warmer", cex= 3, col = "gray30", font = 2)

legend("topright", legend = c("Prokaryotes", "Ph. Euks"), fill = c("#1A85FF", "#D41159"), 
       inset = c(-0.32, 0), xpd = TRUE, box.lwd = 2.5, cex = 2)
par(xpd = TRUE)
xPos <- 3.6
points(x = rep(xPos, 3), y = c(-0.8, 0.2, 1), pch = 1, lwd = 2, cex = c(1200, 600, 120)/100)
text(xPos, -1.8, "Average %", cex = 2, font = 2)
text(xPos, -1.6, "Abundance", cex = 2, font = 2)
text(xPos, -1.2, "10%", cex = 2)
text(xPos, -0.05, "5%", cex = 2)
text(xPos, 0.8, "1%", cex = 2)

dev.off()

