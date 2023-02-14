smallAbund0 <- 
  read.csv("~/FlatheadMicrobes/SubsampledData/small_updatedNames.count_table", 
                        stringsAsFactors = FALSE, row.names = 1)
taxTable <- read.csv("RawData/prok.98.taxonomy")
`%notin%` <- Negate(`%in%`)

#Rename stuff
for (i in 1){
  colnames(smallAbund0) <- gsub("Sept", "Sep", colnames(smallAbund0))
  colnames(smallAbund0) <- gsub("June", "Jun", colnames(smallAbund0))
  colnames(smallAbund0) <- gsub("July", "Jul", colnames(smallAbund0))
  colnames(smallAbund0) <- sub("[^_]*_(.*)", "\\1", colnames(smallAbund0))
}

smallAbund <- smallAbund0
targDf <- smallAbund0

mainVal <- colSums(targDf[,grep("Sequence", colnames(targDf), invert = TRUE)])
targDf <- sweep(targDf[,grep("Sequence", colnames(targDf), invert = TRUE)], 
                2, mainVal, FUN = '/')
targDf$Sequence <- smallAbund0$Sequence


targDf$RelAbund <- rowSums(targDf[,colnames(targDf) != "Sequence"])/
  ((ncol(targDf)-1)*100)
targDf$MinVal <- apply(targDf[,colnames(targDf) != "Sequence"]/100, 1, min)
targDf$MaxVal <- apply(targDf[,colnames(targDf) != "Sequence"]/100, 1, max)
targDf$Sd <- apply(targDf[,colnames(targDf) != "Sequence"]/100, 1, sd)

reorderedSmall <- targDf[order(targDf$RelAbund, decreasing = TRUE),]

allTop10000 <- reorderedSmall
allTop10000[, c("Sd", "RelAbund", "MaxVal", "MinVal")] <- NULL

alldep <- c("05m", "10m", "DCM", "50m", "90m")

frame5m <- allTop10000[,grep("05m|Sequence", colnames(allTop10000))]
twoParts5 <- strsplit(colnames(frame5m), "_")
subDates5 <- sapply(twoParts5, `[[`, 1)
subDates5 <- as.Date(subDates5, format = "%d%b%y")
frame5m <- frame5m[, order(subDates5)]
subDates5 <- sort(subDates5)
frame10m <- allTop10000[,grep("10m|Sequence", colnames(allTop10000))]
twoParts10 <- strsplit(colnames(frame10m), "_")
subDates10 <- sapply(twoParts10, `[[`, 1)
subDates10 <- as.Date(subDates10, format = "%d%b%y")
frame10m <- frame10m[, order(subDates10)]
subDates10 <- sort(subDates10)
frameDCM <- allTop10000[,grep("DCM|Sequence", colnames(allTop10000))]
twoPartsDCM <- strsplit(colnames(frameDCM), "_")
subDatesDCM <- sapply(twoPartsDCM, `[[`, 1)
subDatesDCM <- as.Date(subDatesDCM, format = "%d%b%y")
frameDCM <- frameDCM[, order(subDatesDCM)]
subDatesDCM <- sort(subDatesDCM)
frame50m <- allTop10000[,grep("50m|Sequence", colnames(allTop10000))]
twoParts50 <- strsplit(colnames(frame50m), "_")
subDates50 <- sapply(twoParts50, `[[`, 1)
subDates50 <- as.Date(subDates50, format = "%d%b%y")
frame50m <- frame50m[, order(subDates50)]
subDates50 <- sort(subDates50)
frame90m <- allTop10000[,grep("90|Sequence", colnames(allTop10000))]
twoParts90 <- strsplit(colnames(frame90m), "_")
subDates90 <- sapply(twoParts90, `[[`, 1)
subDates90 <- as.Date(subDates90, format = "%d%b%y")
frame90m <- frame90m[, c(order(subDates90))]
subDates90 <- sort(subDates90)

frameDCM <- cbind(frameDCM[1:13], rep(0, nrow(frameDCM)), frameDCM[14:ncol(frameDCM)])
frameDCM <- cbind(frameDCM[1:31], rep(0, nrow(frameDCM)), rep(0, nrow(frameDCM)),
                  frameDCM[32])


targetPhylum <- c("Actinobacteria", "Bacteroidetes", "Cyanobacteria",
                  "Proteobacteria", "Chloroflexi", "Planctomycetes", "Verrucomicrobia")
colorsDF <- data.frame("Blue" = c("#0c2c84", "#225ea8", "#1d91c0", "#41b6c4",
                                  "#7fcdbb", "#c7e9b4"),
                       "Pink" = c("#91003f", "#ce1256", "#e7298a", "#df65b0",
                                  "#c994c7", "#d4b9da"),
                       "Green" = c("#005824", "#2ca25f", "#66c2a4", "#b2e2e2",
                                   "#99d8c9", "#ccece6"),
                       "Orange" = c("#8c2d04", "#cc4c02", "#ec7014", "#fe9929",
                                    "#fec44f", "#fee391"),
                       "Purple" = c("#6e016b", "#9ebcda", "#8c6bb1", "#8c96c6",
                                    "#9ebcda", "#bfd3e6"), 
                       "Red" = c("#b10026", "#e31a1c", "#fc4e2a", "#fd8d3c",
                                 "#feb24c", "#fed976"),
                       "Blue2" = c("#045a8d", "#2b8cbe", "#74a9cf", "#a6bddb",
                                 "#d0d1e6", "#f1eef6"))
ymaxes <- c(0.4, 0.1, 0.4, 0.141, 0.14, 0.1, 0.05)
legendY <- c(0.4, 0.075, 0.2, 0.11, 0.06, 0.1, 0.03)
axesRange <- matrix(c(seq(-0.4, 0.4, length.out = 5), seq(-0.1, 0.1, length.out = 5), 
                      seq(-0.4, 0.4, length.out = 5), seq(-0.14, 0.14, length.out = 5),
                      seq(-0.14, 0.14, length.out = 5), seq(-0.1, 0.1, length.out = 5),
                      seq(-0.05, 0.05, length.out = 5)), ncol = 7)
labelsOrgs <- list(c("acI-B1", "acI-A6", "acI-A7", "acI-B1", "acI-C2", "acIV-A"),
                   c("Fluviicola", "Sphingobacteriales", "bacI-A1", "Fluviicola", 
                     "bacI-A1", "Sphingobacteriales"),
                   c("Cyanobium", "Cyanobium", "Cyanobium", "Cyanobium"),
                   c("Alpha LD12", "Beta betI-A", "Beta Lhab-A1", "Beta LD28", 
                     "Beta Rhodoferax", "Beta betIV"),
                   c("Thermomarinilinea", "SL56"),
                   c("Gemmata", "Botrimarina", "Botrimarina", "Botrimarina",
                     "Fimbriiglobus", "Gemmataceae"),
                   c("Ereboglobus", "Pedosphaeraceae", "Chthoniobacter", 
                     "Opitutaceae", "LD19"))
nInt <- c(6, 6, 4, 6, 2, 6, 5)
phylTextX <- c(53, 53, 53, 53, 52, 53, 53)
dTop <- 0.01
diffLength <- -0.03
pltLength <- 0.88/(2*length(targetPhylum)) - diffLength/2
png("~/FlatheadMicrobes/FigBinExtras/sixPhyla_5_90_OTU_labels.png", width = 1250, 
    height = 1400)
j <- 1
plot.new()
for(i in 1:length(targetPhylum)){
  allPossibles <- taxTable[grep(targetPhylum[i], taxTable[,"phylum"]), "seqID"]
  
  newDf5m <- frame5m[rownames(frame5m) %in% allPossibles,]
  rownames(newDf5m) <- newDf5m$Sequence
  newDf5m$Sequence <- NULL
  newDf90m <- frame90m[rownames(frame90m) %in% allPossibles,]
  rownames(newDf90m) <- newDf90m$Sequence
  newDf90m$Sequence <- NULL
  
  rowMaxes5m <- apply(newDf5m, 1, FUN = max)
  rowMaxes90m <- apply(newDf90m, 1, FUN = max)
  
  rowMaxes <- pmax(rowMaxes5m, rowMaxes90m)
  newOrder <- order(rowMaxes, decreasing = TRUE)
  
  newDf5m <- newDf5m[newOrder,]
  newDf90m <- newDf90m[newOrder,]
  
  ymaxa <- 1 - dTop - pltLength*(j-1)*2 - diffLength * (j-1) #0.95
  ymina <- 1 - dTop-pltLength*(j*2-1)-diffLength*(j-1) #0.7625 1-dTop-pltLength
  ymaxb <- 1 - dTop-pltLength*(j*2-1)-diffLength*(j-1) #
  yminb <- 1 - dTop-pltLength*j*2-diffLength*(j-1) #0.575 1-0.05-2*pltLength
  
  
  par(new = "TRUE",plt = c(0.09,0.64,ymina,ymaxa),las = 1, xpd = NA)
  barplot(as.matrix(newDf5m[1:nInt[i],]), xaxt = "n", col = colorsDF[1:nInt[i],i], 
          border = NA, ylim = c(-ymaxes[i]/50, ymaxes[i]), yaxt = "n")
  if (i %% 2 == 1){
    axis(2, at = axesRange[,i], labels = paste0(abs(axesRange[,i]*100), "%"), 
         lwd = 3, cex.axis = 2, font.axis = 2)
    text(41.6, ymaxes[i]*0.15, "5 m", font = 2, cex = 2, col = "#636363")
    text(41.9, -ymaxes[i]*0.15, "90 m", font = 2, cex = 2, col = "#636363")
  }else{
    axis(4, at = axesRange[,i], labels = paste0(abs(axesRange[,i]*100), "%"), 
         lwd = 3, cex.axis = 2, font.axis = 2)
    text(-1.6, ymaxes[i]*0.15, "5 m", font = 2, cex = 2, col = "#636363")
    text(-1.9, -ymaxes[i]*0.15, "90 m", font = 2, cex = 2, col = "#636363")
  }
  if (i == 4){
    mtext("Relative Abundance (%)", 2, line = 6, las = 0, font = 2, cex = 2,
          adj = 1)
  }
  par(new = "TRUE",plt = c(0.09,0.64,yminb,ymaxb),las = 1)
  x <- barplot(-as.matrix(newDf90m[1:nInt[i],]), xaxt = "n", 
               col = colorsDF[1:nInt[i],i], 
               border = NA, ylim = c(-ymaxes[i], ymaxes[i]/50), yaxt = "n")
  legend(45, legendY[i], 
         legend = paste0(gsub("_", " ", rownames(newDf90m)[1:nInt[i]]), 
                         " (", labelsOrgs[[i]], ")"), 
         fill = colorsDF[1:nInt[i],i], bty = "n", cex = 2)
  text(phylTextX[i], legendY[i], targetPhylum[i], font = 2, cex = 2)
  j = j+1
}
colChoice <- "#494949"
xa0 <- 3.67
xa1 <- 10.85
xb0 <- 21.82
xb1 <- 28.9
y0 <- -0.05
y1 <- 0.54
lines(x = c(xa0, xa0), y = c(y0, y1), lty = 2, lwd = 3, 
      col = colChoice, xpd = NULL)
lines(x = c(xa1, xa1), y = c(y0, y1), lty = 2, lwd = 3, 
      col = colChoice, xpd = NULL)
lines(x = c(xa0, xa1), y = c(y0, y0), lty = 2, lwd = 3, 
      col = colChoice, xpd = NULL)
lines(x = c(xa0, xa1), y = c(y1, y1), lty = 2, lwd = 3, 
      col = colChoice, xpd = NULL)
lines(x = c(xb0, xb0), y = c(y0, y1), lty = 2, lwd = 3, 
      col = colChoice, xpd = NULL)
lines(x = c(xb1, xb1), y = c(y0, y1), lty = 2, lwd = 3, 
      col = colChoice, xpd = NULL)
lines(x = c(xb0, xb1), y = c(y0, y0), lty = 2, lwd = 3, 
      col = colChoice, xpd = NULL)
lines(x = c(xb0, xb1), y = c(y1, y1), lty = 2, lwd = 3, 
      col = colChoice, xpd = NULL)
text(cex=1.5, x=x, y=-0.07, unique(format(subDates90, "%b %d")), xpd=TRUE, srt=90, 
     font = 2)
text(2.5, -0.092, "2016", font = 2, cex = 2)
text(16, -0.092, "2017", font = 2, cex = 2)
text(34, -0.092, "2018", font = 2, cex = 2)
dev.off()









           #####Gathering stats
avgAbund <- rowMeans(allTop10000[,colnames(allTop10000) != "Sequence"])
avgAbund[1:10]
sum(avgAbund[1:10])
allPossibles <- taxTable[grep("Cyanobacteria", taxTable[,"phylum"]), "seqID"]
cyanoTax <- taxTable[grep("Cyanobacteria", taxTable[,"phylum"]),]
cyanoAbund <- allTop10000[rownames(allTop10000) %in% allPossibles,]
cyanoAbund$Avg <- rowMeans(cyanoAbund[, colnames(cyanoAbund) != "Sequence"])
cyanoAbund$seqID <- rownames(cyanoAbund)
cyanoFull <- merge(cyanoAbund, cyanoTax, by = "seqID", all.x = TRUE)
cyanoFull <- cyanoFull[order(cyanoFull$Avg, decreasing = TRUE),]
