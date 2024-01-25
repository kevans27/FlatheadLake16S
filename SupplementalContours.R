taxTable <- read.csv("RawData/prok.98.taxonomy")
countTableSmall <- read.csv("SubsampledData/small_updatedNames.count_table",
                            stringsAsFactors = FALSE, row.names = 1)
countTableSmall$blastID <- rownames(countTableSmall)
source("~/filled.contour3.R")

DCMDat <- read.csv("~/FlatheadPublic/DCMData.csv", header = TRUE, 
                   stringsAsFactors = FALSE)
DCMDat$Date <- as.Date(DCMDat$Date)
DCMDat$Date.1 <- as.Date(DCMDat$Date.1)
DCMDat$Date.1.1 <- as.Date(DCMDat$Date.1.1)
MLDDat <- read.csv("~/FlatheadPublic/MLDData.csv", header = TRUE, 
                   stringsAsFactors = FALSE)
MLDDat$Date <- as.Date(MLDDat$Date)
#MLDDat$Date.1 <- as.Date(MLDDat$Date.1)
#MLDDat$Date.1.1 <- as.Date(MLDDat$Date.1.1)

#Rename/trim things
for (i in 1){
  colnames(countTableSmall) <- sub("[^_]*_(.*)", "\\1", colnames(countTableSmall))
  colnames(countTableSmall) <- gsub("_[^_]+$", "\\1", colnames(countTableSmall))
  colnames(countTableSmall) <- gsub("Sept", "Sep", colnames(countTableSmall))
  colnames(countTableSmall) <- gsub("June", "Jun", colnames(countTableSmall))
}

`%notin%` <- Negate(`%in%`)
# Function to plot color bar
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), 
                      cexVal = 2.5, title='') {
  scale = (length(lut)-1)/(max-min)
  
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', 
       main=title, 
       cex.main = cexVal)
  axis(2, ticks, las=1, cex.axis = cexVal)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}

colMax <- function(data) sapply(data, max, na.rm = TRUE)
colSort <- function(data, ...) sapply(data, sort, ...)

#Print out sorted vector of abundances of a certain taxa
#-sort(-tribTable["LD12",])

##Contour plot of groups of interest

library(akima)
library(viridis)
taxMesh <- function(group, table){
  subTable <- table[table[, "Sequence"] == group,]
  
  subTable <- subTable[,!(colnames(subTable) %in% c("Total", "Sequence", "blastID"))]
  
  
  twoParts <- strsplit(colnames(subTable), "_")
  subDates <- sapply(twoParts, `[[`, 1)
  subDates <- as.Date(subDates, format = "%d%b%y")
  subDepths <- sapply(twoParts, `[[`, 2)
  
  for (j in which(subDepths == "DCM")){
    str <- colnames(subTable)[j]
    
    dat <- subDates[j]
    rowInd <- 
      which(c(DCMDat$Date, DCMDat$Date.1, DCMDat$Date.1.1) %in% dat) %% nrow(DCMDat)
    if (length(rowInd) == 1){
      if (rowInd == 0){
        rowInd <- nrow(DCMDat)
      }
      depth <- DCMDat[rowInd, "DCM"]
      
      subDepths[j] <- paste0(depth, "m")
    }
  }
  dcmsOut <- which(subDepths == "DCM")
  
  subDats <- subDates[-dcmsOut]
  subDeps <- subDepths[-dcmsOut]
  
  subDeps <- as.numeric(as.character(gsub("m", "", subDeps)))
  
  subValues <- as.numeric(subTable)
  subVals <- subValues[-dcmsOut]
  
  MLDDates <- MLDDat$Date
  MLDDepths <- MLDDat$MLD
  MLDDepths <- MLDDepths[order(MLDDates)]
  MLDDates <- MLDDates[order(MLDDates)]
  
  allDat <- list(subVals, subDats, subDeps, MLDDates, MLDDepths)
  
  return(allDat)
}

levelTableSmall <- countTableSmall
numberOfGroups <- 2
#groupsOfInterest <- c("M02585_71_000000000BT9TY_1_2110_18185_7085", 
#                      "M02585_71_000000000BT9TY_1_1104_13892_19400")
##Top 10 plus:
#sort(unique(c(36, 12, 11, 42, 20, 24, 26, 28, 36, 47, 28)))
#11, 12
#20, 24
#26, 28
#36, 42
#47
#
cyanoInterest <- c("Otu3", "Otu28", "Otu36", "Otu47")
cyanoNames <- c("Cyanobium", "Cyanobium", "Cyanobium", "Cyanobium")
goodNamesCyan <- gsub("tu", "TU", cyanoInterest)
epiInterest <- c("Otu1", "Otu5", "Otu6", "Otu10")
epiNames <- c("acI-B1", "acI-A6", "LD12", "Lhab-A1")
goodNamesEpi <- gsub("tu", "TU", epiInterest)
hypoInterest <- c("Otu4", "Otu11", "Otu12", "Otu42")
hypoNames <- c("Thermomarinilinea", "Gemmata", "LD28", "Fimbriiglobus")
goodNamesHypo <- gsub("tu", "TU", hypoInterest)
groupsOfInterest <- cyanoInterest
names <- cyanoNames

smallCex <- 0.9

##Check taxonomy
for (i in 1:length(groupsOfInterest)){
  groupOIblas <- countTableSmall[countTableSmall$Sequence == groupsOfInterest[i],
                                 "blastID"]
  taxonomy <- taxTable[taxTable$seqID == groupOIblas, ]
  print(c(groupsOfInterest[i], taxonomy))
}

library(png)
givenDates <- seq(as.Date("2016-09-01"), as.Date("2019-01-01"), 
                  by = "month")
thirdMonths <- givenDates[seq(1, length(givenDates), by = 3)]

makePlotOIMan <- function(groupsOI, levelTable, size, i){
  letter <- letters[i]
  actinos <- taxMesh(groupsOI, levelTable)
  actDates <- as.Date(actinos[[2]])
  actVals <- as.numeric(as.character(actinos[[1]]))/100
  actDepths <- as.numeric(as.character(actinos[[3]]))
  mldDates <- as.Date(actinos[[4]])
  mldDepths <- as.numeric(as.character(actinos[[5]]))
  
  fld2 <- interp(actDates, actDepths, actVals, duplicate = "mean", 
                 yo = seq(1, 90, 1), 
                 xo = seq(min(actDates), max(actDates), by="days"), extrap = FALSE)
  
  file <- paste0("./FigBinExtras/ClassDS/", groupsOI, "_Man.png")
  dir.create(dirname(file), showWarnings = FALSE)
  
  png(file, width = 2000, height = 1200)
  
  filled.contour(x = seq(min(actDates), max(actDates), by = "day"), y = fld2$y, 
                 z = fld2$z, color.palette = function(x)viridis(x), cex.main = 4,
                 nlevels = 15, cex = 4, cex.axis = 4, cex.lab = 4,
                 plot.axes={
                   if(i %in% c(1,2)){
                     axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1, 
                          font = 2, at = thirdMonths, labels = NA)
                     axis(1, cex.axis=2.5, labels = NA, 
                          at = givenDates, lwd.ticks = 3, tck = -0.02)
                   } else{
                     axis(1, cex.axis=4, lwd.ticks = 3, tck = -0.03, padj = 1.5, 
                          font = 2, at = thirdMonths, labels = NA)
                     axis(1, cex.axis=2.5, labels = NA, 
                          at = givenDates, lwd.ticks = 3, tck = -0.02)
                   }
                   lines(mldDates, mldDepths, col = "gray", lwd = 3);
                   abline(v=as.Date("2017-01-01"), lty = 2, col = "gray75", lwd = 3);
                   abline(v=as.Date("2018-01-01"), lty = 2, col = "gray75", lwd = 3);
                   abline(v=as.Date("2019-01-01"), lty = 2, col = "gray75", lwd = 3)
                   if(i %in% c(1,3)){
                     axis(2, cex.axis=4, lwd.ticks = 3, at = seq(20, 100, by = 20),
                          labels = NA, tck = -0.03)
                   } else{
                     axis(2, cex.axis=4, lwd.ticks = 3, tck = -0.03,
                          at = seq(20, 100, by = 20), labels = NA)
                   }
                   grid(col=rgb(1,1,1,0))
                   box(lwd = 4)
                 },
                 key.axes={
                   axis(4, cex.axis = 0.001) #####################
                 },
                 ylim = rev(range(fld2$y))
  )
  
  #mtext("% Abundance", 3, line = 1, font = 2, adj = 1, cex = 4)
  #mtext(letter, 3, line = 1, font = 2, adj = 0.03, cex = 5, padj = 1.65)
  
  dev.off()
}

levelTableSmall <- countTableSmall
numberOfGroups <- 4

letters <- c("a", "b", "c", "d")
makePlotOIMan(groupsOfInterest[1], countTableSmall, "small", 1)
makePlotOIMan(groupsOfInterest[2], countTableSmall, "small", 2)
makePlotOIMan(groupsOfInterest[3], countTableSmall, "small", 3)
makePlotOIMan(groupsOfInterest[4], countTableSmall, "small", 4)

library(png)
manualFil <- paste0("FigBinExtras/ClassDS/", groupsOfInterest[1], "_",
                    groupsOfInterest[2], "_PaneledManual.tiff")
tiff(manualFil, width = 7, height = 6, pointsize = 12, units = "in", res = 300)

par(mar = c(1, 1, 1, 0), xpd = NA)
plot(0:2, 0:2, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
rasterImage(readPNG(source=paste0("./FigBinExtras/ClassDS/", groupsOfInterest[1],
                                  "_Man.png")), -0.02, 1.08, 0.95, 2.1)
rasterImage(readPNG(source=paste0("./FigBinExtras/ClassDS/", groupsOfInterest[2],
                                  "_Man.png")), 1.02, 1.08, 1.99, 2.1)
rasterImage(readPNG(source=paste0("./FigBinExtras/ClassDS/", groupsOfInterest[3],
                                  "_Man.png")), -0.02, 0.02, 0.95, 1.04)
rasterImage(readPNG(source=paste0("./FigBinExtras/ClassDS/", groupsOfInterest[4],
                                  "_Man.png")), 1.02, 0.02, 1.99, 1.04)
mtext("Date", 1, line = -0.2, adj = 0.48)
depthSeq <- seq(20, 80, by = 20)
depthAdj <- seq(-20.1, -5.3, length.out = 4)
depthDiff <- 25.5
mtext(depthSeq, 2, line = -1, padj = depthAdj, las = 1, cex = smallCex)
mtext(depthSeq, 2, line = -1, padj = depthAdj+depthDiff, las = 1, cex = smallCex)
dateLabs <- format(thirdMonths, "%b")[seq(2, length(thirdMonths)-1)]
dateAdj <- seq(0.06, 0.41, length.out = length(dateLabs))
dateDiff <- 0.51
mtext(dateLabs, 1, line = -1.8, adj = dateAdj+dateDiff, cex = smallCex)
mtext(dateLabs, 1, line = -1.85, adj = dateAdj, cex = smallCex)
##Negative - raise
##Positive - lower
###-14.5, -12.7, -11
topLeft <- seq(25, 0, by = -5)
tlAdj <- seq(-22.3, -2.8, length.out = length(topLeft))
mtext(topLeft, 2, line = -16, padj = tlAdj, las = 1, adj = 0, cex = smallCex)
botLeft <- seq(5, 0, by = -1)
blAdj <- seq(2.5, 22.3, length.out = length(botLeft))
mtext(botLeft, 2, line = -16, padj = blAdj, las = 1, adj = 0, cex = smallCex)
topRigh <- seq(3, 0, by = -0.5)
trAdj <- seq(-22.1, -2.8, length.out = length(topRigh))
mtext(topRigh, 2, line = -32.3, padj = trAdj, las = 1, adj = 0, cex = smallCex)
botRigh <- seq(2.5, 0, by = -0.5)
brAdj <- seq(1.5, 22.3, length.out = length(botRigh))
mtext(botRigh, 2, line = -32.3, padj = brAdj, las = 1, adj = 0, cex = smallCex)
mtext("Depth (m)", 2, line = 0)
mtext(parse(text = paste0("italic('", names[1], "')~", "(", 
                          goodNamesCyan[1], ")")), 3, adj = 0.17, line = -0.3)
mtext(parse(text = paste0("italic('", names[2], "')~", "(", 
                          goodNamesCyan[2], ")")), 3, adj = 0.81, line = -0.3)
mtext(parse(text = paste0("italic('", names[3], "')~", "(", 
                          goodNamesCyan[3], ")")), 3, adj = 0.17, line = -14.1)
mtext(parse(text = paste0("italic('", names[4], "')~", "(", 
                          goodNamesCyan[4], ")")), 3, adj = 0.81, line = -14.1) 
text(0.04, 2.023, letters[1])
text(1.1, 2.023, letters[2])
text(0.04, 0.964, letters[3])
text(1.1, 0.964, letters[4])
mtext("2016", 1, line = -0.7, adj = 0.05, cex = smallCex)
mtext("2017", 1, line = -0.7, adj = 0.19, cex = smallCex)
mtext("2018", 1, line = -0.7, adj = 0.375, cex = smallCex)
mtext("2016", 1, line = -0.7, adj = 0.57, cex = smallCex)
mtext("2017", 1, line = -0.7, adj = 0.71, cex = smallCex)
mtext("2018", 1, line = -0.7, adj = 0.895, cex = smallCex)
dev.off()





groupsOfInterest <- epiInterest
names <- epiNames

letters <- c("a", "b", "c", "d")
makePlotOIMan(groupsOfInterest[1], countTableSmall, "small", 1)
makePlotOIMan(groupsOfInterest[2], countTableSmall, "small", 2)
makePlotOIMan(groupsOfInterest[3], countTableSmall, "small", 3)
makePlotOIMan(groupsOfInterest[4], countTableSmall, "small", 4)


manualFil <- paste0("FigBinExtras/ClassDS/", groupsOfInterest[1], "_",
                    groupsOfInterest[2], "_PaneledManual.tiff")
tiff(manualFil, width = 7, height = 6, pointsize = 12, units = "in", res = 300)

par(mar = c(1, 1, 1, 0), xpd = NA)
plot(0:2, 0:2, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
rasterImage(readPNG(source=paste0("./FigBinExtras/ClassDS/", groupsOfInterest[1],
                                  "_Man.png")), -0.02, 1.08, 0.95, 2.1)
rasterImage(readPNG(source=paste0("./FigBinExtras/ClassDS/", groupsOfInterest[2],
                                  "_Man.png")), 1.02, 1.08, 1.99, 2.1)
rasterImage(readPNG(source=paste0("./FigBinExtras/ClassDS/", groupsOfInterest[3],
                                  "_Man.png")), -0.02, 0.02, 0.95, 1.04)
rasterImage(readPNG(source=paste0("./FigBinExtras/ClassDS/", groupsOfInterest[4],
                                  "_Man.png")), 1.02, 0.02, 1.99, 1.04)
mtext("Date", 1, line = -0.2,adj = 0.48)
depthSeq <- seq(20, 80, by = 20)
depthAdj <- seq(-20.1, -5.3, length.out = 4)
depthDiff <- 25.5
mtext(depthSeq, 2, line = -1, padj = depthAdj, las = 1, cex = smallCex)
mtext(depthSeq, 2, line = -1, padj = depthAdj+depthDiff, las = 1, cex = smallCex)
dateLabs <- format(thirdMonths, "%b")[seq(2, length(thirdMonths)-1)]
dateAdj <- seq(0.06, 0.41, length.out = length(dateLabs))
dateDiff <- 0.51
mtext(dateLabs, 1, line = -1.8, adj = dateAdj+dateDiff, cex = smallCex)
mtext(dateLabs, 1, line = -1.85, adj = dateAdj, cex = smallCex)
##Negative - raise
##Positive - lower
###-14.5, -12.7, -11
topLeft <- seq(15, 0, by = -5)
tlAdj <- seq(-23.5, -2.8, length.out = length(topLeft))
mtext(topLeft, 2, line = -16, padj = tlAdj, las = 1, adj = 0, cex = smallCex)
botLeft <- seq(8, 0, by = -2)
blAdj <- seq(3, 22.3, length.out = length(botLeft))
mtext(botLeft, 2, line = -16, padj = blAdj, las = 1, adj = 0, cex = smallCex)
topRigh <- seq(12, 0, by = -2)
trAdj <- seq(-23, -2.8, length.out = length(topRigh))
mtext(topRigh, 2, line = -32.3, padj = trAdj, las = 1, adj = 0, cex = smallCex)
botRigh <- seq(4, 0, by = -1)
brAdj <- seq(0.8, 22.3, length.out = length(botRigh))
mtext(botRigh, 2, line = -32.3, padj = brAdj, las = 1, adj = 0, cex = smallCex)
mtext("Depth (m)", 2, line = 0)
mtext(paste0(names[1], " ", "(", goodNamesEpi[1], ")"), 3, adj = 0.17, line = -0.3)
mtext(paste0(names[2], " (", goodNamesEpi[2], ")"), 3, adj = 0.81, line = -0.3)
mtext(paste0(names[3], " (", goodNamesEpi[3], ")"), 3, adj = 0.17, line = -14.1)
mtext(paste0(names[4], " (", goodNamesEpi[4], ")"), 3, adj = 0.81, line = -14.1) 
text(0.04, 2.023, letters[1])
text(1.1, 2.023, letters[2])
text(0.04, 0.964, letters[3])
text(1.1, 0.964, letters[4])
mtext("2016", 1, line = -0.7, adj = 0.05, cex = smallCex)
mtext("2017", 1, line = -0.7, adj = 0.19, cex = smallCex)
mtext("2018", 1, line = -0.7, adj = 0.375, cex = smallCex)
mtext("2016", 1, line = -0.7, adj = 0.57, cex = smallCex)
mtext("2017", 1, line = -0.7, adj = 0.71, cex = smallCex)
mtext("2018", 1, line = -0.7, adj = 0.895, cex = smallCex)
dev.off()




groupsOfInterest <- hypoInterest
names <- hypoNames

letters <- c("a", "b", "c", "d")
makePlotOIMan(groupsOfInterest[1], countTableSmall, "small", 1)
makePlotOIMan(groupsOfInterest[2], countTableSmall, "small", 2)
makePlotOIMan(groupsOfInterest[3], countTableSmall, "small", 3)
makePlotOIMan(groupsOfInterest[4], countTableSmall, "small", 4)


manualFil <- paste0("FigBinExtras/ClassDS/", groupsOfInterest[1], "_",
                    groupsOfInterest[2], "_PaneledManual.tiff")
tiff(manualFil, width = 7, height = 6, pointsize = 12, units = "in", res = 300)

par(mar = c(1, 1, 1, 0), xpd = NA)
plot(0:2, 0:2, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
rasterImage(readPNG(source=paste0("./FigBinExtras/ClassDS/", groupsOfInterest[1],
                                  "_Man.png")), -0.02, 1.08, 0.95, 2.1)
rasterImage(readPNG(source=paste0("./FigBinExtras/ClassDS/", groupsOfInterest[2],
                                  "_Man.png")), 1.02, 1.08, 1.99, 2.1)
rasterImage(readPNG(source=paste0("./FigBinExtras/ClassDS/", groupsOfInterest[3],
                                  "_Man.png")), -0.02, 0.02, 0.95, 1.04)
rasterImage(readPNG(source=paste0("./FigBinExtras/ClassDS/", groupsOfInterest[4],
                                  "_Man.png")), 1.02, 0.02, 1.99, 1.04)
mtext("Date", 1, line = -0.2,adj = 0.48)
depthSeq <- seq(20, 80, by = 20)
depthAdj <- seq(-20.1, -5.3, length.out = 4)
depthDiff <- 25.5
mtext(depthSeq, 2, line = -1, padj = depthAdj, las = 1, cex = smallCex)
mtext(depthSeq, 2, line = -1, padj = depthAdj+depthDiff, las = 1, cex = smallCex)
dateLabs <- format(thirdMonths, "%b")[seq(2, length(thirdMonths)-1)]
dateAdj <- seq(0.06, 0.41, length.out = length(dateLabs))
dateDiff <- 0.51
mtext(dateLabs, 1, line = -1.8, adj = dateAdj+dateDiff, cex = smallCex)
mtext(dateLabs, 1, line = -1.85, adj = dateAdj, cex = smallCex)
##Negative - raise
##Positive - lower
###-14.5, -12.7, -11
topLeft <- seq(12, 0, by = -2)
tlAdj <- seq(-23.3, -2.8, length.out = length(topLeft))
mtext(topLeft, 2, line = -16, padj = tlAdj, las = 1, adj = 0, cex = smallCex)
botLeft <- seq(2.5, 0, by = -0.5)
blAdj <- seq(3, 22.3, length.out = length(botLeft))
mtext(botLeft, 2, line = -16, padj = blAdj, las = 1, adj = 0, cex = smallCex)
topRigh <- seq(5, 0, by = -1)
trAdj <- seq(-25, -2.8, length.out = length(topRigh))
mtext(topRigh, 2, line = -32.3, padj = trAdj, las = 1, adj = 0, cex = smallCex)
botRigh <- seq(2, 0, by = -0.5)
brAdj <- seq(4.1, 22.3, length.out = length(botRigh))
mtext(botRigh, 2, line = -32.3, padj = brAdj, las = 1, adj = 0, cex = smallCex)
mtext("Depth (m)", 2, line = 0)
mtext(parse(text = paste0("italic('", names[1], "')~", "(", 
                          goodNamesHypo[1], ")")), 3, adj = 0.17, line = -0.3)
mtext(parse(text = paste0("italic('", names[2], "')~", "(", 
                          goodNamesHypo[2], ")")), 3, adj = 0.81, line = -0.3)
mtext(paste0(names[3], " (", goodNamesHypo[3], ")"), 3, adj = 0.17, line = -14.1)
mtext(parse(text = paste0("italic('", names[4], "')~", "(", 
                          goodNamesHypo[4], ")")), 3, adj = 0.81, line = -14.1) 
text(0.04, 2.023, letters[1])
text(1.1, 2.023, letters[2])
text(0.04, 0.964, letters[3])
text(1.1, 0.964, letters[4])
mtext("2016", 1, line = -0.7, adj = 0.05, cex = smallCex)
mtext("2017", 1, line = -0.7, adj = 0.19, cex = smallCex)
mtext("2018", 1, line = -0.7, adj = 0.375, cex = smallCex)
mtext("2016", 1, line = -0.7, adj = 0.57, cex = smallCex)
mtext("2017", 1, line = -0.7, adj = 0.71, cex = smallCex)
mtext("2018", 1, line = -0.7, adj = 0.895, cex = smallCex)
dev.off()