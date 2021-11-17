taxTable <- read.csv("RawData/prok.98.taxonomy")
countTableLarge <- read.csv("~/FlatheadMicrobes/SubsampledData/rrBothLarge.count_table", 
                            stringsAsFactors = FALSE, row.names = 1)
countTableLarge$Representative_Sequence <- NULL
countTableSmall <- read.csv("~/FlatheadMicrobes/SubsampledData/rrBothSmall.count_table", 
                            stringsAsFactors = FALSE, row.names = 1)
countTableSmall$Representative_Sequence <- NULL
transKey <- read.csv("RawData/prok.name.key")
source("~/filled.contour3.R")

DCMDat <- read.csv("~/FlatheadPublic/DCMData.csv", header = TRUE, stringsAsFactors = FALSE)
DCMDat$Date <- as.Date(DCMDat$Date)
DCMDat$Date.1 <- as.Date(DCMDat$Date.1)
DCMDat$Date.1.1 <- as.Date(DCMDat$Date.1.1)
MLDDat <- read.csv("~/FlatheadPublic/MLDData.csv", header = TRUE, stringsAsFactors = FALSE)
MLDDat$Date <- as.Date(MLDDat$Date)
MLDDat$Date.1 <- as.Date(MLDDat$Date.1)
MLDDat$Date.1.1 <- as.Date(MLDDat$Date.1.1)

#Rename/trim things
for (i in 1){
  countTableLarge$blastname <- rownames(countTableLarge)
  colnames(countTableLarge) <- sub("[^_]*_(.*)", "\\1", colnames(countTableLarge))
  colnames(countTableLarge) <- gsub("_[^_]+$", "\\1", colnames(countTableLarge))
  colnames(countTableLarge) <- gsub("Sept", "Sep", colnames(countTableLarge))
  colnames(countTableLarge) <- gsub("June", "Jun", colnames(countTableLarge))
  
  countTableSmall$blastname <- rownames(countTableSmall)
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
relabundTax <- function(counts, taxes, level){
  subTax <- taxes[,c("seqID", level)]
  subTax[,"replaced"] <- gsub("\\s*\\([^\\)]+\\)","",as.character(subTax[,level]))
  counts$level <- subTax[match(counts$blastname, subTax$seqID),
                         "replaced"]
  
  #Get rid of the sequence names when aggregating, then sum up every sample's thing
  targs <- 
    colnames(counts)[colnames(counts) %notin% c("blastname", "seqID", "level")]
  levelCounts <- aggregate(counts[targs], counts["level"], FUN=sum)
  
  levelCounts$Total <- rowSums(levelCounts[,-grep("level", colnames(levelCounts))])
  levelCounts <- levelCounts[order(-levelCounts$Total),]
  
  rownames(levelCounts) <- levelCounts[,"level"]
  levelCounts <- levelCounts[,-1]
  
  levelCounts$Total <- NULL
  mainVal <- colSums(levelCounts)
  miniMat <- sweep(levelCounts, 2, mainVal, FUN = '/')
  
  
  return(miniMat)
}

colMax <- function(data) sapply(data, max, na.rm = TRUE)
colSort <- function(data, ...) sapply(data, sort, ...)

kingTableLarge <- relabundTax(countTableLarge, taxTable, "kingdom")
phylTableLarge <- relabundTax(countTableLarge, taxTable, "phylum")
classTableLarge <- relabundTax(countTableLarge, taxTable, "class")
orderTableLarge <- relabundTax(countTableLarge, taxTable, "order")
linTableLarge <- relabundTax(countTableLarge, taxTable, "lineage")
cladTableLarge <- relabundTax(countTableLarge, taxTable, "clade")
tribTableLarge <- relabundTax(countTableLarge, taxTable, "tribe")

kingTableSmall <- relabundTax(countTableSmall, taxTable, "kingdom")
phylTableSmall <- relabundTax(countTableSmall, taxTable, "phylum")
classTableSmall <- relabundTax(countTableSmall, taxTable, "class")
orderTableSmall <- relabundTax(countTableSmall, taxTable, "order")
linTableSmall <- relabundTax(countTableSmall, taxTable, "lineage")
cladTableSmall <- relabundTax(countTableSmall, taxTable, "clade")
tribTableSmall <- relabundTax(countTableSmall, taxTable, "tribe")

#Print out sorted vector of abundances of a certain taxa
#-sort(-tribTable["LD12",])

##Contour plot of groups of interest

library(akima)
library(viridis)
taxMesh <- function(group, table){
  subTable <- table[group,]
  
  subTable <- subTable[colnames(subTable) != "Total"]
  
  
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

levelTableLarge <- classTableLarge
levelTableSmall <- classTableSmall
numberOfGroups <- 2
groupsOfInterest <- c("Oxyphotobacteria", "Phycisphaerae")

makePlotOI <- function(groupsOI, levelTable, size, i){
  letter <- letters[i]
  actinos <- taxMesh(groupsOI, levelTable)
  actDates <- as.Date(actinos[[2]])
  actVals <- as.numeric(as.character(actinos[[1]]))
  actDepths <- as.numeric(as.character(actinos[[3]]))
  mldDates <- as.Date(actinos[[4]])
  mldDepths <- as.numeric(as.character(actinos[[5]]))
  
  fld2 <- interp(actDates, actDepths, actVals*100, duplicate = "mean", yo = seq(1, 90, 1), 
                 xo = seq(min(actDates), max(actDates), by="days"), extrap = TRUE)
  
  file <- paste0("./FigBinExtras/ClassDS/", groupsOI, size, ".png")
  dir.create(dirname(file), showWarnings = FALSE)
  
  png(file, width = 1470, height = 1050)
  
  filled.contour(x = seq(min(actDates), max(actDates), by = "day"),
                 y = fld2$y, z = fld2$z,
                 color.palette = function(x)viridis(x),
                 cex.main = 2.5, nlevels = 15,
                 cex = 2.5, cex.axis = 2.5, cex.lab = 2.5,
                 plot.axes={
                   if(i %in% c(3,4)){
                     axis(1, cex.axis=2.5, lwd.ticks = 1.5, tck = -0.015, padj = 0.5,
                          labels = format(seq(min(actDates), max(actDates),
                                              by = "month"),
                                          '%m-%y'),
                          at = as.numeric(seq(min(actDates), max(actDates),
                                              by = "month")));
                   } else{
                     axis(1, cex.axis=2.5, lwd.ticks = 1.5, tck = -0.015, padj = 0.5,
                          labels = NA,
                          at = as.numeric(seq(min(actDates), max(actDates),
                                              by = "month")));
                   }
                   lines(mldDates, mldDepths, col = "gray", lwd = 3);
                   if(i %in% c(1,3)){
                     axis(2, cex.axis=2.5, lwd.ticks = 1.5, at = seq(20, 100, by = 20))
                   } else{
                     axis(2, cex.axis=2.5, lwd.ticks = 1.5, at = seq(20, 100, by = 20),
                          labels = NA)
                   }
                   grid(col=rgb(1,1,1,0))
                   box(lwd = 1.5)
                 },
                 key.axes={
                   axis(4, cex.axis = 1.9) #####################
                 },
                 ylim = rev(range(fld2$y))
  )
  
  mtext("% Abundance", 3, line = 1, font = 2, adj = 1, cex = 1.5)
  mtext(letter, 3, line = 1, font = 2, adj = 0, cex = 2.5)
  
  dev.off()
}

letters <- c("A", "B", "C", "D")
makePlotOI(groupsOfInterest[1], levelTableSmall, "small", 1)
makePlotOI(groupsOfInterest[1], levelTableLarge, "large", 2)

makePlotOI(groupsOfInterest[2], levelTableSmall, "small", 3)
makePlotOI(groupsOfInterest[2], levelTableLarge, "large", 4)


library(png)

png("FigBinExtras/ClassDS/OxyPhcisPaneled.png", width = 2100, height = 1500)

par(mar = c(1, 3, 1, 0))
plot(0:2, 0:2, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
rasterImage(readPNG(source="FigBinExtras/ClassDS/Oxyphotobacteriasmall.png"), 0, 1, 1, 2)
rasterImage(readPNG(source="FigBinExtras/ClassDS/Oxyphotobacterialarge.png"), 1, 1, 2, 2)
rasterImage(readPNG(source="FigBinExtras/ClassDS/Phycisphaeraesmall.png"), 0, 0, 1, 1)
rasterImage(readPNG(source="FigBinExtras/ClassDS/Phycisphaeraelarge.png"), 1, 0, 2, 1)
mtext("Date", 1, line = -2, font = 2, cex = 2.5)
mtext("Oxyphotobacteria", 2, font = 2, cex = 2.5, adj = 0.8, line = -4)
mtext("Phycisphaerae", 2, font = 2, cex = 2.5, adj = 0.2, line = -4)
mtext("Depth (m)", 2, font = 2, cex = 2.5, line = -1)
mtext("Small (0.2-3 μm)", 3, font = 2, cex = 2.5, adj = 0.23, line = -4)
mtext("Large (>3 μm)", 3, font = 2, cex = 2.5, adj = 0.77, line = -4)
dev.off()