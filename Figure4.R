taxTable <- read.csv("RawData/prok.98.taxonomy")
countTableSmall <- read.csv("SubsampledData/small_updatedNames.count_table",
                            stringsAsFactors = FALSE, 
                            row.names = 1)
transKey <- read.csv("RawData/prok.name.key")
`%notin%` <- Negate(`%in%`)
library(cetcolor)

otuNames <- countTableSmall$Representative_Sequence

#Rename/trim things
for (i in 1){
  countTableSmall$Representative_Sequence <- NULL
  
  countTableSmall$blastname <- rownames(countTableSmall)
  colnames(countTableSmall) <- sub("[^_]*_(.*)", "\\1", colnames(countTableSmall))
  colnames(countTableSmall) <- gsub("_[^_]+$", "\\1", colnames(countTableSmall))
  colnames(countTableSmall) <- gsub("Sept", "Sep", colnames(countTableSmall))
  colnames(countTableSmall) <- gsub("June", "Jun", colnames(countTableSmall))
}

`%notin%` <- Negate(`%in%`)
# Function to plot color bar
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), 
                      title='') {
  scale = (length(lut)-1)/(max-min)
  
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', 
       main=title)
  axis(2, ticks, las=1)
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

phylTableSmall <- relabundTax(countTableSmall, taxTable, "phylum")


cyano5 <- phylTableSmall["Cyanobacteria", grep("5m", colnames(phylTableSmall))]
dates5 <- as.Date(sapply(strsplit(colnames(cyano5), "_"), `[[`, 1), 
                    format = "%d%b%y")
cyano5 <- cyano5[order(dates5)]
dates5 <- dates5[order(dates5)]

cyano10 <- phylTableSmall["Cyanobacteria", grep("10m", colnames(phylTableSmall))]
dates10 <- as.Date(sapply(strsplit(colnames(cyano10), "_"), `[[`, 1), 
                  format = "%d%b%y")
cyano10 <- cyano10[order(dates10)]
dates10 <- dates10[order(dates10)]

cyano50 <- phylTableSmall["Cyanobacteria", grep("50m", colnames(phylTableSmall))]
dates50 <- as.Date(sapply(strsplit(colnames(cyano50), "_"), `[[`, 1), 
                  format = "%d%b%y")
cyano50 <- cyano50[order(dates50)]
dates50 <- dates50[order(dates50)]

cyano90 <- phylTableSmall["Cyanobacteria", grep("90m", colnames(phylTableSmall))]
dates90 <- as.Date(sapply(strsplit(colnames(cyano90), "_"), `[[`, 1), 
                   format = "%d%b%y")
cyano90 <- cyano90[order(dates90)]
dates90 <- dates90[order(dates90)]


targDf <- countTableSmall
targDf$blastname <- NULL


mainVal <- colSums(targDf[,grep("Sequence", colnames(targDf), invert = TRUE)])
targDf <- sweep(targDf[,grep("Sequence", colnames(targDf), invert = TRUE)], 
                2, mainVal, FUN = '/')

targDf$Sequence <- otuNames

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

allPossibles <- taxTable[grep("Cyanobacteria", taxTable[,"phylum"]), "seqID"]

newDf5m <- frame5m[rownames(frame5m) %in% allPossibles,]
rownames(newDf5m) <- newDf5m$Sequence
newDf5m$Sequence <- NULL
colnames(newDf5m) <- subDates5
newDf10m <- frame10m[rownames(frame10m) %in% allPossibles,]
rownames(newDf10m) <- newDf10m$Sequence
newDf10m$Sequence <- NULL
colnames(newDf10m) <- subDates10
newDfDCM <- frameDCM[rownames(frameDCM) %in% allPossibles,]
rownames(newDfDCM) <- newDfDCM$Sequence
newDfDCM$Sequence <- NULL
colnames(newDfDCM) <- subDatesDCM
newDf50m <- frame50m[rownames(frame50m) %in% allPossibles,]
rownames(newDf50m) <- newDf50m$Sequence
newDf50m$Sequence <- NULL
colnames(newDf50m) <- subDates50
newDf90m <- frame90m[rownames(frame90m) %in% allPossibles,]
rownames(newDf90m) <- newDf90m$Sequence
newDf90m$Sequence <- NULL
colnames(newDf90m) <- subDates90


library(viridis)
cols <- viridis(6)[c(1:5)]
depthLegend <- c("5 m", "10 m", "50 m", "90 m")

allMonths <- seq.Date(as.Date("2016-03-01"), as.Date("2020-03-01"), by = "month")
allThirds <- allMonths[c(TRUE, FALSE, FALSE)]
allTabs <- allMonths[c(FALSE, TRUE, TRUE)]



earlyFirst <- as.Date("2016-12-01")
earlyLast <- as.Date("2017-05-01")
lateFirst <- as.Date("2017-12-01")
lateLast <- as.Date("2018-05-15")


targetPhylum <- "Cyanobacteria"
png("~/FlatheadMicrobes/FigBinExtras/cyanobacteria_5_90_Rel.tiff", width = 7, 
    height = 6, pointsize = 12, units = "in", res = 1200)
plot.new()
par(new = "TRUE",plt = c(0.1,0.95,0.25,0.9),las = 1, xpd = FALSE)
plot(dates5, cyano5, col = cols[5], ylim = c(0, 0.4), xlab = "", ylab = "", 
     xaxt = "n", yaxt = "n", type = 'b')
rect(earlyFirst, -10, earlyLast, 200, col = "lightgray", border = NA)
rect(lateFirst, -10, lateLast, 200, col = "lightgray", border = NA)
lines(dates5, cyano5, col = cols[5])
lines(dates10, cyano10, col = cols[4])
lines(dates50, cyano50, col = cols[2])
lines(dates90, cyano90, col = cols[1])
points(dates5, cyano5, col = cols[5])
points(dates10, cyano10, col = cols[4])
points(dates50, cyano50, col = cols[2])
points(dates90, cyano90, col = cols[1])
axis(1, tck = -0.05, padj = 1,
     labels = format(allThirds, "%b"), at = allThirds)
axis(1, tck = -0.03, padj = 1, 
     at = allTabs, labels = NA)
axis(2, tck = -0.02, padj = 1, 
     at = seq(0.1, 0.5, by = 0.2), labels = NA)
axis(2, tck = -0.04, hadj = 1.2, las = 1,
     labels = seq(0, 40, by = 20), at = seq(0, 0.4, by = 0.2))
box(lwd = 1)
mtext(bquote(bold("Relative abundance (%)")), 
      side = 2, line = 3.5, las = 0)
text(as.Date("2016-10-01"), -0.092, "2016", xpd = NA)
text(as.Date("2017-07-01"), -0.092, "2017", xpd = NA)
text(as.Date("2018-07-01"), -0.092, "2018", xpd = NA)
abline(v = seq.Date(as.Date("2016-01-01"), as.Date("2021-01-01"), by = "year"),
       lty = 2, col = "darkgrey")

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

legend("bottom", legend = depthLegend, bty = "n",
       fill = cols[c(5,4,2,1)], horiz = TRUE, box.lwd = 0)

dev.off()


mldData <- 
  read.csv("~/FlatheadPublic/MLDData.csv", stringsAsFactors = FALSE, row.names = 1)
mldData$Date <- as.Date(mldData$Date)
mldData <- mldData[order(mldData$Date),]

spaceBottom <- 0.12
spaceTop <- 0.98
spaceL <- 0.1
spaceR <- 0.95
nOTUs <- 3
#mldSpace <- 0.13
vertGap <- 0.02
cyanoVert <- (spaceTop - spaceBottom - vertGap * (nOTUs + 1)) / (nOTUs + 1)
xLim <- as.Date(c("2016-09-01", "2019-01-01"))
ymaxes <- c(30, 5.2, 5.2, 5.2)

labels <- c("Aggregate Cyanos", "OTU3", "OTU28", "OTU36")
#labels <- c("b", "c", "d", "e")

ptCex <- 1
cexSmall <- 1

tiff("~/FlatheadMicrobes/FigBinExtras/cyanobacteria_5_90_Rel_OTUs.tiff", 
     width = 7, height = 7, pointsize = 12, units = "in", res = 1200)
plot.new()


###Mixed layer
#par(new = "TRUE",plt = c(spaceL,spaceR,spaceTop - mldSpace,spaceTop),las = 1)
#plot(mldData$Date, -mldData$MLD, xlab = "", xaxt = "n", ylim = c(-100, 0), 
#     type = "l", yaxt = "n", ylab = "", cex = ptCex, xlim = xLim)
#rect(earlyFirst, -120, earlyLast, 20, col = "lightgray", border = NA)
#rect(lateFirst, -120, lateLast, 20, col = "lightgray", border = NA)
#abline(v=as.Date("2017-01-01"), lty = 2, col = "gray25")
#abline(v=as.Date("2018-01-01"), lty = 2, col = "gray25")
#lines(mldData$Date, -mldData$MLD, type = "l", cex = ptCex)
#points(mldData$Date, -mldData$MLD)
#axis(1, tck = -0.05, padj = 1,
#     labels = NA, at = allThirds)
#axis(1, tck = -0.03, padj = 1, 
#     at = allTabs, labels = NA)
#axis(4, tck = -0.03, 
#     at = seq(-20, -100, by = -40), labels = NA)
#axis(4, tck = -0.05,  adj = 0, hadj = 0, las = 1,
#     labels = NA, at = seq(0, -80, by = -40))
#axLab <- seq(80, 0, by = -40)
#axAdj <- seq(3, -3.5, length.out = length(axLab))
#mtext(axLab, side = 4, line = 0.5, padj = axAdj, cex = cexSmall)
#box(lwd = 1)
#text(as.Date("2019-05-01"), -46, "Mixed Layer (m)", 
#     srt = 270, xpd = NA)



###Agg cyanos
par(new = "TRUE",plt = c(spaceL,spaceR,
                         spaceTop - cyanoVert,
                         spaceTop),las = 1, xpd = FALSE)
plot(dates5, cyano5, col = cols[5], ylim = c(0, 0.4), xlab = "", ylab = "", 
     xaxt = "n", yaxt = "n", type = 'b', xlim = xLim, cex = ptCex)
rect(earlyFirst, -10, earlyLast, 200, col = "lightgray", border = NA)
rect(lateFirst, -10, lateLast, 200, col = "lightgray", border = NA)
lines(dates5, cyano5, col = cols[5])
lines(dates90, cyano90, col = cols[1])
points(dates5, cyano5, col = cols[5], cex = ptCex)
points(dates90, cyano90, col = cols[1], cex = ptCex)
axis(1, tck = -0.05, padj = 1,
     labels = NA, at = allThirds)
axis(1, tck = -0.03, padj = 1, 
     at = allTabs, labels = NA)
axis(2, tck = -0.03, padj = 1, 
     at = seq(0.1, 0.5, by = 0.2), labels = NA)
axis(2, tck = -0.05, adj = 0, hadj = 0.4, las = 1,
     labels = NA, at = seq(0, 0.4, by = 0.2))
axLab <- seq(0, 40, by = 20)
axAdj <- seq(5.8, -4.8, length.out = length(axLab))
mtext(axLab, side = 2, line = 0.5, padj = axAdj, cex = cexSmall)
box(lwd = 1)
text(as.Date("2016-07-25"), 0.35, labels[1], pos = 4)
abline(v = seq.Date(as.Date("2016-01-01"), as.Date("2021-01-01"), by = "year"),
       lty = 2, col = "gray25")


for (i in 1:nOTUs){
  ###Otu3
  par(new = "TRUE",plt = c(spaceL,spaceR,
                           spaceTop - vertGap*(i) - cyanoVert*(i+1),
                           spaceTop - vertGap*(i) - cyanoVert*(i)),
      las = 1, xpd = FALSE)
  plot(as.Date(colnames(newDf5m)), newDf5m[i,]*100, yaxt = "n", xaxt = "n", 
       ylab = "", xlab = "", type = "b", col = cols[5], cex = ptCex, 
       ylim = c(0, ymaxes[i]), xlim = xLim)
  rect(earlyFirst, -10, earlyLast, 200, col = "lightgray", border = NA)
  rect(lateFirst, -10, lateLast, 200, col = "lightgray", border = NA)
  lines(as.Date(colnames(newDf5m)), newDf5m[i,]*100, col = cols[5],)
  lines(as.Date(colnames(newDf90m)), newDf90m[i,]*100, col = cols[1])
  points(as.Date(colnames(newDf5m)), newDf5m[i,]*100, col = cols[5], cex = ptCex)
  points(as.Date(colnames(newDf90m)), newDf90m[i,]*100, col = cols[1], 
         cex = ptCex)
  if (i == nOTUs){
  axis(1, tck = -0.05, padj = -1, labels = format(allThirds, "%b"), at = allThirds,
       cex = cexSmall)
    text(as.Date("2016-10-01"), -1.9, "2016", xpd = NA, cex = cexSmall)
    text(as.Date("2017-07-01"), -1.9, "2017", xpd = NA, cex = cexSmall)
    text(as.Date("2018-07-01"), -1.9, "2018", xpd = NA, cex = cexSmall)
    text(as.Date("2017-10-01"), -2.6, "Date", xpd = NA)
  }else{
    axis(1, tck = -0.05, padj = 1, labels = NA, at = allThirds)
  }
  axis(1, tck = -0.03, padj = 1, at = allTabs, labels = NA)
  if (i == 1){
    axis(4, tck = -0.03, padj = 1, at = seq(7.5, 37.5, by = 15), labels = NA)
    axis(4, tck = -0.05, adj = 0, hadj = 0.6, las = 1,
         labels = NA, at = seq(0, 30, by = 15))
    axLab <- seq(0, 30, by = 15)
    axAdj <- seq(5.8, -4.8, length.out = length(axLab))
    mtext(axLab, side = 4, line = 0.5, padj = axAdj, cex = cexSmall)
    text(as.Date("2016-07-25"), 26.25, labels[i+1], pos = 4)
  }
  if (i == 2){
    axis(2, tck = -0.03, padj = 1, at = seq(1.25, 10, by = 2.5), labels = NA)
    axis(2, tck = -0.05, adj = 0, hadj = 0.4, las = 1,
         labels = NA, at = seq(0, 5, by = 2.5))
    axLab <- seq(0, 5, by = 2.5)
    axAdj <- seq(5.8, -4.4, length.out = length(axLab))
    mtext(axLab, side = 2, line = 0.5, padj = axAdj, cex = cexSmall)
    text(as.Date("2016-07-25"), 4.375, labels[i+1], pos = 4)
  }
  if (i == 3){
    axis(4, tck = -0.03, padj = 1, at = seq(1.25, 10, by = 2.5), labels = NA)
    axis(4, tck = -0.05, adj = 0, hadj = 0.5, las = 1,
         labels = NA, at = seq(0, 5, by = 2.5))
    axLab <- seq(0, 5, by = 2.5)
    axAdj <- seq(5.8, -4.5, length.out = length(axLab))
    mtext(axLab, side = 4, line = 0.5, padj = axAdj, cex = cexSmall)
    text(as.Date("2016-07-25"), 4.375, labels[i+1], pos = 4)
  }
  box(lwd = 1)
  if (i == 2){
    mtext(bquote("Relative abundance (%)"), side = 2, line = 2, las = 0, adj = -1.5)
  }
  abline(v = seq.Date(as.Date("2016-01-01"), as.Date("2021-01-01"), by = "year"),
         lty = 2, col = "gray25")
}

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

legend(-0.3, -0.98, legend = depthLegend[c(1,4)], bty = "n",
       fill = cols[c(5,1)], horiz = TRUE, box.lwd = 0, xpd = NA)

dev.off()






png("~/FlatheadMicrobes/FigBinExtras/cyanobacteria_Four_Rel_OTUs.png", width = 800, 
    height = 1200)
plot.new()


###Mixed layer
par(new = "TRUE",plt = c(spaceL,spaceR,spaceTop - mldSpace,spaceTop),las = 1)
plot(mldData$Date, -mldData$MLD, xlab = "", xaxt = "n", ylim = c(-100, 0), 
     type = "b", yaxt = "n", ylab = "", cex = ptCex, xlim = xLim)
rect(earlyFirst, -120, earlyLast, 20, col = "lightgray", border = NA)
rect(lateFirst, -120, lateLast, 20, col = "lightgray", border = NA)
abline(v=as.Date("2017-01-01"), lty = 2, col = "gray75")
abline(v=as.Date("2018-01-01"), lty = 2, col = "gray75")
lines(mldData$Date, -mldData$MLD, type = "b", cex = ptCex)
axis(1, tck = -0.05, padj = 1,
     labels = NA, at = allThirds)
axis(1, tck = -0.03, padj = 1, 
     at = allTabs, labels = NA)
axis(2, tck = -0.02, padj = 1, 
     at = seq(-20, -100, by = -40), labels = NA)
axis(2, tck = -0.04, hadj = 1.2, las = 1,
     labels = seq(0, 80, by = 40), at = seq(0, -80, by = -40))
box(lwd = 1)
mtext("Mixed Layer (m)", side = 2, line = 3.5, las = 0)


###Agg cyanos
par(new = "TRUE",plt = c(spaceL,spaceR,
                         spaceTop - mldSpace - vertGap - cyanoVert,
                         spaceTop - mldSpace - vertGap),las = 1, xpd = FALSE)
plot(dates5, cyano5, col = cols[5], ylim = c(0, 0.4), xlab = "", ylab = "", 
     xaxt = "n", yaxt = "n", type = 'b', xlim = xLim, cex = ptCex)
rect(earlyFirst, -10, earlyLast, 200, col = "lightgray", border = NA)
rect(lateFirst, -10, lateLast, 200, col = "lightgray", border = NA)
lines(dates5, cyano5, col = cols[5])
lines(dates10, cyano10, col = cols[4])
lines(dates50, cyano50, col = cols[2])
lines(dates90, cyano90, col = cols[1])
points(dates5, cyano5, col = cols[5], cex = ptCex)
points(dates10, cyano10, col = cols[4], cex = ptCex)
points(dates50, cyano50, col = cols[2], cex = ptCex)
points(dates90, cyano90, col = cols[1], cex = ptCex)
axis(1, tck = -0.05, padj = 1,
     labels = NA, at = allThirds)
axis(1, tck = -0.03, padj = 1, 
     at = allTabs, labels = NA)
axis(2, tck = -0.02, padj = 1, 
     at = seq(0.1, 0.5, by = 0.2), labels = NA)
axis(2, tck = -0.04, hadj = 1.2, las = 1,
     labels = seq(0, 40, by = 20), at = seq(0, 0.4, by = 0.2))
box(lwd = 1)
text(as.Date("2016-07-25"), 0.35, labels[1], pos = 4)
abline(v = seq.Date(as.Date("2016-01-01"), as.Date("2021-01-01"), by = "year"),
       lty = 2, col = "darkgrey")


for (i in 1:nOTUs){
  ###Otu3
  par(new = "TRUE",plt = c(spaceL,spaceR,
                           spaceTop - mldSpace - vertGap*(i+1) - cyanoVert*(i+1),
                           spaceTop - mldSpace - vertGap*(i+1) - cyanoVert*i),
      las = 1, xpd = FALSE)
  plot(as.Date(colnames(newDf5m)), newDf5m[i,]*100, yaxt = "n", xaxt = "n", 
       ylab = "", xlab = "", type = "b", col = cols[5], cex = ptCex, 
       ylim = c(0, ymaxes[i]), xlim = xLim)
  rect(earlyFirst, -10, earlyLast, 200, col = "lightgray", border = NA)
  rect(lateFirst, -10, lateLast, 200, col = "lightgray", border = NA)
  lines(as.Date(colnames(newDf5m)), newDf5m[i,]*100, col = cols[5],)
  lines(as.Date(colnames(newDf10m)), newDf10m[i,]*100, col = cols[4])
  lines(as.Date(colnames(newDf50m)), newDf50m[i,]*100, col = cols[2])
  lines(as.Date(colnames(newDf90m)), newDf90m[i,]*100, col = cols[1])
  points(as.Date(colnames(newDf5m)), newDf5m[i,]*100, col = cols[5], 
         cex = ptCex)
  points(as.Date(colnames(newDf10m)), newDf10m[i,]*100, col = cols[4], 
         cex = ptCex)
  points(as.Date(colnames(newDf50m)), newDf50m[i,]*100, col = cols[2], 
         cex = ptCex)
  points(as.Date(colnames(newDf90m)), newDf90m[i,]*100, col = cols[1], 
         cex = ptCex)
  if (i == nOTUs){
    axis(1, tck = -0.05, padj = 1,
         labels = format(allThirds, "%b"), at = allThirds)
    text(as.Date("2016-10-01"), -2.3, "2016", xpd = NA)
    text(as.Date("2017-07-01"), -2.3, "2017", xpd = NA)
    text(as.Date("2018-07-01"), -2.3, "2018", xpd = NA)
  }else{
    axis(1, tck = -0.05, padj = 1,
         labels = NA, at = allThirds)
  }
  axis(1, tck = -0.03, padj = 1, 
       at = allTabs, labels = NA)
  if (i == 1){
    axis(2, tck = -0.02, padj = 1, 
         at = seq(7.5, 37.5, by = 15), labels = NA)
    axis(2, tck = -0.04, hadj = 1.2, las = 1,
         labels = seq(0, 30, by = 15), at = seq(0, 30, by = 15))
    text(as.Date("2016-07-25"), 26.25, labels[i+1], pos = 4)
  }else{
    axis(2, tck = -0.02, padj = 1, 
         at = seq(1.25, 10, by = 2.5), labels = NA)
    axis(2, tck = -0.04, hadj = 1.2, las = 1,
         labels = seq(0, 5, by = 2.5), at = seq(0, 5, by = 2.5))
    text(as.Date("2016-07-25"), 4.375, labels[i+1], pos = 4)
  }
  box(lwd = 1)
  if (i == 2){
    mtext(bquote(bold("Relative abundance (%)")), 
          side = 2, line = 3.5, las = 0)
  }
  abline(v = seq.Date(as.Date("2016-01-01"), as.Date("2021-01-01"), by = "year"),
         lty = 2, col = "darkgrey")
}

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

legend("bottom", legend = depthLegend, bty = "n",
       fill = cols[c(5,4,2,1)], horiz = TRUE, box.lwd = 0)

dev.off()