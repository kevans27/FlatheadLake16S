smallAbund0 <- read.csv("~/FlatheadMicrobes/SubsampledData/small_updatedNames.count_table", 
                        stringsAsFactors = FALSE, row.names = 1)

library(vegan)
library(cetcolor)

smallAbund <- smallAbund0

#Rename stuff
for (i in 1){
  colnames(smallAbund) <- gsub("Sept", "Sep", colnames(smallAbund))
  colnames(smallAbund) <- gsub("June", "Jun", colnames(smallAbund))
  colnames(smallAbund) <- gsub("July", "Jul", colnames(smallAbund))
  colnames(smallAbund) <- sub("[^_]*_(.*)", "\\1", colnames(smallAbund))
}


#Find dates for all
for (i in 1){
  uniqueSmall <- unique(gsub("_[^_]+$", "\\1", colnames(smallAbund)))[
    grep("_", unique(gsub("_[^_]+$", "\\1", colnames(smallAbund))))]
  
  datesSmall <- unique(as.Date(sapply(strsplit(uniqueSmall, "_"), `[[`, 1), format = "%d%b%y"))
  datesSmall <- datesSmall[order(datesSmall)]
}

smalls <- data.frame("Date" = datesSmall, "m10" = NA, "DCM" = NA, "m50" = NA, "m90" = NA)
for (i in 1:length(datesSmall)){
  samps <- grep(paste0("^", gsub("(?<![0-9])0+", "", format(datesSmall[i], "%d%b%y"), 
                                 perl = TRUE)), colnames(smallAbund))
  subFram <- t(smallAbund[, samps])
  subFram <- subFram[, colSums(subFram != 0) > 0]
  
  dfram <- vegdist(subFram, method = "jaccard", upper = TRUE)
  fMRow <- grep("05m", labels(dfram))
  
  if (length(fMRow) != 0){
    compRow <- as.matrix(dfram)[fMRow,]
    depths <- sapply(strsplit(names(compRow), "_"), `[[`, 2)
    
    if("10m" %in% depths){
      smalls[i, "m10"] <- compRow[depths == "10m"]
    }
    if("DCM" %in% depths){
      smalls[i, "DCM"] <- compRow[depths == "DCM"]
    }
    if("50m" %in% depths){
      smalls[i, "m50"] <- compRow[depths == "50m"]
    }
    if("90m" %in% depths){
      smalls[i, "m90"] <- compRow[depths == "90m"]
    }
  }
}

all5ms <- t(smallAbund[, grep("5m", colnames(smallAbund))])
rowlyDates <- unique(as.Date(sapply(strsplit(rownames(all5ms), "_"), `[[`, 1), format = "%d%b%y"))
all5ms <- all5ms[order(rowlyDates),]
rowlyDatesSorted <- rowlyDates[order(rowlyDates)]
daysSince <- as.numeric(rowlyDatesSorted - rowlyDatesSorted[1])
all5ms <- all5ms[, colSums(all5ms != 0) > 0]
jaccAll5ms <- as.matrix(vegdist(all5ms, method = "jaccard"))
dayDiff <- c()
jaccVal <- c()
for (i in 1:(nrow(jaccAll5ms)-1)){
  targRow <- unname(jaccAll5ms[i,(i+1):ncol(jaccAll5ms)])
  startVal <- daysSince[i]
  targDateDiff <- daysSince[(i+1):ncol(jaccAll5ms)] - startVal
  
  dayDiff <- c(dayDiff, targDateDiff)
  jaccVal <- c(jaccVal, targRow)
}

orderTarg <- order(dayDiff)
sortedDays <- dayDiff[orderTarg]
sortedJaccs <- jaccVal[orderTarg]

modelFit <- lm(sortedJaccs~sin(2*pi*sortedDays/365)+cos(2*pi*sortedDays/365))


smallAbundTrim <- smallAbund[, grep("Sequence", colnames(smallAbund), invert = TRUE)]

mCTS <- t(smallAbundTrim)
shannonSmall <- diversity(mCTS, "shannon")

smallShan <- data.frame("Date" = datesSmall, "m5_shann" = NA, "m10_shann" = NA, 
                        "DCM_shann" = NA, "m50_shann" = NA, "m90_shann" = NA)
for (i in 1:length(datesSmall)){
samps <- grep(paste0("^", gsub("(?<![0-9])0+", "", format(datesSmall[i], "%d%b%y"), 
                               perl = TRUE)), names(shannonSmall))
subFram <- shannonSmall[samps]

depths <- sapply(strsplit(names(subFram), "_"), `[[`, 2)

if("05m" %in% depths){
  smallShan[i, "m5_shann"] <- subFram[depths == "05m"]
}
if("10m" %in% depths){
  smallShan[i, "m10_shann"] <- subFram[depths == "10m"]
}
if("DCM" %in% depths){
  smallShan[i, "DCM_shann"] <- subFram[depths == "DCM"]
}
if("50m" %in% depths){
  smallShan[i, "m50_shann"] <- subFram[depths == "50m"]
}
if("90m" %in% depths){
  smallShan[i, "m90_shann"] <- subFram[depths == "90m"]
}
}

shannEpi <- c(smallShan$m5_shann, smallShan$m10_shann, smallShan$DCM_shann)
shannHyp <- c(smallShan$m50_shann, smallShan$m90_shann)
library(stats)
t.test(x = shannEpi, y = shannHyp, alternative = "less")


#smalls <- smalls[smalls$Date != as.Date("2017-05-15"),]

mldData <- read.csv("~/FlatheadPublic/MLDData.csv", stringsAsFactors = FALSE, row.names = 1)
mldData$Date <- as.Date(mldData$Date)
mldData <- mldData[order(mldData$Date),]

fullDF_0 <- merge(smalls, smallShan, by = "Date")
fullDF <- merge(fullDF_0, mldData, by = "Date")

mixStart <- as.Date(c("2016-12-15", "2017-11-29"))
mixEnd <- as.Date(c("2017-04-30", "2018-05-08"))

viridis4 <- cet_pal(7, name = "l16")[1:5]
pchs <- c(25, 18, 17, 19, 15)
ltys <- c(5, 1, 2, 3, 6)
cexSmall <- 0.9

tiff("FigBinExtras/jacMixingSmall_Double.tiff", 
     width = 7, height = 6, pointsize = 12, units = "in", res = 300)
plot.new()

par(new = "TRUE",plt = c(0.08,0.98,0.6,0.98),las = 1)
plot(fullDF$Date, fullDF$DCM_shann, xlab = "", xaxt = "n", ylim = c(4.5, 7), type = "b", 
     yaxt = "n", ylab = "", col = viridis4[3], lwd = 1, pch = pchs[2], lty = ltys[2])
rect(mixStart[1], -120, mixEnd[1], 20, col = "#B4B4B4", border = NA)
rect(mixStart[2], -120, mixEnd[2], 20, col = "#B4B4B4", border = NA)
abline(v=as.Date("2017-01-01"), lty = 2, col = "gray25", lwd = 1)
abline(v=as.Date("2018-01-01"), lty = 2, col = "gray25", lwd = 1)
lines(fullDF$Date, fullDF$m90_shann, col = viridis4[1], lwd = 1, pch = pchs[5], lty = ltys[5])
lines(fullDF$Date, fullDF$m50_shann, col = viridis4[2], lwd = 1, pch = pchs[4], lty = ltys[4])
lines(fullDF$Date, fullDF$DCM_shann, col = viridis4[3], lwd = 1, pch = pchs[3], lty = ltys[3])
lines(fullDF$Date, fullDF$m10_shann, col = viridis4[4], lwd = 1, pch = pchs[2], lty = ltys[2])
lines(fullDF$Date, fullDF$m5_shann, col = viridis4[5], lwd = 1, 
      pch = pchs[1], lty = ltys[1], bg = viridis4[5])
points(fullDF$Date, fullDF$m90_shann, col = viridis4[1], lwd = 1, pch = pchs[5], lty = ltys[5])
points(fullDF$Date, fullDF$m50_shann, col = viridis4[2], lwd = 1, pch = pchs[4], lty = ltys[4])
points(fullDF$Date, fullDF$DCM_shann, col = viridis4[3], lwd = 1, pch = pchs[3], lty = ltys[3])
points(fullDF$Date, fullDF$m10_shann, col = viridis4[4], lwd = 1, pch = pchs[2], lty = ltys[2])
points(fullDF$Date, fullDF$m5_shann, col = viridis4[5], lwd = 1, 
      pch = pchs[1], lty = ltys[1], bg = viridis4[5])
axis(1, lwd.ticks = 1, tck = -0.05, padj = 1, labels = NA, 
     at = as.numeric(seq(as.Date("2016-09-01"), max(smalls$Date), by = "month"), '%m-%y')
     [seq(1, length(format(seq(as.Date("2016-09-01"), max(smalls$Date), 
                               by = "month"), '%m-%y')), 3)])
axis(1, lwd.ticks = 1, tck = -0.02, padj = 1, 
     at = as.numeric(seq(as.Date("2016-09-01"), max(smalls$Date), by = "month"), 
                     '%m-%y'), labels = NA)
axis(2, tck = -0.03, lwd.ticks = 1, hadj = 1.2, labels = NA)
axLab <- seq(4.5, 7, by = 0.5)
axAdj <- seq(10.2, -9.2, length.out = length(axLab))
mtext(axLab, side = 2, line = 0.5, padj = axAdj, cex = cexSmall, xpd = NA)
box(lwd = 1)
mtext("Diversity (H´)", side = 2, line = 1.8, las = 0)
#text(as.Date("09-15-2016", format = "%m-%d-%Y"), 6.85, "a")


par(new = "TRUE",plt = c(0.08,0.98,0.18,0.56),las = 1)
plot(smalls$Date, smalls$DCM, xlab = "", xaxt = "n", ylim = c(0.5, 1), yaxt = "n", ylab = "", 
     col = viridis4[3], lwd = 1, pch = pchs[2], lty = ltys[2])
rect(mixStart[1], -120, mixEnd[1], 20, col = "#B4B4B4", border = NA)
rect(mixStart[2], -120, mixEnd[2], 20, col = "#B4B4B4", border = NA)
abline(v=as.Date("2017-01-01"), lty = 2, col = "gray25", lwd = 1)
abline(v=as.Date("2018-01-01"), lty = 2, col = "gray25", lwd = 1)
lines(smalls$Date, smalls$m90, col = viridis4[1], lwd = 1, pch = pchs[5], lty = ltys[5])
lines(smalls$Date, smalls$m50, col = viridis4[2], lwd = 1, pch = pchs[4], lty = ltys[4])
lines(smalls$Date, smalls$DCM, col = viridis4[3], lwd = 1, pch = pchs[3], lty = ltys[3])
lines(smalls$Date, smalls$m10, col = viridis4[4], lwd = 1, pch = pchs[2], lty = ltys[2])
points(smalls$Date, smalls$m90, col = viridis4[1], lwd = 1, pch = pchs[5], lty = ltys[5])
points(smalls$Date, smalls$m50, col = viridis4[2], lwd = 1, pch = pchs[4], lty = ltys[4])
points(smalls$Date, smalls$DCM, col = viridis4[3], lwd = 1, pch = pchs[3], lty = ltys[3])
points(smalls$Date, smalls$m10, col = viridis4[4], lwd = 1, pch = pchs[2], lty = ltys[2])
axis(1, lwd.ticks = 1, tck = -0.05, padj = -0.5,
     labels = format(seq(as.Date("2016-09-01"), max(smalls$Date), 
                         by = "month"), '%b')
     [seq(1, length(format(seq(as.Date("2016-09-01"), max(smalls$Date), 
                               by = "month"), '%m-%y')), 3)], 
     at = as.numeric(seq(as.Date("2016-09-01"), max(smalls$Date), 
                         by = "month"), '%m-%y')
     [seq(1, length(format(seq(as.Date("2016-09-01"), max(smalls$Date), 
                               by = "month"), '%m-%y')), 3)], cex = cexSmall)
axis(1, lwd.ticks = 1, tck = -0.02, padj = 1, 
     at = as.numeric(seq(as.Date("2016-09-01"), max(smalls$Date), by = "month"), 
                     '%m-%y'), labels = NA)
axis(2, tck = -0.03, lwd.ticks = 1, hadj = 0.5, cex = cexSmall)
box(lwd = 1)
mtext("Jaccard distance", side = 2, line = 1.8, las = 0)
#text(as.Date("09-15-2016", format = "%m-%d-%Y"), 0.98, "b")

mtext("2016", side = 1, line = 2, adj = 0.06, cex = cexSmall)
mtext("2017", side = 1, line = 2, adj = 0.375, cex = cexSmall)
mtext("2018", side = 1, line = 2, adj = 0.83, cex = cexSmall)
mtext("Date", side = 1, line = 3, adj = 0.5)
box(lwd = 1)

legend("bottom", inset=c(0,-0.5), legend = c("5 m", "10 m", "Chl max  ", "50 m", "90 m"), 
       col = rev(viridis4)[1:5], horiz = TRUE, pch = pchs[1:5], lty = ltys[1:5], xpd = NA, 
       bty = "n", text.width = c(40, 60, 100, 60), lwd = 1, box.lwd = 0, adj = 0.2, pt.bg = viridis4[5])

dev.off()



