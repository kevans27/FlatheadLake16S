smallAbund0 <- read.csv(
  "~/FlatheadMicrobes/SubsampledData/small_updatedNames.count_table", 
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
  
  datesSmall <- unique(as.Date(sapply(strsplit(uniqueSmall, "_"), `[[`, 1), 
                               format = "%d%b%y"))
  datesSmall <- datesSmall[order(datesSmall)]
}

smalls <- data.frame("Date" = datesSmall, "m10" = NA, "DCM" = NA, "m50" = NA, 
                     "m90" = NA)
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
rowlyDates <- unique(as.Date(sapply(strsplit(rownames(all5ms), "_"), `[[`, 1), 
                             format = "%d%b%y"))
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

plot(sortedDays, sortedJaccs, ylim = c(0.5, 1), xlab = "Days between samples",
     ylab = "Jaccard dissimilarity", main = "5 m samples")
abline(v = 365, lty = 2, col = "red")
abline(v = 730, lty = 2, col = "red")
lines(sortedDays,modelFit$fitted.values,col="blue", lty = 2)

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


#smalls <- smalls[smalls$Date != as.Date("2017-05-15"),]

mldData <- 
  read.csv("~/FlatheadPublic/MLDData.csv", stringsAsFactors = FALSE, row.names = 1)
mldData$Date <- as.Date(mldData$Date)
mldData <- mldData[order(mldData$Date),]

fullDF_0 <- merge(smalls, smallShan, by = "Date")
fullDF <- merge(fullDF_0, mldData, by = "Date")

mixStart <- as.Date(c("2016-12-15", "2017-11-29"))
mixEnd <- as.Date(c("2017-04-30", "2018-05-08"))

viridis4 <- cet_pal(7, name = "l16")[1:5]
pchs <- c(25, 18, 17, 19, 15)
ltys <- c(5, 1, 2, 3, 6)

png("FigBinExtras/jacMixingSmall.png", width = 1200, height = 1200)
plot.new()

par(new = "TRUE",plt = c(0.14,0.99,0.8,0.98),las = 1, cex.axis = 2.5)
plot(fullDF$Date, -fullDF$MLD, xlab = "", xaxt = "n", ylim = c(-100, 0), type = "b",
     yaxt = "n", ylab = "", cex = 2.5, lwd = 2, pch = pchs[2], 
     lty = ltys[2])
rect(mixStart[1], -120, mixEnd[1], 20, col = "#B4B4B4")
rect(mixStart[2], -120, mixEnd[2], 20, col = "#B4B4B4")
abline(v=as.Date("2017-01-01"), lty = 2, col = "gray75", lwd = 3)
abline(v=as.Date("2018-01-01"), lty = 2, col = "gray75", lwd = 3)
lines(fullDF$Date, -fullDF$MLD, type = "b", cex = 2.5, lwd = 2, 
      pch = pchs[2], lty = ltys[2])
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, padj = 1,
     labels = NA, 
     at = as.numeric(seq(as.Date("2016-09-01"), max(smalls$Date), 
                         by = "month"), '%m-%y')
     [seq(1, length(format(seq(as.Date("2016-09-01"), max(smalls$Date), by = "month"), 
                           '%m-%y')), 3)], font = 2)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, padj = 1, 
     at = as.numeric(seq(as.Date("2016-09-01"), max(smalls$Date), by = "month"), '%m-%y'),
     labels = NA)
axis(2, cex.axis=2.5, tck = -0.03, lwd.ticks = 3, hadj = 1.2, font = 2, 
     at = seq(-80, 0, by = 20), labels = seq(80, 0, by = -20))
box(lwd = 3)
mtext("Mixed Layer", side = 2, line = 9, las = 0, 
      font = 2, cex = 3)
mtext("(m)", side = 2, line = 6, las = 0, 
      font = 2, cex = 3)
text(as.Date("09-15-2016", format = "%m-%d-%Y"), -6, "A", font = 2, cex = 3)

box(lwd = 3)

par(new = "TRUE",plt = c(0.14,0.99,0.48,0.77),las = 1, cex.axis = 2.5)
plot(fullDF$Date, fullDF$DCM_shann, xlab = "", xaxt = "n", ylim = c(4.5, 7), type = "b",
     yaxt = "n", ylab = "", col = viridis4[3], cex = 2.5, lwd = 2, pch = pchs[2], 
     lty = ltys[2])
rect(mixStart[1], -120, mixEnd[1], 20, col = "#B4B4B4")
rect(mixStart[2], -120, mixEnd[2], 20, col = "#B4B4B4")
abline(v=as.Date("2017-01-01"), lty = 2, col = "gray75", lwd = 3)
abline(v=as.Date("2018-01-01"), lty = 2, col = "gray75", lwd = 3)
lines(fullDF$Date, fullDF$m90_shann, type = "b", col = viridis4[1], cex = 2.5, lwd = 2, 
      pch = pchs[5], lty = ltys[5])
lines(fullDF$Date, fullDF$m50_shann, type = "b", col = viridis4[2], cex = 2.5, lwd = 2, 
      pch = pchs[4], lty = ltys[4])
lines(fullDF$Date, fullDF$DCM_shann, type = "b", col = viridis4[3], cex = 2.5, lwd = 2, 
      pch = pchs[3], lty = ltys[3])
lines(fullDF$Date, fullDF$m10_shann, type = "b", col = viridis4[4], cex = 2.5, lwd = 2, 
      pch = pchs[2], lty = ltys[2])
lines(fullDF$Date, fullDF$m5_shann, type = "b", col = viridis4[5], cex = 1.5, lwd = 3, 
      pch = pchs[1], lty = ltys[1], bg = viridis4[5])
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, padj = 1,
     labels = NA, 
     at = as.numeric(seq(as.Date("2016-09-01"), max(smalls$Date), 
                         by = "month"), '%m-%y')
     [seq(1, length(format(seq(as.Date("2016-09-01"), max(smalls$Date), by = "month"), 
                           '%m-%y')), 3)], font = 2)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, padj = 1, 
     at = as.numeric(seq(as.Date("2016-09-01"), max(smalls$Date), by = "month"), '%m-%y'),
     labels = NA)
axis(2, cex.axis=2.5, tck = -0.03, lwd.ticks = 3, hadj = 1.2, font = 2)
box(lwd = 3)
mtext("Diversity", side = 2, line = 9, las = 0, 
      font = 2, cex = 3)
mtext("(H)", side = 2, line = 6, las = 0, 
      font = 2, cex = 3)
text(as.Date("09-15-2016", format = "%m-%d-%Y"), 6.85, "B", font = 2, cex = 3)


par(new = "TRUE",plt = c(0.14,0.99,0.16,0.45),las = 1, cex.axis = 2.5)
plot(smalls$Date, smalls$DCM, xlab = "", xaxt = "n", ylim = c(0.5, 1), type = "b",
     yaxt = "n", ylab = "", col = viridis4[3], cex = 2.5, lwd = 2, pch = pchs[2], 
     lty = ltys[2])
rect(mixStart[1], -120, mixEnd[1], 20, col = "#B4B4B4")
rect(mixStart[2], -120, mixEnd[2], 20, col = "#B4B4B4")
abline(v=as.Date("2017-01-01"), lty = 2, col = "gray75", lwd = 3)
abline(v=as.Date("2018-01-01"), lty = 2, col = "gray75", lwd = 3)
lines(smalls$Date, smalls$m90, type = "b", col = viridis4[1], cex = 2.5, lwd = 2, 
      pch = pchs[5], lty = ltys[5])
lines(smalls$Date, smalls$m50, type = "b", col = viridis4[2], cex = 2.5, lwd = 2, 
      pch = pchs[4], lty = ltys[4])
lines(smalls$Date, smalls$DCM, type = "b", col = viridis4[3], cex = 2.5, lwd = 2, 
      pch = pchs[3], lty = ltys[3])
lines(smalls$Date, smalls$m10, type = "b", col = viridis4[4], cex = 2.5, lwd = 2, 
      pch = pchs[2], lty = ltys[2])
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, padj = 1,
     labels = format(seq(as.Date("2016-09-01"), max(smalls$Date), 
                         by = "month"), '%b')
     [seq(1, length(format(seq(as.Date("2016-09-01"), max(smalls$Date), by = "month"), 
                           '%m-%y')), 3)], 
     at = as.numeric(seq(as.Date("2016-09-01"), max(smalls$Date), 
                         by = "month"), '%m-%y')
     [seq(1, length(format(seq(as.Date("2016-09-01"), max(smalls$Date), by = "month"), 
                           '%m-%y')), 3)], font = 2)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.02, padj = 1, 
     at = as.numeric(seq(as.Date("2016-09-01"), max(smalls$Date), by = "month"), '%m-%y'),
     labels = NA)
axis(2, cex.axis=2.5, tck = -0.03, lwd.ticks = 3, hadj = 1.2, font = 2)
box(lwd = 3)
mtext("Jaccard Distance", side = 2, line = 7.5, las = 0, 
      font = 2, cex = 3)
text(as.Date("09-15-2016", format = "%m-%d-%Y"), 0.98, "C", font = 2, cex = 3)

mtext("2016", side = 1, line = 5.1, cex = 2.5, font = 2, adj = 0.06)
mtext("2017", side = 1, line = 5.1, cex = 2.5, font = 2, adj = 0.375)
mtext("2018", side = 1, line = 5.1, cex = 2.5, font = 2, adj = 0.83)
box(lwd = 3)

legend("bottom", inset=c(0,-0.5), 
       legend = c("5 m", "10 m", "Chl max  ", "50 m", "90 m"), 
       col = rev(viridis4)[1:5], horiz = TRUE, pch = pchs[1:5], lty = ltys[1:5], 
       xpd = NA, pt.cex = 2.5, pt.bg = viridis4[5], text.width = c(40, 60, 100, 60),
       cex = 2.5, lwd = 2.5, box.lwd = 0, text.font = 2, adj = 0.2)

dev.off()






richnessSmall <- specnumber(mCTS)
shannonSmall <- diversity(mCTS, "shannon")

smallDiv <- data.frame("Date" = datesSmall, "m5_rich" = NA, "m10_rich" = NA, 
                       "DCM_rich" = NA, "m50_rich" = NA, "m90_rich" = NA, 
                       "m5_even" = NA, "m10_even" = NA, "DCM_even" = NA, 
                       "m50_even" = NA, "m90_even" = NA, "m5_shann" = NA, 
                       "m10_shann" = NA, "DCM_shann" = NA, "m50_shann" = NA,
                       "m90_shann" = NA)
for (i in 1:length(datesSmall)){
  samps <- grep(paste0("^", gsub("(?<![0-9])0+", "", format(datesSmall[i], "%d%b%y"), 
                                 perl = TRUE)), names(shannonSmall))
  subFram <- shannonSmall[samps]
  
  depths <- sapply(strsplit(names(subFram), "_"), `[[`, 2)
  
  if("05m" %in% depths){
    smallDiv[i, "m5_shann"] <- subFram[depths == "05m"]
  }
  if("10m" %in% depths){
    smallDiv[i, "m10_shann"] <- subFram[depths == "10m"]
  }
  if("DCM" %in% depths){
    smallDiv[i, "DCM_shann"] <- subFram[depths == "DCM"]
  }
  if("50m" %in% depths){
    smallDiv[i, "m50_shann"] <- subFram[depths == "50m"]
  }
  if("90m" %in% depths){
    smallDiv[i, "m90_shann"] <- subFram[depths == "90m"]
  }
}
for (i in 1:length(datesSmall)){
  samps <- grep(paste0("^", gsub("(?<![0-9])0+", "", format(datesSmall[i], "%d%b%y"),
                                 perl = TRUE)), names(richnessSmall))
  subFram <- richnessSmall[samps]
  
  depths <- sapply(strsplit(names(subFram), "_"), `[[`, 2)
  
  if("05m" %in% depths){
    smallDiv[i, "m5_rich"] <- subFram[depths == "05m"]
  }
  if("10m" %in% depths){
    smallDiv[i, "m10_rich"] <- subFram[depths == "10m"]
  }
  if("DCM" %in% depths){
    smallDiv[i, "DCM_rich"] <- subFram[depths == "DCM"]
  }
  if("50m" %in% depths){
    smallDiv[i, "m50_rich"] <- subFram[depths == "50m"]
  }
  if("90m" %in% depths){
    smallDiv[i, "m90_rich"] <- subFram[depths == "90m"]
  }
}
for (i in 1){
  smallDiv$m5_even <- smallDiv$m5_shann/log(smallDiv$m5_rich)
  smallDiv$m10_even <- smallDiv$m10_shann/log(smallDiv$m10_rich)
  smallDiv$DCM_even <- smallDiv$DCM_shann/log(smallDiv$DCM_rich)
  smallDiv$m50_even <- smallDiv$m50_shann/log(smallDiv$m50_rich)
  smallDiv$m90_even <- smallDiv$m90_shann/log(smallDiv$m90_rich)
}

rownames(smallDiv) <- smallDiv$Date
smallDiv$Date <- NULL
colMeans(smallDiv, na.rm = TRUE)
