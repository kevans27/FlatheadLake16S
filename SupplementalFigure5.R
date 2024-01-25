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

smalls5 <- data.frame("Date" = datesSmall, "m10" = NA, "DCM" = NA, "m50" = NA, 
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
      smalls5[i, "m10"] <- compRow[depths == "10m"]
    }
    if("DCM" %in% depths){
      smalls5[i, "DCM"] <- compRow[depths == "DCM"]
    }
    if("50m" %in% depths){
      smalls5[i, "m50"] <- compRow[depths == "50m"]
    }
    if("90m" %in% depths){
      smalls5[i, "m90"] <- compRow[depths == "90m"]
    }
  }
}

smalls10 <- data.frame("Date" = datesSmall, "DCM" = NA, "m50" = NA, 
                      "m90" = NA)
for (i in 1:length(datesSmall)){
  samps <- grep(paste0("^", gsub("(?<![0-9])0+", "", format(datesSmall[i], "%d%b%y"), 
                                 perl = TRUE)), colnames(smallAbund))
  subFram <- t(smallAbund[, samps])
  subFram <- subFram[, colSums(subFram != 0) > 0]
  
  dfram <- vegdist(subFram, method = "jaccard", upper = TRUE)
  fMRow <- grep("10m", labels(dfram))
  
  if (length(fMRow) != 0){
    compRow <- as.matrix(dfram)[fMRow,]
    depths <- sapply(strsplit(names(compRow), "_"), `[[`, 2)
    
    if("DCM" %in% depths){
      smalls10[i, "DCM"] <- compRow[depths == "DCM"]
    }
    if("50m" %in% depths){
      smalls10[i, "m50"] <- compRow[depths == "50m"]
    }
    if("90m" %in% depths){
      smalls10[i, "m90"] <- compRow[depths == "90m"]
    }
  }
}

smallsDCM <- data.frame("Date" = datesSmall, "m50" = NA,  "m90" = NA)
for (i in 1:length(datesSmall)){
  samps <- grep(paste0("^", gsub("(?<![0-9])0+", "", format(datesSmall[i], "%d%b%y"), 
                                 perl = TRUE)), colnames(smallAbund))
  subFram <- t(smallAbund[, samps])
  subFram <- subFram[, colSums(subFram != 0) > 0]
  
  dfram <- vegdist(subFram, method = "jaccard", upper = TRUE)
  fMRow <- grep("DCM", labels(dfram))
  
  if (length(fMRow) != 0){
    compRow <- as.matrix(dfram)[fMRow,]
    depths <- sapply(strsplit(names(compRow), "_"), `[[`, 2)
    
    if("50m" %in% depths){
      smallsDCM[i, "m50"] <- compRow[depths == "50m"]
    }
    if("90m" %in% depths){
      smallsDCM[i, "m90"] <- compRow[depths == "90m"]
    }
  }
}

smallsDCM <- smallsDCM[smallsDCM$Date != as.Date("2017-08-07"),]

smalls50 <- data.frame("Date" = datesSmall, "m90" = NA)
for (i in 1:length(datesSmall)){
  samps <- grep(paste0("^", gsub("(?<![0-9])0+", "", format(datesSmall[i], "%d%b%y"), 
                                 perl = TRUE)), colnames(smallAbund))
  subFram <- t(smallAbund[, samps])
  subFram <- subFram[, colSums(subFram != 0) > 0]
  
  dfram <- vegdist(subFram, method = "jaccard", upper = TRUE)
  fMRow <- grep("50m", labels(dfram))
  
  if (length(fMRow) != 0){
    compRow <- as.matrix(dfram)[fMRow,]
    depths <- sapply(strsplit(names(compRow), "_"), `[[`, 2)
    
    if("90m" %in% depths){
      smalls50[i, "m90"] <- compRow[depths == "90m"]
    }
  }
}




mixStart <- as.Date(c("2016-12-15", "2017-11-29"))
mixEnd <- as.Date(c("2017-04-30", "2018-05-08"))

viridis4 <- cet_pal(7, name = "l16")[1:5]
pchs <- c(25, 18, 17, 19, 15)
ltys <- c(5, 1, 2, 3, 6)
cexSmall <- 0.9

tiff("FigBinExtras/jacMixing_allDepths.tiff", 
     width = 7, height = 8, pointsize = 12, units = "in", res = 300)
plot.new()

par(new = "TRUE",plt = c(0.09,0.95,0.79,0.99),las = 1)
plot(smalls5$Date, smalls5$DCM, xlab = "", xaxt = "n", ylim = c(0.5, 1),
     yaxt = "n", ylab = "", col = viridis4[3], lwd = 1, pch = pchs[2], 
     lty = ltys[2])
rect(mixStart[1], -120, mixEnd[1], 20, col = "#B4B4B4", border = NA)
rect(mixStart[2], -120, mixEnd[2], 20, col = "#B4B4B4", border = NA)
abline(v=as.Date("2017-01-01"), lty = 2, col = "gray75", lwd = 1)
abline(v=as.Date("2018-01-01"), lty = 2, col = "gray75", lwd = 1)
lines(smalls5$Date, smalls5$m90, col = viridis4[1], lwd = 1, 
      pch = pchs[5], lty = ltys[5])
lines(smalls5$Date, smalls5$m50, col = viridis4[2], lwd = 1, 
      pch = pchs[4], lty = ltys[4])
lines(smalls5$Date, smalls5$DCM, col = viridis4[3], lwd = 1, 
      pch = pchs[3], lty = ltys[3])
lines(smalls5$Date, smalls5$m10, col = viridis4[4], lwd = 1, 
      pch = pchs[2], lty = ltys[2])
points(smalls5$Date, smalls5$m90, col = viridis4[1], lwd = 1, 
      pch = pchs[5], lty = ltys[5])
points(smalls5$Date, smalls5$m50, col = viridis4[2], lwd = 1, 
      pch = pchs[4], lty = ltys[4])
points(smalls5$Date, smalls5$DCM, col = viridis4[3], lwd = 1, 
      pch = pchs[3], lty = ltys[3])
points(smalls5$Date, smalls5$m10, col = viridis4[4], lwd = 1, 
      pch = pchs[2], lty = ltys[2])
axis(1, lwd.ticks = 1, tck = -0.05, padj = -0.5,
     labels = NA, 
     at = as.numeric(seq(as.Date("2016-09-01"), max(smalls5$Date), 
                         by = "month"), '%m-%y')
     [seq(1, length(format(seq(as.Date("2016-09-01"), max(smalls5$Date), 
                               by = "month"), '%m-%y')), 3)], cex = cexSmall)
axis(1, lwd.ticks = 1, tck = -0.02, padj = 1, 
     at = as.numeric(seq(as.Date("2016-09-01"), max(smalls5$Date), by = "month"), 
                     '%m-%y'), labels = NA)
axis(2, tck = -0.03, lwd.ticks = 1, hadj = 0.5, cex = cexSmall)
box(lwd = 1)
text(as.Date("09-25-2016", format = "%m-%d-%Y"), 0.98, "5 m (a)")


###10 m
par(new = "TRUE",plt = c(0.09,0.95,0.57,0.77),las = 1)
plot(smalls10$Date, smalls10$DCM, xlab = "", xaxt = "n", ylim = c(0.5, 1),
     yaxt = "n", ylab = "", col = viridis4[3], lwd = 1, pch = pchs[2], 
     lty = ltys[2])
rect(mixStart[1], -120, mixEnd[1], 20, col = "#B4B4B4", border = NA)
rect(mixStart[2], -120, mixEnd[2], 20, col = "#B4B4B4", border = NA)
abline(v=as.Date("2017-01-01"), lty = 2, col = "gray75", lwd = 1)
abline(v=as.Date("2018-01-01"), lty = 2, col = "gray75", lwd = 1)
lines(smalls10$Date, smalls10$m90, col = viridis4[1], lwd = 1, 
      pch = pchs[5], lty = ltys[5])
lines(smalls10$Date, smalls10$m50, col = viridis4[2], lwd = 1, 
      pch = pchs[4], lty = ltys[4])
lines(smalls10$Date, smalls10$DCM, col = viridis4[3], lwd = 1, 
      pch = pchs[3], lty = ltys[3])
points(smalls10$Date, smalls10$m90, col = viridis4[1], lwd = 1, 
       pch = pchs[5], lty = ltys[5])
points(smalls10$Date, smalls10$m50, col = viridis4[2], lwd = 1, 
       pch = pchs[4], lty = ltys[4])
points(smalls10$Date, smalls10$DCM, col = viridis4[3], lwd = 1, 
       pch = pchs[3], lty = ltys[3])
axis(1, lwd.ticks = 1, tck = -0.05, padj = -0.5,
     labels = NA, 
     at = as.numeric(seq(as.Date("2016-09-01"), max(smalls10$Date), 
                         by = "month"), '%m-%y')
     [seq(1, length(format(seq(as.Date("2016-09-01"), max(smalls10$Date), 
                               by = "month"), '%m-%y')), 3)], cex = cexSmall)
axis(1, lwd.ticks = 1, tck = -0.02, padj = 1, 
     at = as.numeric(seq(as.Date("2016-09-01"), max(smalls10$Date), by = "month"), 
                     '%m-%y'), labels = NA)
axis(2, tck = -0.03, lwd.ticks = 1, hadj = 0.5, cex = cexSmall)
box(lwd = 1)
mtext("Jaccard distance", side = 2, line = 2, las = 0, adj = -1.3)
text(as.Date("09-27-2016", format = "%m-%d-%Y"), 0.98, "10 m (b)")



###DCM
par(new = "TRUE",plt = c(0.09,0.95,0.35,0.55),las = 1)
plot(smallsDCM$Date, smallsDCM$m50, xlab = "", xaxt = "n", ylim = c(0.5, 1),
     yaxt = "n", ylab = "", col = viridis4[2], lwd = 1, pch = pchs[4], 
     lty = ltys[4])
rect(mixStart[1], -120, mixEnd[1], 20, col = "#B4B4B4", border = NA)
rect(mixStart[2], -120, mixEnd[2], 20, col = "#B4B4B4", border = NA)
abline(v=as.Date("2017-01-01"), lty = 2, col = "gray75", lwd = 1)
abline(v=as.Date("2018-01-01"), lty = 2, col = "gray75", lwd = 1)
lines(smallsDCM$Date, smallsDCM$m90, col = viridis4[1], lwd = 1, 
      pch = pchs[5], lty = ltys[5])
lines(smallsDCM$Date, smallsDCM$m50, col = viridis4[2], lwd = 1, 
      pch = pchs[4], lty = ltys[4])
points(smallsDCM$Date, smallsDCM$m90, col = viridis4[1], lwd = 1, 
       pch = pchs[5], lty = ltys[5])
points(smallsDCM$Date, smallsDCM$m50, col = viridis4[2], lwd = 1, 
       pch = pchs[4], lty = ltys[4])
axis(1, lwd.ticks = 1, tck = -0.05, padj = -0.5,
     labels = NA, 
     at = as.numeric(seq(as.Date("2016-09-01"), max(smallsDCM$Date), 
                         by = "month"), '%m-%y')
     [seq(1, length(format(seq(as.Date("2016-09-01"), max(smallsDCM$Date), 
                               by = "month"), '%m-%y')), 3)], cex = cexSmall)
axis(1, lwd.ticks = 1, tck = -0.02, padj = 1, 
     at = as.numeric(seq(as.Date("2016-09-01"), max(smallsDCM$Date), by = "month"), 
                     '%m-%y'), labels = NA)
axis(2, tck = -0.03, lwd.ticks = 1, hadj = 0.5, cex = cexSmall)
box(lwd = 1)
text(as.Date("10-13-2016", format = "%m-%d-%Y"), 0.98, "Chl max (c)")



###50
par(new = "TRUE",plt = c(0.09,0.95,0.13,0.33),las = 1)
plot(smalls50$Date, smalls50$m90, xlab = "", xaxt = "n", ylim = c(0.5, 1),
     yaxt = "n", ylab = "", col = viridis4[1], lwd = 1, pch = pchs[5], 
     lty = ltys[5])
rect(mixStart[1], -120, mixEnd[1], 20, col = "#B4B4B4", border = NA)
rect(mixStart[2], -120, mixEnd[2], 20, col = "#B4B4B4", border = NA)
abline(v=as.Date("2017-01-01"), lty = 2, col = "gray75", lwd = 1)
abline(v=as.Date("2018-01-01"), lty = 2, col = "gray75", lwd = 1)
lines(smalls50$Date, smalls50$m90, col = viridis4[1], lwd = 1, 
      pch = pchs[5], lty = ltys[5])
points(smalls50$Date, smalls50$m90, col = viridis4[1], lwd = 1, 
       pch = pchs[5], lty = ltys[5])
axis(1, lwd.ticks = 1, tck = -0.05, padj = -0.5,
     labels = format(seq(as.Date("2016-09-01"), max(smalls50$Date), 
                         by = "month"), '%b')
     [seq(1, length(format(seq(as.Date("2016-09-01"), max(smalls50$Date), 
                               by = "month"), '%m-%y')), 3)], 
     at = as.numeric(seq(as.Date("2016-09-01"), max(smalls50$Date), 
                         by = "month"), '%m-%y')
     [seq(1, length(format(seq(as.Date("2016-09-01"), max(smalls50$Date), 
                               by = "month"), '%m-%y')), 3)], cex = cexSmall)
axis(1, lwd.ticks = 1, tck = -0.02, padj = 1, 
     at = as.numeric(seq(as.Date("2016-09-01"), max(smalls50$Date), by = "month"), 
                     '%m-%y'), labels = NA)
axis(2, tck = -0.03, lwd.ticks = 1, hadj = 0.5, cex = cexSmall)
box(lwd = 1)
text(as.Date("09-28-2016", format = "%m-%d-%Y"), 0.98, "50 m (d)")

mtext("2016", side = 1, line = 2, adj = 0.06, cex = cexSmall)
mtext("2017", side = 1, line = 2, adj = 0.375, cex = cexSmall)
mtext("2018", side = 1, line = 2, adj = 0.83, cex = cexSmall)
mtext("Date", side = 1, line = 3, adj = 0.5)
box(lwd = 1)

legend("bottom", inset=c(0,-0.7), 
       legend = c("10 m", "Chl max  ", "50 m", "90 m"), 
       col = rev(viridis4)[2:5], horiz = TRUE, pch = pchs[2:5], lty = ltys[2:5], 
       xpd = NA, bty = "n", text.width = c(60, 100, 60),
       lwd = 1, box.lwd = 0, adj = 0.2, pt.bg = viridis4[5])

dev.off()




