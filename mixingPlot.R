library(vegan)
viridis4 <- c("#440154FF", "#33638DFF", "#20A387FF", "#95D840FF")
pchs <- c(18, 17, 19, 15)
ltys <- c(1, 2, 3, 6)


largeAbund <- read.csv("SubsampledData/rrBothLarge.count_table", stringsAsFactors = FALSE)
smallAbund <- read.csv("SubsampledData/rrBothSmall.count_table", stringsAsFactors = FALSE)
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
}


#Find dates for all
for (i in 1){
  uniqueLarge <- unique(gsub("_[^_]+$", "\\1", colnames(largeAbund)))[
    grep("_", unique(gsub("_[^_]+$", "\\1", colnames(largeAbund))))]
  uniqueSmall <- unique(gsub("_[^_]+$", "\\1", colnames(smallAbund)))[
    grep("_", unique(gsub("_[^_]+$", "\\1", colnames(smallAbund))))]
  
  datesLarge <- unique(as.Date(sapply(strsplit(uniqueLarge, "_"), `[[`, 1), 
                               format = "%d%b%y"))
  datesSmall <- unique(as.Date(sapply(strsplit(uniqueSmall, "_"), `[[`, 1), 
                               format = "%d%b%y"))
  datesLarge <- datesLarge[order(datesLarge)]
  datesSmall <- datesSmall[order(datesSmall)]
}

larges <- data.frame("Date" = datesLarge, "m10" = NA, "DCM" = NA, "m50" = NA, "m90" = NA)
smalls <- data.frame("Date" = datesSmall, "m10" = NA, "DCM" = NA, "m50" = NA, "m90" = NA)
for (i in 1:length(datesLarge)){
  samps <- grep(paste0("^", gsub("(?<![0-9])0+", "", format(datesLarge[i], "%d%b%y"), 
                               perl = TRUE)), colnames(largeAbund))
  subFram <- t(largeAbund[, samps])
  subFram <- subFram[, colSums(subFram != 0) > 0]
  
  dfram <- vegdist(subFram, method = "jaccard", upper = TRUE)
  fMRow <- grep("05m", labels(dfram))
  
  if (length(fMRow) != 0){
    compRow <- as.matrix(dfram)[fMRow,]
    depths <- sapply(strsplit(names(compRow), "_"), `[[`, 2)
    
    if("10m" %in% depths){
      larges[i, "m10"] <- compRow[depths == "10m"]
    }
    if("DCM" %in% depths){
      larges[i, "DCM"] <- compRow[depths == "DCM"]
    }
    if("50m" %in% depths){
      larges[i, "m50"] <- compRow[depths == "50m"]
    }
    if("90m" %in% depths){
      larges[i, "m90"] <- compRow[depths == "90m"]
    }
  }
}
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

#smalls <- smalls[smalls$Date != as.Date("2017-05-15"),]
larges <- larges[larges$Date != as.Date("2018-10-15"),]

mixStart <- as.Date(c("2016-12-15", "2017-11-29"))
mixEnd <- as.Date(c("2017-04-30", "2018-05-08"))

png("FigBin/jacMixingAll.png", width = 1200, height = 1000)
plot.new()

par(new = "TRUE",plt = c(0.12,0.99,0.60,0.98),las = 1, cex.axis = 2.5)
plot(smalls$Date, smalls$DCM, xlab = "", xaxt = "n", ylim = c(0, 1.1), type = "b",
     yaxt = "n", ylab = "", col = viridis4[3], cex = 2.5, lwd = 2, pch = pchs[2], 
     lty = ltys[2])
rect(mixStart[1], -1, mixEnd[1], 2, col = "#B4B4B4")
rect(mixStart[2], -1, mixEnd[2], 2, col = "#B4B4B4")
lines(smalls$Date, smalls$m10, type = "b", col = viridis4[4], cex = 2.5, lwd = 2, 
      pch = pchs[1], lty = ltys[1])
lines(smalls$Date, smalls$DCM, type = "b", col = viridis4[3], cex = 2.5, lwd = 2, 
      pch = pchs[2], lty = ltys[2])
lines(smalls$Date, smalls$m50, type = "b", col = viridis4[2], cex = 2.5, lwd = 2, 
      pch = pchs[3], lty = ltys[3])
lines(smalls$Date, smalls$m90, type = "b", col = viridis4[1], cex = 2.5, lwd = 2, 
      pch = pchs[4], lty = ltys[4])
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
abline(v=as.Date("2017-01-01"), lty = 2, col = "gray75", lwd = 3)
abline(v=as.Date("2018-01-01"), lty = 2, col = "gray75", lwd = 3)
text(min(range(larges$Date)), 0, adj = c(0, 0), labels = "Small (0.2-3 μm)", cex = 2.5, 
     font = 2)
text(min(range(larges$Date)), 1.02, adj = c(0, 0), labels = "A", cex = 2.5, font = 2)
box(lwd = 3)
mtext("Jaccard Distance", side = 2, line = 6, adj = -2.3, las = 0, font = 2, cex = 3)


par(new = "TRUE",plt = c(0.12,0.99,0.19,0.57),las = 1, cex.axis = 2.5)
plot(larges$Date, larges$DCM, xlab = "", xaxt = "n", ylim = c(0, 1.1), type = "b",
     yaxt = "n", ylab = "", col = viridis4[3], cex = 2.5, lwd = 2, pch = pchs[2], 
     lty = ltys[2])
rect(mixStart[1], -1, mixEnd[1], 2, col = "#B4B4B4")
rect(mixStart[2], -1, mixEnd[2], 2, col = "#B4B4B4")
lines(larges$Date, larges$m10, type = "b", col = viridis4[4], cex = 2.5, lwd = 2, 
      pch = pchs[1], lty = ltys[1])
lines(larges$Date, larges$DCM, type = "b", col = viridis4[3], cex = 2.5, lwd = 2, 
      pch = pchs[2], lty = ltys[2])
lines(larges$Date, larges$m50, type = "b", col = viridis4[2], cex = 2.5, lwd = 2, 
      pch = pchs[3], lty = ltys[3])
lines(larges$Date, larges$m90, type = "b", col = viridis4[1], cex = 2.5, lwd = 2, 
      pch = pchs[4], lty = ltys[4])
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
abline(v=as.Date("2017-01-01"), lty = 2, col = "gray75", lwd = 3)
abline(v=as.Date("2018-01-01"), lty = 2, col = "gray75", lwd = 3)
text(min(range(larges$Date)), 0, adj = c(0, 0), labels = "Large (>3 μm)", cex = 2.5, 
     font = 2)
text(min(range(larges$Date)), 1.02, adj = c(0, 0), labels = "B", cex = 2.5, font = 2)
mtext("2016", side = 1, line = 5.1, cex = 2.5, font = 2, adj = 0.06)
mtext("2017", side = 1, line = 5.1, cex = 2.5, font = 2, adj = 0.375)
mtext("2018", side = 1, line = 5.1, cex = 2.5, font = 2, adj = 0.83)
box(lwd = 3)

legend("bottom", inset=c(0,-0.43), legend = c("10 m", "Chl. Max  ", "50 m", "90 m"), 
       col = rev(viridis4), horiz = TRUE, pch = pchs, lty = ltys, xpd = NA, pt.cex = 2.5, 
       cex = 2.5, lwd = 2.5, box.lwd = 3, text.font = 2, adj = 0.2)

dev.off()


