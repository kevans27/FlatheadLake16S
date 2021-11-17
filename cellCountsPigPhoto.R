library(readxl)    
read_excel_allsheets <- function(filename, tibble = FALSE){
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}
DCMDat <- read.csv("~/FlatheadPublic/DCMData.csv", header = TRUE, stringsAsFactors = FALSE)
DCMDat$Date <- as.Date(DCMDat$Date)
DCMDat$Date.1 <- as.Date(DCMDat$Date.1)
DCMDat$Date.1.1 <- as.Date(DCMDat$Date.1.1)
DCMDat$X <- NULL

mySheets <- read_excel_allsheets("~/FlatheadPPs/FlowBase.xlsx")
totals <- mySheets[grep("SYBR", names(mySheets))]
startDate <- as.Date("2016-09-09")
endDate  <- as.Date("2018-11-13")

dates <- as.Date(unlist(lapply(strsplit(names(totals), " "), `[[`, 1)), format = "%m.%d.%y")

relTotals <- totals[dates <= endDate & dates >= startDate]
allCounts <- NULL
for (i in seq(1, length(relTotals))){
  df <- relTotals[[i]]
  df$`Sample Name` <- gsub(" 14 2500", "", df$`Sample Name`)
  df$Depth <- gsub('m','', gsub(' ', '', gsub('I', '', 
                                              gsub('\\s*\\([^\\)]+\\)','', df$`Sample Name`))))
  date <- as.Date(unlist(lapply(strsplit(names(totals), " "), `[[`, 1)), format = "%m.%d.%y")[i]
  
  subdf <- df[grep("NA|Pop", df$`Gate Name`),]
  subdf$Concentration <- as.numeric(subdf$Concentration)
  
  agAll <- aggregate(Concentration~`Sample Name`+Depth, subdf, FUN = sum)
  agAll$Depth <- gsub("-", "", agAll$Depth)
  agAll$Depth <- gsub("SYBR", "", agAll$Depth)
  
  meanAll <- aggregate(Concentration~Depth, agAll, FUN = mean)
  
  meanAll$Date <- date
  colnames(meanAll) <- c("Depth", "CellConc", "Date")
  
  if(!exists("allCounts") || is.null(allCounts)){
    allCounts <- meanAll
  } else{
    allCounts <- rbind(allCounts, meanAll)
  }
}
allCounts[tolower(allCounts$Depth) == "chlaax", "Depth"] <- "DCM"
allCounts[tolower(allCounts$Depth) == "chla", "Depth"] <- "DCM"

allCounts[allCounts$Date==as.Date("2016-09-09") & allCounts$Depth == 18, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2016-10-26") & allCounts$Depth == 17, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2016-11-16") & allCounts$Depth == 14, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2016-12-22") & allCounts$Depth == 17, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2017-01-20") & allCounts$Depth == 15, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2017-02-21") & allCounts$Depth == 11, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2017-03-22") & allCounts$Depth == 15, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2017-04-11") & allCounts$Depth == 14, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2017-04-24") & allCounts$Depth == 14, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2017-06-07") & allCounts$Depth == 14, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2017-06-30") & allCounts$Depth == 13, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2018-05-21") & allCounts$Depth == 14, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2018-07-05") & allCounts$Depth == 14, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2018-07-25") & allCounts$Depth == 14, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2018-08-13") & allCounts$Depth == 20, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2018-09-06") & allCounts$Depth == 26, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2018-09-24") & allCounts$Depth == 20, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2018-11-13") & allCounts$Depth == 15, "Depth"] <- "DCM"


goodDepthsCounts <- allCounts[which(allCounts$Depth %in% c(5, 10, "10DCM", "DCM", 50, 90)),]

allDepths <- rbind(goodDepthsCounts, goodDepthsCounts[goodDepthsCounts$Depth == "10DCM",])
allDepths[nrow(allDepths), "Depth"] <- "DCM"
allDepths[allDepths$Depth == "10DCM", "Depth"] <- 10

reOrgDepths <- data.frame("Date" = unique(allDepths$Date))
reOrgDepths <- merge(reOrgDepths, allDepths[allDepths$Depth==5,], by = "Date", all.x = TRUE)
colnames(reOrgDepths)[colnames(reOrgDepths) == "CellConc"] <- "m5"
reOrgDepths$Depth <- NULL
reOrgDepths <- merge(reOrgDepths, allDepths[allDepths$Depth==10,], by = "Date", all.x = TRUE)
colnames(reOrgDepths)[colnames(reOrgDepths) == "CellConc"] <- "m10"
reOrgDepths$Depth <- NULL
reOrgDepths <- merge(reOrgDepths, allDepths[allDepths$Depth=="DCM",], by = "Date", 
                     all.x = TRUE)
colnames(reOrgDepths)[colnames(reOrgDepths) == "CellConc"] <- "DCM"
reOrgDepths$Depth <- NULL
reOrgDepths <- merge(reOrgDepths, allDepths[allDepths$Depth==50,], by = "Date", all.x = TRUE)
colnames(reOrgDepths)[colnames(reOrgDepths) == "CellConc"] <- "m50"
reOrgDepths$Depth <- NULL
reOrgDepths <- merge(reOrgDepths, allDepths[allDepths$Depth==90,], by = "Date", all.x = TRUE)
colnames(reOrgDepths)[colnames(reOrgDepths) == "CellConc"] <- "m90"
reOrgDepths$Depth <- NULL

totalCellCounts <- reOrgDepths



totals <- mySheets[grep("AUTO", names(mySheets))]
startDate <- as.Date("2016-09-09")
endDate  <- as.Date("2018-11-13")

dates <- as.Date(unlist(lapply(strsplit(names(totals), " "), `[[`, 1)), format = "%m.%d.%y")

relTotals <- totals[dates <= endDate & dates >= startDate]
allCounts <- NULL
for (i in seq(1, length(relTotals))){
  df <- relTotals[[i]]
  df$Depth <- gsub('m','', gsub(' ', '', gsub('I', '', 
                                              gsub('\\s*\\([^\\)]+\\)','', df$`Sample Name`))))
  date <- as.Date(unlist(lapply(strsplit(names(totals), " "), `[[`, 1)), format = "%m.%d.%y")[i]
  
  subdfAll <- df[grep("All Events", df$`Gate Name`),]
  subdfAll$Concentration <- as.numeric(subdfAll$Concentration)
  subdfNoise <- df[grep("Noise", df$`Gate Name`),]
  subdfNoise$Concentration <- as.numeric(subdfNoise$Concentration)
  
  subdfAll$Concentration <- subdfAll$Concentration-subdfNoise$Concentration
  
  agAll <- aggregate(Concentration~`Sample Name`+Depth, subdfAll, FUN = sum)
  agAll$Depth <- gsub("-", "", agAll$Depth)
  
  meanAll <- aggregate(Concentration~Depth, agAll, FUN = mean)
  
  meanAll$Date <- date
  colnames(meanAll) <- c("Depth", "CellConc", "Date")
  
  if(!exists("allCounts") || is.null(allCounts)){
    allCounts <- meanAll
  } else{
    allCounts <- rbind(allCounts, meanAll)
  }
}
allCounts[tolower(allCounts$Depth) == "chlaax", "Depth"] <- "DCM"
allCounts[tolower(allCounts$Depth) == "chla", "Depth"] <- "DCM"

allCounts[allCounts$Date==as.Date("2016-09-09") & allCounts$Depth == 18, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2016-10-26") & allCounts$Depth == 17, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2016-11-16") & allCounts$Depth == 14, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2016-12-22") & allCounts$Depth == 17, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2017-01-20") & allCounts$Depth == 15, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2017-02-21") & allCounts$Depth == 11, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2017-03-22") & allCounts$Depth == 15, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2017-04-11") & allCounts$Depth == 14, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2017-04-24") & allCounts$Depth == 14, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2017-06-07") & allCounts$Depth == 14, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2017-06-30") & allCounts$Depth == 13, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2018-05-21") & allCounts$Depth == 14, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2018-07-05") & allCounts$Depth == 14, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2018-07-25") & allCounts$Depth == 14, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2018-08-13") & allCounts$Depth == 20, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2018-09-06") & allCounts$Depth == 26, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2018-09-24") & allCounts$Depth == 20, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2018-11-13") & allCounts$Depth == 15, "Depth"] <- "DCM"


goodDepthsCounts <- allCounts[which(allCounts$Depth %in% c(5, 10, "10DCM", "DCM", 50, 90)),]

allDepths <- rbind(goodDepthsCounts, goodDepthsCounts[goodDepthsCounts$Depth == "10DCM",])
allDepths[nrow(allDepths), "Depth"] <- "DCM"
allDepths[allDepths$Depth == "10DCM", "Depth"] <- 10

reOrgDepths <- data.frame("Date" = unique(allDepths$Date))
reOrgDepths <- merge(reOrgDepths, allDepths[allDepths$Depth==5,], by = "Date", all.x = TRUE)
colnames(reOrgDepths)[colnames(reOrgDepths) == "CellConc"] <- "m5"
reOrgDepths$Depth <- NULL
reOrgDepths <- merge(reOrgDepths, allDepths[allDepths$Depth==10,], by = "Date", all.x = TRUE)
colnames(reOrgDepths)[colnames(reOrgDepths) == "CellConc"] <- "m10"
reOrgDepths$Depth <- NULL
reOrgDepths <- merge(reOrgDepths, allDepths[allDepths$Depth=="DCM",], by = "Date", 
                     all.x = TRUE)
colnames(reOrgDepths)[colnames(reOrgDepths) == "CellConc"] <- "DCM"
reOrgDepths$Depth <- NULL
reOrgDepths <- merge(reOrgDepths, allDepths[allDepths$Depth==50,], by = "Date", all.x = TRUE)
colnames(reOrgDepths)[colnames(reOrgDepths) == "CellConc"] <- "m50"
reOrgDepths$Depth <- NULL
reOrgDepths <- merge(reOrgDepths, allDepths[allDepths$Depth==90,], by = "Date", all.x = TRUE)
colnames(reOrgDepths)[colnames(reOrgDepths) == "CellConc"] <- "m90"
reOrgDepths$Depth <- NULL

photoCellCounts <- reOrgDepths
photoCellCounts <- na.omit(photoCellCounts)

nonPigCellCounts <- totalCellCounts
nonPigCellCounts[, seq(2, 6)] <- totalCellCounts[, seq(2, 6)]-photoCellCounts[, seq(2, 6)]

mixStart <- as.Date(c("2016-12-15", "2017-11-29"))
mixEnd <- as.Date(c("2017-04-30", "2018-05-08"))

png("FigBin/cellCountsPigPhoto.png", width = 1400, height = 1600)
plot.new()

###5m###
par(new = "TRUE",plt = c(0.12,0.9,0.83,0.99),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(nonPigCellCounts$Date, nonPigCellCounts$m5, ylim = c(360, 1200),
     type = "b", pch = 18, yaxt = "n", xaxt = "n", xlab = "", ylab = "",
     cex.main = 2.5, cex = 2.5, lwd = 2.5, cex.lab = 2.5,
     main = "", col = "#AA4499")
rect(mixStart[1], 0, mixEnd[1], 1500, col = "#B4B4B4")
rect(mixStart[2], 0, mixEnd[2], 1500, col = "#B4B4B4")
points(nonPigCellCounts$Date, nonPigCellCounts$m5,
       type = "b", pch = 18, cex = 2.5, lwd = 2.5, col = "#AA4499")
axis(2, cex.axis=2.5, tck = -0.05, lwd.ticks = 3, hadj = 1.2, font = 2, 
     at = seq(0, 1200, by = 400), labels = seq(0, 1200, by = 400))
axis(2, cex.axis=2.5, tck = -0.03, lwd.ticks = 3, hadj = 1.2, font = 2, 
     at = seq(200, 1400, by = 400), labels = NA)
par(new = "TRUE")
plot(photoCellCounts$Date, photoCellCounts$m5, ylim = c(0, 90), yaxt = "n",
     type = "b", xaxt = "n", xlab = "", ylab = "",
     cex.main = 2.5, cex = 2.5, lwd = 2.5, cex.lab = 2.5, 
     main = "", col = "#117733",
     xlim = c(min(nonPigCellCounts$Date, na.rm = TRUE), max(nonPigCellCounts$Date, 
                                                           na.rm = TRUE)))
box(lwd = 3)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, padj = 1,
     labels = NA, 
     at = as.numeric(seq(as.Date("2016-09-01"), max(nonPigCellCounts$Date, na.rm = TRUE), 
                         by = "month"), '%m-%y')
     [seq(1, length(format(seq(as.Date("2016-09-01"), 
                               max(nonPigCellCounts$Date, na.rm = TRUE), by = "month"), 
                           '%m-%y')), 3)], font = 2)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1, 
     at = as.numeric(seq(as.Date("2016-09-01"), 
                         max(nonPigCellCounts$Date, na.rm = TRUE), by = "month"), '%m-%y'),
     labels = NA)
axis(4, cex.axis=2.5, tck = -0.05, lwd.ticks = 3, font = 2, hadj = -0.2,
     at = seq(0, 90, by = 45), labels = seq(0, 90, by = 45))
axis(4, cex.axis=2.5, tck = -0.03, lwd.ticks = 3, hadj = 1.2, font = 2, 
     at = seq(22.5, 112.5, by = 45), labels = NA)
abline(v=as.Date("2017-01-01"), lty = 2, col = "black", lwd = 3)
abline(v=as.Date("2018-01-01"), lty = 2, col = "black", lwd = 3)
text(min(nonPigCellCounts$Date, na.rm = TRUE), 85, "5 m", adj = c(0, 1), font = 2, cex = 2.5)

###10###
par(new = "TRUE",plt = c(0.12,0.9,0.65,0.81),las = 1, cex.axis = 2.5)
plot(nonPigCellCounts$Date, nonPigCellCounts$m10, ylim = c(360, 1200),
     type = "b", pch = 18, yaxt = "n", xaxt = "n", xlab = "", ylab = "",
     cex.main = 2.5, cex = 2.5, lwd = 2.5, cex.lab = 2.5,
     main = "", col = "#AA4499")
rect(mixStart[1], 0, mixEnd[1], 1500, col = "#B4B4B4")
rect(mixStart[2], 0, mixEnd[2], 1500, col = "#B4B4B4")
points(nonPigCellCounts$Date, nonPigCellCounts$m10,
       type = "b", pch = 18, cex = 2.5, lwd = 2.5, col = "#AA4499")
axis(2, cex.axis=2.5, tck = -0.05, lwd.ticks = 3, hadj = 1.2, font = 2, 
     at = seq(0, 1200, by = 400), labels = seq(0, 1200, by = 400))
axis(2, cex.axis=2.5, tck = -0.03, lwd.ticks = 3, hadj = 1.2, font = 2, 
     at = seq(200, 1400, by = 400), labels = NA)
par(new = "TRUE")
plot(photoCellCounts$Date, photoCellCounts$m10, ylim = c(0, 90), yaxt = "n",
     type = "b", xaxt = "n", xlab = "", ylab = "",
     cex.main = 2.5, cex = 2.5, lwd = 2.5, cex.lab = 2.5, 
     main = "", col = "#117733",
     xlim = c(min(nonPigCellCounts$Date, na.rm = TRUE), max(nonPigCellCounts$Date, 
                                                           na.rm = TRUE)))
box(lwd = 3)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, padj = 1,
     labels = NA, 
     at = as.numeric(seq(as.Date("2016-09-01"), max(nonPigCellCounts$Date, na.rm = TRUE), 
                         by = "month"), '%m-%y')
     [seq(1, length(format(seq(as.Date("2016-09-01"), 
                               max(nonPigCellCounts$Date, na.rm = TRUE), by = "month"), 
                           '%m-%y')), 3)], font = 2)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1, 
     at = as.numeric(seq(as.Date("2016-09-01"), 
                         max(nonPigCellCounts$Date, na.rm = TRUE), by = "month"), '%m-%y'),
     labels = NA)
axis(4, cex.axis=2.5, tck = -0.05, lwd.ticks = 3, font = 2, hadj = -0.2,
     at = seq(0, 90, by = 45), labels = seq(0, 90, by = 45))
axis(4, cex.axis=2.5, tck = -0.03, lwd.ticks = 3, hadj = 1.2, font = 2, 
     at = seq(22.5, 112.5, by = 45), labels = NA)
abline(v=as.Date("2017-01-01"), lty = 2, col = "black", lwd = 3)
abline(v=as.Date("2018-01-01"), lty = 2, col = "black", lwd = 3)
text(min(nonPigCellCounts$Date, na.rm = TRUE), 85, "10 m", adj = c(0, 1), font = 2, cex = 2.5)


#Chl Max
par(new = "TRUE",plt = c(0.12,0.9,0.47,0.63),las = 1, cex.axis = 2.5)
plot(nonPigCellCounts$Date, nonPigCellCounts$DCM, ylim = c(360, 1200),
     type = "b", pch = 18, yaxt = "n", xaxt = "n", xlab = "", ylab = "",
     cex.main = 2.5, cex = 2.5, lwd = 2.5, cex.lab = 2.5,
     main = "", col = "#AA4499")
rect(mixStart[1], 0, mixEnd[1], 1500, col = "#B4B4B4")
rect(mixStart[2], 0, mixEnd[2], 1500, col = "#B4B4B4")
points(nonPigCellCounts$Date, nonPigCellCounts$DCM,
       type = "b", pch = 18, cex = 2.5, lwd = 2.5, col = "#AA4499")
mtext(bquote(bold("Total Cells (10"^{"3"}~"Cells mL"^{"-1"}~")")), 2, font = 2, line = 7.5, cex = 2.5, las = 0)
axis(2, cex.axis=2.5, tck = -0.05, lwd.ticks = 3, hadj = 1.2, font = 2, 
     at = seq(0, 1200, by = 400), labels = seq(0, 1200, by = 400))
axis(2, cex.axis=2.5, tck = -0.03, lwd.ticks = 3, hadj = 1.2, font = 2, 
     at = seq(200, 1400, by = 400), labels = NA)
par(new = "TRUE")
plot(photoCellCounts$Date, photoCellCounts$DCM, yaxt = "n", ylim = c(0, 180),
     type = "b", xaxt = "n", xlab = "", ylab = "",
     cex.main = 2.5, cex = 2.5, lwd = 2.5, cex.lab = 2.5, 
     main = "", col = "#117733",
     xlim = c(min(nonPigCellCounts$Date, na.rm = TRUE), max(nonPigCellCounts$Date, 
                                                           na.rm = TRUE)))
par(xpd = TRUE)
text(max(nonPigCellCounts$Date)+120, 105, 
     bquote(bold("Cyanobacteria (10"^{"3"}~"Cells mL"^{"-1"}~")")),
     font = 2, cex = 2.5, las = 0, srt = 270)
par(xpd = FALSE)
box(lwd = 3)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, padj = 1,
     labels = NA, 
     at = as.numeric(seq(as.Date("2016-09-01"), max(nonPigCellCounts$Date, na.rm = TRUE), 
                         by = "month"), '%m-%y')
     [seq(1, length(format(seq(as.Date("2016-09-01"), 
                               max(nonPigCellCounts$Date, na.rm = TRUE), by = "month"), 
                           '%m-%y')), 3)], font = 2)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1, 
     at = as.numeric(seq(as.Date("2016-09-01"), 
                         max(nonPigCellCounts$Date, na.rm = TRUE), by = "month"), '%m-%y'),
     labels = NA)
axis(4, cex.axis=2.5, tck = -0.05, lwd.ticks = 3, font = 2, hadj = -0.2,
     at = seq(0, 180, by = 90))
axis(4, cex.axis=2.5, tck = -0.03, lwd.ticks = 3, hadj = 1.2, font = 2, 
     at = seq(45, 225, by = 90), labels = NA)
abline(v=as.Date("2017-01-01"), lty = 2, col = "black", lwd = 3)
abline(v=as.Date("2018-01-01"), lty = 2, col = "black", lwd = 3)
text(min(nonPigCellCounts$Date, na.rm = TRUE), 170, "Chl. Max", adj = c(0, 1), font = 2, 
     cex = 2.5)


##50 m
par(new = "TRUE",plt = c(0.12,0.9,0.29,0.45),las = 1, cex.axis = 2.5)
plot(nonPigCellCounts$Date, nonPigCellCounts$m50, ylim = c(360, 1200),
     type = "b", pch = 18, yaxt = "n", xaxt = "n", xlab = "", ylab = "",
     cex.main = 2.5, cex = 2.5, lwd = 2.5, cex.lab = 2.5,
     main = "", col = "#AA4499")
rect(mixStart[1], 0, mixEnd[1], 1500, col = "#B4B4B4")
rect(mixStart[2], 0, mixEnd[2], 1500, col = "#B4B4B4")
points(nonPigCellCounts$Date, nonPigCellCounts$m50,
       type = "b", pch = 18, cex = 2.5, lwd = 2.5, col = "#AA4499")
axis(2, cex.axis=2.5, tck = -0.05, lwd.ticks = 3, hadj = 1.2, font = 2, 
     at = seq(0, 1200, by = 400), labels = seq(0, 1200, by = 400))
axis(2, cex.axis=2.5, tck = -0.03, lwd.ticks = 3, hadj = 1.2, font = 2, 
     at = seq(200, 1400, by = 400), labels = NA)
par(new = "TRUE")
plot(photoCellCounts$Date, photoCellCounts$m50, ylim = c(0, 90), yaxt = "n",
     type = "b", xaxt = "n", xlab = "", ylab = "",
     cex.main = 2.5, cex = 2.5, lwd = 2.5, cex.lab = 2.5, 
     main = "", col = "#117733",
     xlim = c(min(nonPigCellCounts$Date, na.rm = TRUE), max(nonPigCellCounts$Date, 
                                                           na.rm = TRUE)))
box(lwd = 3)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, padj = 1,
     labels = NA, 
     at = as.numeric(seq(as.Date("2016-09-01"), max(nonPigCellCounts$Date, na.rm = TRUE), 
                         by = "month"), '%m-%y')
     [seq(1, length(format(seq(as.Date("2016-09-01"), 
                               max(nonPigCellCounts$Date, na.rm = TRUE), by = "month"), 
                           '%m-%y')), 3)], font = 2)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1, 
     at = as.numeric(seq(as.Date("2016-09-01"), 
                         max(nonPigCellCounts$Date, na.rm = TRUE), by = "month"), '%m-%y'),
     labels = NA)
axis(4, cex.axis=2.5, tck = -0.05, lwd.ticks = 3, font = 2, hadj = -0.2,
     at = seq(0, 90, by = 45), labels = seq(0, 90, by = 45))
axis(4, cex.axis=2.5, tck = -0.03, lwd.ticks = 3, hadj = 1.2, font = 2, 
     at = seq(22.5, 112.5, by = 45), labels = NA)
abline(v=as.Date("2017-01-01"), lty = 2, col = "black", lwd = 3)
abline(v=as.Date("2018-01-01"), lty = 2, col = "black", lwd = 3)
text(min(nonPigCellCounts$Date, na.rm = TRUE), 85, "50 m", adj = c(0, 1), font = 2, 
     cex = 2.5)


#90
par(new = "TRUE",plt = c(0.12,0.9,0.11,0.27),las = 1, cex.axis = 2.5)
plot(nonPigCellCounts$Date, nonPigCellCounts$m90, ylim = c(360, 1200),
     type = "b", pch = 18, yaxt = "n", xaxt = "n", xlab = "", ylab = "",
     cex.main = 2.5, cex = 2.5, lwd = 2.5, cex.lab = 2.5,
     main = "", col = "#AA4499")
rect(mixStart[1], 0, mixEnd[1], 1500, col = "#B4B4B4")
rect(mixStart[2], 0, mixEnd[2], 1500, col = "#B4B4B4")
points(nonPigCellCounts$Date, nonPigCellCounts$m90,
       type = "b", pch = 18, cex = 2.5, lwd = 2.5, col = "#AA4499")
axis(2, cex.axis=2.5, tck = -0.05, lwd.ticks = 3, hadj = 1.2, font = 2, 
     at = seq(0, 1200, by = 400), labels = seq(0, 1200, by = 400))
axis(2, cex.axis=2.5, tck = -0.03, lwd.ticks = 3, hadj = 1.2, font = 2, 
     at = seq(200, 1400, by = 400), labels = NA)
par(new = "TRUE")
plot(photoCellCounts$Date, photoCellCounts$m90, ylim = c(0, 90), yaxt = "n",
     type = "b", xaxt = "n", xlab = "", ylab = "",
     cex.main = 2.5, cex = 2.5, lwd = 2.5, cex.lab = 2.5, 
     main = "", col = "#117733",
     xlim = c(min(nonPigCellCounts$Date, na.rm = TRUE), max(nonPigCellCounts$Date, 
                                                           na.rm = TRUE)))
box(lwd = 3)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, padj = 1,
     labels = format(seq(as.Date("2016-09-01"), max(nonPigCellCounts$Date, na.rm = TRUE), 
                         by = "month"), '%b')
     [seq(1, length(format(seq(as.Date("2016-09-01"), 
                               max(nonPigCellCounts$Date, na.rm = TRUE), by = "month"), 
                           '%m-%y')), 3)], 
     at = as.numeric(seq(as.Date("2016-09-01"), max(nonPigCellCounts$Date, na.rm = TRUE), 
                         by = "month"), '%m-%y')
     [seq(1, length(format(seq(as.Date("2016-09-01"), 
                               max(nonPigCellCounts$Date, na.rm = TRUE), by = "month"), 
                           '%m-%y')), 3)], font = 2)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1, 
     at = as.numeric(seq(as.Date("2016-09-01"), 
                         max(nonPigCellCounts$Date, na.rm = TRUE), by = "month"), '%m-%y'),
     labels = NA)
axis(4, cex.axis=2.5, tck = -0.05, lwd.ticks = 3, font = 2, hadj = -0.2,
     at = seq(0, 90, by = 45), labels = seq(0, 90, by = 45))
axis(4, cex.axis=2.5, tck = -0.03, lwd.ticks = 3, hadj = 1.2, font = 2, 
     at = seq(22.5, 112.5, by = 45), labels = NA)
abline(v=as.Date("2017-01-01"), lty = 2, col = "black", lwd = 3)
abline(v=as.Date("2018-01-01"), lty = 2, col = "black", lwd = 3)
text(min(nonPigCellCounts$Date, na.rm = TRUE), 85, "90 m", adj = c(0, 1), font = 2, cex = 2.5)
mtext("2016", side = 1, line = 5, cex = 2.5, font = 2, adj = 0.06)
mtext("2017", side = 1, line = 5, cex = 2.5, font = 2, adj = 0.375)
mtext("2018", side = 1, line = 5, cex = 2.5, font = 2, adj = 0.83)



par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

legend("bottom", 
       legend = c("Non-pigmented Cells", "Photosynthetic Cells"), 
       text.font = 2,
       fill = c("#AA4499", "#117733"), horiz = TRUE, cex = 2.5, box.lwd = 0)

dev.off()




#Nonpig
#5 - 375 to 987
#10 - 438 to 1106
#Chl - 438 - 958
#50 - 420 - 690
#90 - 439 - 768
#
#Photo
#5 - 16 - 68
#10 - 18 - 80
#Chl - 17 - 179
#50 - 4 - 62
#90 - 2 - 59