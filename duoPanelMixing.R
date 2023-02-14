library(readxl)
library(akima)
library(viridis)

mainCol <- "#DC3220"
col2 <- "#005AB5"

fullset <- read.csv('~/FlatheadPublic/hydrolab.csv', na.strings="")
fullset <- fullset[fullset$Site == 'Flathead Lake, Midlake Deep',]
units <- fullset[1,]
fullset <- fullset[-1,]
tempStart <- as.Date("2016-01-01")
tempEnd  <- as.Date("2019-12-31")

DCMDat <- read.csv("~/FlatheadPublic/DCMData.csv", header = TRUE, stringsAsFactors = FALSE)
DCMDat$Date <- as.Date(DCMDat$Date)
DCMDat$Date.1 <- as.Date(DCMDat$Date.1)
DCMDat$Date.1.1 <- as.Date(DCMDat$Date.1.1)
DCMDat <- DCMDat[order(DCMDat$Date),]
DCMDat <- aggregate(DCM~Date, DCMDat, FUN = mean)

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

#df <- df[df$Depth %in% c(1, 5, 10, 15, 20, 50, 90),]

tempDF <- df[df$Date >= tempStart & df$Date <= tempEnd,]
tempDF <- tempDF[!is.na(tempDF$Time),]

fldTemp <- interp(df$Month, df$Depth, df$Temp, 
                  duplicate = "mean", yo = seq(1, 90, 1), 
                  xo = seq(1, 12), extrap = TRUE, jitter = TRUE)

allDatSeq <- seq(min(tempDF$Date), max(tempDF$Date), by = "day")
fldTempAllDate <- interp(tempDF$Date, tempDF$Depth, tempDF$Temp, duplicate = "mean", 
                         xo = allDatSeq, yo = seq(1, 90, 1), extrap = TRUE, 
                         jitter = TRUE)

# library(mgcv)

tempDF$Month <- as.numeric(format(tempDF$Date, '%m'))
tempDF$Year <- as.numeric(format(tempDF$Date, '%Y'))
tempDF$Day <- as.numeric(format(tempDF$Date, '%d'))
tempDF$Date <- as.numeric(tempDF$Date)
tempDF <- tempDF[order(tempDF$Date),]

# testTempAllDate <- gam(tempDF$Temp ~ te(tempDF$Date, tempDF$Depth))
# persp(tempDF$Date, tempDF$Depth, tempDF$Temp)
# vis.gam(testTempAllDate)
# 
library(tgp)

fldTempAllDate <- interp.loess(tempDF$Date, tempDF$Depth, tempDF$Temp, 
                               gridlen = c(400,400), span = 0.12)

p <- filled.contour(x = fldTempAllDate$x, bty = "n",
                    y = fldTempAllDate$y, z = fldTempAllDate$z,
                    color.palette = function(x)viridis(x), nlevels = 10, 
                    plot.title={
                      mtext("",1,line=5,las=1,cex=2.5)
                      mtext("Depth (m)",2,line=5.5, las = 0,cex=2.5, font = 2)
                    },
                    plot.axes={
                      axis(1, cex.axis=2.5, labels = NA, 
                           at = as.numeric(seq(as.Date("2016-09-01"), 
                                               as.Date("2018-12-31"), 
                                               by = "month")), lwd.ticks = 3, 
                           tck = -0.015);
                      lines(mldData$Date, mldData$td05, col = "gray", lwd = 5); 
                      abline(v=as.Date("2017-01-01"), lty = 2, col = "gray75", 
                             lwd = 3);
                      abline(v=as.Date("2018-01-01"), lty = 2, col = "gray75", 
                             lwd = 3);
                      abline(v=as.Date("2019-01-01"), lty = 2, col = "gray75", 
                             lwd = 5)
                      axis(2, cex.axis=2.5, hadj = 1.2, lwd.ticks = 3, 
                           tck = -0.025, font = 2);
                      axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, 
                           padj = 1, font = 2,
                           at = datSeq[seq(1, length(datSeq), 3)],
                           labels = format(datSeq[seq(1, length(datSeq), 3)], "%b"));
                      box(lwd = 3)
                    },
                    key.axes={
                      axis(4, cex.axis=2.5) #####################
                    },
                    ylim = rev(range(fldTempAllDate$y)),
                    xlim = c(as.Date("2016-09-01"), as.Date("2018-12-31"))
)




mldData <- read.csv("~/TMNT/MLDAll.csv", stringsAsFactors = FALSE)
mldData$Date <- as.Date(mldData$Date)
mldData <- mldData[order(mldData$Date),]

datSeq <- seq(as.Date("2016-09-01"), as.Date("2018-12-31"), by = "month")

startDate <- as.Date("2016-09-01")
endDate  <- as.Date("2018-12-31")
pchS <- c(1, 2, 0, 8, 6)


kVals <- read.csv("~/FlatheadPublic/lightAttenuation_Jan2023_lm.csv", 
                  stringsAsFactors = FALSE)
kVals$Date <- as.Date(kVals$Date)
kValsSub <- kVals[kVals$Date >= tempStart & kVals$Date <= tempEnd,]

parDat <- read_excel("~/FlatheadPublic/RoundingWOFirst_FMPPar2011_2019.xls")
parDat$Date <- as.Date(parDat$Date)

flbsHill <- read.csv("~/FlatheadPPs/FLBSHill_PARTotal_full.csv", stringsAsFactors = FALSE)
flbsHill$Date <- as.Date.POSIXct(as.POSIXct(flbsHill$timestamp, 
                                            format = "%m/%d/%Y %H:%M"))
goodDates <- names(table(flbsHill$Date)[table(flbsHill$Date) == 96])
dailyPARHill <- data.frame("Date" = as.Date(unique(flbsHill$Date)))
summedStuff <- aggregate(parameterValue~Date, flbsHill, FUN = sum)
colnames(summedStuff)[2] <- "dailyPAR"
dailyPARHill <- merge(dailyPARHill, summedStuff, by = "Date", all.x = TRUE)
dailyPARHill <- merge(dailyPARHill, kVals, by = "Date", all.x = TRUE)
approx(dailyPARHill$kVal)


lightDat <- data.frame("Date" = rep(unique(kValsSub$Date), each = 181), 
                       "Depth" = seq(0, 90, by = 0.5), "Light" = NA)
lightDat$Light <- unlist(parDat[match(lightDat$Date, parDat$Date), "DailySum"]*
                           exp(-kValsSub[match(lightDat$Date, kValsSub$Date), "kVal"]*lightDat$Depth))

lightDat <- lightDat[grep(as.Date("2018-08-13"), lightDat$Date, invert = TRUE),]
lightDat <- na.omit(lightDat)

allDatSeq <- seq(min(lightDat$Date), max(lightDat$Date), by = "day")
fldLightAllDate <- interp(lightDat$Date, lightDat$Depth, 
                          lightDat$Light, jitter = TRUE,
                          duplicate = "mean", 
                          xo = allDatSeq,
                          yo = seq(0, 90, 0.5), extrap = TRUE)



fldLightMiniDat <- interp(lightDat$Date, lightDat$Depth, 
                          lightDat$Light, jitter = TRUE,
                          duplicate = "mean", 
                          xo = unique(lightDat$Date),
                          yo = seq(0, 90, 0.5), extrap = TRUE)


MLDPar <- read.csv("~/FlatheadPublic/NPPIncubationPARRatios_2011-2022.csv", 
                   stringsAsFactors = FALSE)

dailyTotals <- data.frame("Date" = as.Date(MLDPar$Date), 
                          "DailyPar" = MLDPar$DailySum)
subTotals <- dailyTotals[dailyTotals$Date %in% unique(mldData$Date),]
subTotals$DailyPar <- subTotals$DailyPar/1000
subTotals <- subTotals[subTotals$Date != as.Date("2018-08-13"),]

png("FigBinExtras/surfaceLightMLD.png", width = 1600, height = 492)
plot.new()
par(mar = c(2, 4.5, 4, 8))
plot(subTotals$Date, subTotals$DailyPar, xaxt = "n", ylab = "", yaxt = "n", 
     xlab = "", pch = 18, cex = 4, lwd = 2.5, col = mainCol, ylim = c(0, 63), 
     type = "b", xlim = c(as.Date("2016-09-01"), as.Date("2018-12-31")))
box(lwd = 2.5)
abline(v=as.Date("2017-01-01"), lty = 2, col = "gray75", lwd = 3)
abline(v=as.Date("2018-01-01"), lty = 2, col = "gray75", lwd = 3)
abline(v=as.Date("2019-01-01"), lty = 2, col = "gray75", lwd = 3)

axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.07, padj = 1, font = 2,
     at = datSeq[seq(1, length(datSeq), 3)],
     labels = NA)
axis(1, cex.axis=2.5, labels = NA, 
     at = as.numeric(seq(as.Date("2016-09-01"), as.Date("2018-12-31"), 
                         by = "month")), lwd.ticks = 3, tck = -0.04)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.025, hadj = 1, 
     at = seq(10, 70, by = 20), labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, hadj = 1.2, las = 1,
     at = seq(0, 60, by = 20), labels = seq(0, 60, by = 20), font = 2)
mtext("A", font = 2, side = 3, line = 1, cex = 2.5, adj = 0)
dev.off()



###Temp
png("FigBinExtras/tempProfileSub.png", width = 1600, height = 800)
plot.new()

#####Temperature
p <- filled.contour(x = fldTempAllDate$x, bty = "n",
                    y = fldTempAllDate$y, z = fldTempAllDate$z,
                    color.palette = function(x)viridis(x), nlevels = 10, 
                    plot.title={
                      mtext("",1,line=5,las=1,cex=2.5)
                      mtext("Depth (m)",2,line=5.5, las = 0,cex=2.5, font = 2)
                    },
                    plot.axes={
                      axis(1, cex.axis=2.5, labels = NA, 
                           at = as.numeric(seq(as.Date("2016-09-01"), 
                                               as.Date("2018-12-31"), 
                                               by = "month")), lwd.ticks = 3, 
                           tck = -0.015);
                      lines(mldData$Date, mldData$td05, col = "gray", lwd = 5); 
                      points(DCMDat$Date, DCMDat$DCM, col = "#FFF8EF",
                             pch = 18, cex = 2.5); 
                      abline(v=as.Date("2017-01-01"), lty = 2, col = "gray75", 
                             lwd = 3);
                      abline(v=as.Date("2018-01-01"), lty = 2, col = "gray75", 
                             lwd = 3);
                      abline(v=as.Date("2019-01-01"), lty = 2, col = "gray75", 
                             lwd = 5)
                      axis(2, cex.axis=2.5, hadj = 1.2, lwd.ticks = 3, 
                           tck = -0.025, font = 2, at = seq(80, 0, by = -20));
                      axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, 
                           padj = 1, font = 2,
                           at = datSeq[seq(1, length(datSeq), 3)],
                           labels = format(datSeq[seq(1, length(datSeq), 3)], "%b"));
                      box(lwd = 3)
                    },
                    key.axes={
                      axis(4, cex.axis=2.5) #####################
                    },
                    ylim = rev(range(fldTempAllDate$y)),
                    xlim = c(as.Date("2016-09-01"), as.Date("2018-12-31"))
)
mtext("B", 3, line = 1, font = 2, adj = 0, cex = 2.5)


dev.off()

library(png)
mainCol <- "#DC3220"
col2 <- "#005AB5"


png("FigBinExtras/SurfLightandTemp_Intro.png", width = 1600, height = 1200)

par(mar = c(1, 3, 1, 0))
plot(0:2, 0:2, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
rasterImage(readPNG(source="FigBinExtras/surfaceLightMLD.png"), 0, 1.3, 2, 2.1)
rasterImage(readPNG(source="FigBinExtras/tempProfileSub.png"), 0, -0.025, 2, 1.3)
mtext("Depth (m)", 2, font = 2, cex = 2.5, line = -2, adj = 0.33)
mtext(expression(bold("PAR")), 
      2, font = 2, cex = 2.5, line = -1, adj = 0.84)
mtext(expression(bold("(mol quanta m"^{"-2"}~"d"^{"-1"}~")")), 
      2, font = 2, cex = 2.5, line = -4, adj = 0.92)
text(1.93, 1.25, "Temp (°C)", font = 2, cex=2.5)
mtext("2016", 1, font = 2, cex = 2.5, line = -1, adj = 0.12)
mtext("2017", 1, font = 2, cex = 2.5, line = -1, adj = 0.35)
mtext("2018", 1, font = 2, cex = 2.5, line = -1, adj = 0.74)
dev.off()

