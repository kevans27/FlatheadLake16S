library(readxl)
library(akima)
library(viridis)
source("~/filled.contour3.R")

fullset <- read.csv('~/FlatheadPublic/hydrolab.csv', na.strings="")
fullset <- fullset[fullset$Site == 'Flathead Lake, Midlake Deep',]
units <- fullset[1,]
fullset <- fullset[-1,]
tempStart <- as.Date("2016-01-01")
tempEnd  <- as.Date("2019-01-31")

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

df <- df[df$Depth %in% c(1, 5, 10, 15, 20, 50, 90),]

tempDF <- df[df$Date >= tempStart & df$Date <= tempEnd,]
tempDF <- tempDF[!is.na(tempDF$Time),]

fldTemp <- interp(df$Month, df$Depth, df$Temp, 
                  duplicate = "mean", yo = seq(1, 90, 1), 
                  xo = seq(1, 12), extrap = TRUE)

allDatSeq <- seq(min(tempDF$Date), max(tempDF$Date), by = "day")
fldTempAllDate <- interp(tempDF$Date, tempDF$Depth, tempDF$Temp, duplicate = "mean", 
                         xo = allDatSeq,
                         yo = seq(1, 90, 1), extrap = TRUE)


DCMDat <- read.csv("~/FlatheadPublic/DCMData.csv", header = TRUE, stringsAsFactors = FALSE)
DCMDat$Date <- as.Date(DCMDat$Date)
DCMDat$Date.1 <- as.Date(DCMDat$Date.1)
DCMDat$Date.1.1 <- as.Date(DCMDat$Date.1.1)
DCMDat <- DCMDat[order(DCMDat$Date),]
DCMDat <- aggregate(DCM~Date, DCMDat, FUN = mean)

mldData <- 
  read.csv("~/FlatheadPublic/MLDData.csv", stringsAsFactors = FALSE)
mldData$Date <- as.Date(mldData$Date)
mldData <- mldData[order(mldData$Date),]

datSeq <- seq(as.Date("2016-09-01"), as.Date("2020-01-01"), by = "month")

startDate <- as.Date("2016-01-01")
endDate  <- as.Date("2019-01-31")
pchS <- c(1, 2, 0, 8, 6)

png("~/FlatheadMicrobes/figBin/tempProfile.png", width = 1600, height = 800)
plot.new()

#####Temperature
par(new = "TRUE",plt = c(0.1,0.8,0.15,0.93),las = 1, cex.axis = 2.5)
p <- filled.contour3(x = fldTempAllDate$x, bty = "n",
                     y = fldTempAllDate$y, z = fldTempAllDate$z,
                     color.palette = function(x)viridis(x), nlevels = 10, 
                     plot.title={
                       mtext("",1,line=5,las=1,cex=2.5)
                       mtext("Depth (m)",2,line=5.5, las = 0,cex=2.5, font = 2)
                     },
                     plot.axes={
                       axis(1, cex.axis=2.5, labels = NA, 
                            at = as.numeric(seq(as.Date("2016-09-01"), as.Date("2020-01-01"), 
                                                by = "month")), lwd.ticks = 3, tck = -0.015); 
                       points(DCMDat$Date, DCMDat$DCM, col = "black",
                              pch = 18, cex = 2.5); 
                       lines(mldData$Date, mldData$MLD, col = "gray", lwd = 3.5); 
                       abline(v=as.Date("2017-01-01"), lty = 2, col = "gray75", lwd = 3);
                       abline(v=as.Date("2018-01-01"), lty = 2, col = "gray75", lwd = 3);
                       abline(v=as.Date("2019-01-01"), lty = 2, col = "gray75", lwd = 3)
                       axis(2, cex.axis=2.5, hadj = 1.2, lwd.ticks = 3, tck = -0.025, 
                            font = 2);
                       axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1, font = 2,
                            at = datSeq[seq(1, length(datSeq), 3)],
                            labels = format(datSeq[seq(1, length(datSeq), 3)], "%b"))
                     },
                     key.axes={
                       axis(4, cex.axis = 1.9) #####################
                     },
                     ylim = rev(range(fldTempAllDate$y)),
                     xlim = c(as.Date("2016-09-01"), as.Date("2019-01-01"))
)
mtext("2016", 1, line = 6, font = 2, cex = 2.5, las = 1, adj = 0.05)
mtext("2017", 1, line = 6, font = 2, cex = 2.5, las = 1, adj = 0.35)
mtext("2018", 1, line = 6, font = 2, cex = 2.5, las = 1, adj = 0.8)
box(lwd = 3)

par(new = "TRUE",plt = c(0.83,0.9,0.15,0.93), cex.axis = 2.5, font = 2)
filled.legend(seq(1, 12), fldTemp$y, fldTemp$z, color = viridis, xlab = "", ylab = "",
              xlim = c(min(xintercepts),max(xintercepts)), ylim = c(min(slopes),max(slopes)),
              zlim = c(min(fldTemp$z), max(fldTemp$z)), nlevels = 10, cex.lab = 2.5, bty = "n",
              lwd.ticks = 3, tck = -0.025, font = 2)
mtext("Temperature (°C)", 3, line = 1, cex = 2.5, adj = -0, font = 2)
box(lwd = 3)

dev.off()
