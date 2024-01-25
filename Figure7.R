source("~/FlatheadMicrobes/ready5mProfilerData.R")
source("~/FlatheadMicrobes/readInWindData.R")
datadump <- "~/FlatheadPublic/Profiler_LMPDat/DescendingProfiles"
allfiles <- list.files(path = datadump, recursive = TRUE, 
                       pattern = "C-Vars.csv$", 
                       full.names = TRUE)
allfiles <- allfiles[grep("LMP46", allfiles)]

rm(fullDf)
for (i in 1:length(allfiles)){
  dfNew <- read.csv(allfiles[i], stringsAsFactors = FALSE)
  dfNew$dateTime <- as.POSIXct(dfNew$datetime_z, format = "%m/%d/%y %H:%M:%S")
  
  if(exists("fullDf")){
    fullDf <- rbind(fullDf, dfNew)
  }else{
    fullDf <- dfNew
  }
}

fullDf$Month <- as.numeric(as.character(format(fullDf$dateTime, "%m")))
fullDf$Year <- as.numeric(as.character(format(fullDf$dateTime, "%Y")))
fullDf$DOY <- as.numeric(as.character(format(fullDf$dateTime, "%j")))

library(tgp)
library(viridis)

copyFullDf <- fullDf

copyFullDf$Depth <- round(copyFullDf$depth_m)
targDepths <- c(11, 20, 50, 90)
copyFullDf <- copyFullDf[!is.na(copyFullDf$Depth),]
copyFullDf[copyFullDf$Depth == 13 | copyFullDf$Depth == 12, "Depth"] <- 11
copyFullDf[copyFullDf$Depth == 89 | copyFullDf$Depth == 88, "Depth"] <- 90
good <- copyFullDf[copyFullDf$Depth %in% targDepths,]
good <- good[good$dateTime > as.POSIXct("2011-10-01"),]
good$DTH <- format(good$dateTime, "%m/%d/%Y %H")
good$Date <- as.Date(good$dateTime)


td <- aggregate(temp_c ~ Date + Depth, good, FUN = mean)
colnames(td) <- c("Date", "Depth", "Temp")
dailyAverages5m$Temp.sd <- NULL
dailyAverages5m$Depth <- 5
colnames(dailyAverages5m) <- c("Date", "Temp", "Depth")

merged <- rbind(td, dailyAverages5m)

topPlotInterest <- merged[merged$Date > as.Date("2013-02-28") & 
                            merged$Date < as.Date("2013-06-01"),]
topRightPI <- merged[merged$Date > as.Date("2013-09-30") & 
                       merged$Date < as.Date("2014-01-01"),]

#td$Date <- as.POSIXct(td$DTH, format = "%m/%d/%Y %H")
colors <- rev(viridis(6)[1:5])

zoomMinDate <- as.Date("2013-04-01")
zoomMaxDate <- as.Date("2013-04-25")
zoomMinY <- 3.7
zoomMaxY <- 4.85

xminL <- 0.08
xmaxL <- 0.39
xminR <- 0.56
xmaxR <- 0.87

miniDf <- topPlotInterest[topPlotInterest$Date >= zoomMinDate & 
                            topPlotInterest$Date <= zoomMaxDate,]

miniWind <- fulldf[as.Date(fulldf$DateTime) >= zoomMinDate & 
                     as.Date(fulldf$DateTime) <= zoomMaxDate,]

letters <- c("a", "b", "c", "d")

miniWind$kph <- miniWind$wspeed*1.609344
#miniDf <- miniDf[miniDf$Depth == 11,]

tiff("~/FlatheadMicrobes/FigBinExtras/profilerTempsBiSeason.tiff", 
     width = 9, height = 5.5, pointsize = 12, units = "in", res = 300)
par(plt = c(xminL, xmaxL, 0.55, 0.99), las = 1)
plot(topPlotInterest$Date, topPlotInterest$Temp, 
     col = colors[as.factor(topPlotInterest$Depth)],
     xlab = "", ylab = "", yaxt = "n", xaxt = "n")
mtext("Temperature (°C)", side = 2, line = 2.6, las = 0)
axis(1, at = seq(as.Date("2013-03-01"), as.Date("2013-06-01"), by = "month"),
     labels = format(seq(as.Date("2013-03-01"), as.Date("2013-06-01"), 
                         by = "month"), "%b"))
#text(as.Date("2013-03-05"), 11.8, letters[1])
axis(2, at = seq(0, 25, by = 2.5))
box(lwd = 1)
rect(zoomMinDate, zoomMinY, zoomMaxDate, zoomMaxY, lwd = 1)

lines(c(zoomMinDate, min(topPlotInterest$Date)-3.5), c(zoomMaxY, 0.5), lty = 2, 
      lwd = 1, xpd = NA)
lines(c(zoomMaxDate, max(topPlotInterest$Date)+3.5), c(zoomMaxY, 0.5), lty = 2, 
      lwd = 1, xpd = NA)

par(new = "TRUE", plt = c(xminL, xmaxL, 0.12, 0.44), las = 1)

plot(miniWind$DateTime, miniWind$kph, xlab = "", ylab = "", yaxt = "n", xaxt = "n",
     lwd = 2, type = "l", col = "lightgray", ylim = c(0, 65))
text(as.POSIXct("2013-05-02"), 30, bquote("Wind speed (km h"^{"-1"}~")"), 
     srt = 270, xpd = TRUE)
axis(1, at = seq(zoomMinDate, zoomMaxDate, by = "day")[seq(1, 31, by = 6)],
     labels = NA, tck = -0.06)
axis(1, at = seq(zoomMinDate, zoomMaxDate, by = "day"), labels = NA, tck = -0.03)
axis(4, at = seq(0, 75, by = 15))
box(lwd = 1)

par(new = "TRUE", plt = c(xminL, xmaxL, 0.12, 0.44), las = 1)
plot(miniDf$Date, miniDf$Temp, col = colors[as.factor(miniDf$Depth)], xlab = "", 
     ylab = "", xaxt = "n", yaxt = "n", ylim = c(zoomMinY, zoomMaxY))
#text(as.Date("2013-04-02"), 4.75, letters[3])
axis(1, at = seq(zoomMinDate, zoomMaxDate, by = "day")[seq(1, 31, by = 6)],
     labels = format(as.Date(seq(zoomMinDate, zoomMaxDate, by = "day")), 
                     "%b %d")[seq(1, 31, by = 6)], tck = -0.06)
axis(1, at = seq(zoomMinDate, zoomMaxDate, by = "day"), labels = NA, 
     tck = -0.03)
axis(2, at = seq(3, 6, by = 0.5))
mtext("Temperature (°C)", side = 2, line = 2.6, las = 0)
box(lwd = 1)
mtext("Date", side = 1, line = 2.4)




zoomMinDate <- as.Date("2013-12-01")
zoomMaxDate <- as.Date("2013-12-25")
zoomMinY <- 3.7
zoomMaxY <- 7

miniDf <- topRightPI[topRightPI$Date >= zoomMinDate & 
                       topRightPI$Date <= zoomMaxDate,]

miniWind <- fulldf[as.Date(fulldf$DateTime) >= zoomMinDate & 
                     as.Date(fulldf$DateTime) <= zoomMaxDate,]
miniWind$kph <- miniWind$wspeed*1.609344

par(new = "TRUE", plt = c(xminR, xmaxR, 0.55, 0.99), las = 1)
plot(topRightPI$Date, topRightPI$Temp, col = colors[as.factor(topRightPI$Depth)],
     xlab = "", ylab = "", yaxt = "n", xaxt = "n")
mtext("Temperature (°C)", side = 2, line = 2.5, las = 0)
#text(as.Date("2013-10-05"), 13.6, letters[2])
axis(1, at = seq(as.Date("2013-10-01"), as.Date("2014-01-01"), by = "month"),
     labels = format(seq(as.Date("2013-10-01"), as.Date("2014-01-01"), 
                         by = "month"), "%b"))
axis(2, at = seq(0, 25, by = 2.5))
box(lwd = 1)
rect(zoomMinDate, zoomMinY, zoomMaxDate, zoomMaxY, lwd = 1)
lines(c(zoomMinDate, min(topRightPI$Date)-3.5), c(zoomMaxY, 0.8), lty = 2, xpd = NA)
lines(c(zoomMaxDate, max(topRightPI$Date)+3.5), c(zoomMaxY, 0.8), lty = 2, xpd = NA)
par(new = "TRUE", plt = c(xminR, xmaxR, 0.12, 0.44), las = 1)
plot(miniWind$DateTime, miniWind$kph, xlab = "", ylab = "", yaxt = "n", 
     xaxt = "n", lwd = 2, type = "l", col = "lightgray", ylim = c(0, 65))
text(as.POSIXct("2014-01-01"), 30, bquote("Wind speed (km h"^{"-1"}~")"), 
     srt = 270, xpd = TRUE)
axis(1, at = seq(zoomMinDate, zoomMaxDate, by = "day")[seq(1, 31, by = 6)],
     labels = NA, tck = -0.06)
axis(1, at = seq(zoomMinDate, zoomMaxDate, by = "day"), labels = NA, 
     tck = -0.03)
axis(4, at = seq(0, 75, by = 15))
box(lwd = 1)

par(new = "TRUE", plt = c(xminR, xmaxR, 0.12, 0.44), las = 1)
plot(miniDf$Date, miniDf$Temp, col = colors[as.factor(miniDf$Depth)], xlab = "", 
     ylab = "", xaxt = "n", yaxt = "n",
     ylim = c(zoomMinY, zoomMaxY))
#text(as.Date("2013-12-02"), 6.8, letters[4])
axis(1, at = seq(zoomMinDate, zoomMaxDate, by = "day")[seq(1, 31, by = 6)],
     labels = format(as.Date(seq(zoomMinDate, zoomMaxDate, by = "day")), 
                     "%b %d")[seq(1, 31, by = 6)], tck = -0.06)
axis(1, at = seq(zoomMinDate, zoomMaxDate, by = "day"), labels = NA, tck = -0.03)
axis(2, at = seq(3, 7, by = 1), labels = sprintf("%.1f",seq(3, 7, by = 1)))
mtext("Temperature (°C)", side = 2, line = 2.5, las = 0)
legend(as.Date("2013-12-26"), 12, legend = c(5, targDepths), 
       fill = c(colors), 
       title = "Depth (m)", xpd = NA, box.lwd = 0)
box(lwd = 1)
mtext("Date", side = 1, line = 2.4)
legend(as.Date("2013-12-25"), 9.5, legend = c("Wind speed"), 
       fill = c("lightgray"), xpd = NA, box.lwd = 0, bg = NA, x.intersp = 0.2)
dev.off()


miniWind$Date <- as.Date(miniWind$DateTime)
dailyWind <- aggregate(kph ~ Date, data = miniWind, FUN = mean)

rangeVals <- c()
for(i in unique(miniDf$Date)){
  mldCheckDf <- miniDf[miniDf$Date == i,]
  val <- max(mldCheckDf$Temp) - min(mldCheckDf$Temp)
  rangeVals <- c(rangeVals, val)
}

dailyWind$Range <- rangeVals

tiff("~/FlatheadMicrobes/FigBinExtras/windVsTempRange.tiff", 
     width = 7, height = 5.5, pointsize = 12, units = "in", res = 1200)
plot(dailyWind$kph, dailyWind$Range)
dev.off()