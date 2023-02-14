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


#td$Date <- as.POSIXct(td$DTH, format = "%m/%d/%Y %H")
colors <- rev(viridis(6)[1:5])

zoomMinDate <- as.Date("2013-04-01")
zoomMaxDate <- as.Date("2013-04-25")
zoomMinY <- 3.7
zoomMaxY <- 4.85

xmin <- 0.15
xmax <- 0.8

miniDf <- topPlotInterest[topPlotInterest$Date >= zoomMinDate & 
                            topPlotInterest$Date <= zoomMaxDate,]

miniWind <- fulldf[as.Date(fulldf$DateTime) >= zoomMinDate & 
                     as.Date(fulldf$DateTime) <= zoomMaxDate,]

miniWind$kph <- miniWind$wspeed*1.609344
#miniDf <- miniDf[miniDf$Depth == 11,]

png("~/FlatheadMicrobes/FigBinExtras/profilerTempsFiveDepthsSpring_wind.png", 
    width = 600, height = 600)
par(plt = c(xmin, xmax, 0.55, 0.99), las = 1, cex = 1.5)
plot(topPlotInterest$Date, topPlotInterest$Temp, 
     col = colors[as.factor(topPlotInterest$Depth)],
     xlab = "", ylab = "", yaxt = "n", xaxt = "n", cex = 1, lwd = 2)
mtext("Temperature (°C)", side = 2, line = 2.5, cex = 2, font = 2, las = 0,
      adj = 0.5)
axis(1, at = seq(as.Date("2013-03-01"), as.Date("2013-06-01"), by = "month"),
     labels = format(seq(as.Date("2013-03-01"), as.Date("2013-06-01"), 
                         by = "month"), "%b"),
     font = 2, lwd.ticks = 3)
axis(2, at = seq(0, 25, by = 2.5), font = 2, lwd.ticks = 3)
legend(as.Date("2013-06-05"), 11, legend = c(5, targDepths), fill = colors, 
       title = "Depth (m)", xpd = NA, box.lwd = 0)
box(lwd = 3)
rect(zoomMinDate, zoomMinY, zoomMaxDate, zoomMaxY, lwd = 3)

lines(c(zoomMinDate, min(topPlotInterest$Date)-3.5), c(zoomMaxY, 0.5), lty = 2, 
      lwd = 3, xpd = NA)
lines(c(zoomMaxDate, max(topPlotInterest$Date)+3.5), c(zoomMaxY, 0.5), lty = 2, 
      lwd = 3, xpd = NA)

par(new = "TRUE", plt = c(xmin, xmax, 0.12, 0.44), las = 1, cex = 1.5)

plot(miniWind$DateTime, miniWind$kph, 
     xlab = "", ylab = "", yaxt = "n", xaxt = "n", cex = 1, lwd = 2, type = "l",
     col = "lightgray")
text(as.POSIXct("2013-05-01"), 30, bquote(bold("Wind speed (km h"^{"-1"}~")")), 
     cex = 1.3, font = 2, srt = 270, xpd = TRUE)
axis(1, at = seq(zoomMinDate, zoomMaxDate, by = "day")[seq(1, 31, by = 6)],
     labels = NA, font = 2, lwd.ticks = 3, tck = -0.06)
axis(1, at = seq(zoomMinDate, zoomMaxDate, by = "day"),
     labels = NA, font = 2, lwd.ticks = 3, tck = -0.03)
axis(4, at = seq(0, 75, by = 15), font = 2, lwd.ticks = 3)
box(lwd = 3)

par(new = "TRUE", plt = c(xmin, xmax, 0.12, 0.44), las = 1, cex = 1.5)
plot(miniDf$Date, miniDf$Temp, col = colors[as.factor(miniDf$Depth)],
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", lwd = 2.5,
     ylim = c(zoomMinY, zoomMaxY))
mtext("Temperature (°C)", side = 2, line = 2.5, cex = 2, font = 2, las = 0,
      adj = -1)
axis(1, at = seq(zoomMinDate, zoomMaxDate, by = "day")[seq(1, 31, by = 6)],
     labels = 
       format(as.Date(seq(zoomMinDate, zoomMaxDate, by = "day")), 
              "%b %d")[seq(1, 31, by = 6)], 
     font = 2, lwd.ticks = 3, tck = -0.06)
axis(1, at = seq(zoomMinDate, zoomMaxDate, by = "day"),
     labels = NA, font = 2, lwd.ticks = 3, tck = -0.03)
axis(2, at = seq(3, 6, by = 0.5), font = 2, lwd.ticks = 3)
legend(as.Date("2013-06-05"), 11, legend = c(5, targDepths), fill = colors, 
       title = "Depth (m)", xpd = NA, box.lwd = 0)
box(lwd = 3)
mtext("Date", side = 1, line = 2.4, cex = 2, font = 2)

dev.off()






##High wind events

hw <- miniWind[miniWind$kph > 30,]
plot(hw$kph, hw$wdir, xlab = "Wind speed (kph)", ylab = "Wind direction (Deg)")
plot(miniWind$DateTime, miniWind$wdir, xlab = "Date", ylab = "Wind direction (Deg)")
