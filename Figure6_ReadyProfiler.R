df0912 <- 
  read.csv("~/FlatheadPublic/Profiler_LMPDat/Stationary5m/MLD_5m_091212_Trim.csv",
           stringsAsFactors = FALSE)
df0912 <- df0912[!is.na(df0912$Temp...C),]
df0912$DTH <- as.POSIXct(df0912$Date.Time..GMT.07.00, format = "%m/%d/%y %H:%M:%S")
df0912$X. <- NULL
colnames(df0912) <- c("OldDTH", "Temp", "DTH")
df0912$OldDTH <- NULL
df1213 <- 
  read.csv("~/FlatheadPublic/Profiler_LMPDat/Stationary5m/MLD_5m_121013_Trim.csv",
           stringsAsFactors = FALSE)
df1213 <- df1213[!is.na(df1213$Temperature...C.),]
df1213$DTH <- as.POSIXct(paste(df1213$Date.yyyy.mm.dd., df1213$Time.hh.mm.ss.),
                         format = "%m/%d/%Y %H:%M:%S")
colnames(df1213) <- c("Date", "Time", "Temp", "DTH")
df1213$Date <- NULL
df1213$Time <- NULL
df1014 <- 
  read.csv("~/FlatheadPublic/Profiler_LMPDat/Stationary5m/MLD_5m_100314_Trim.csv",
           stringsAsFactors = FALSE)
df1014 <- df1014[!is.na(df1014$Temperature...C.),]
df1014$DTH <- as.POSIXct(paste(df1014$Date.yyyy.mm.dd., df1014$Time.hh.mm.ss.),
                         format = "%Y-%m-%d %H:%M:%S")
colnames(df1014) <- c("Date", "Time", "Temp", "DTH")
df1014$Date <- NULL
df1014$Time <- NULL

#Suspicious jump in temp during the last 48 hours. Guessing someone didn't turn the sensor off after retrieving it
df1014 <- df1014[c(1:101, 125:3411),]
df1213 <- df1213[c(1:4532, 4632:5125, 5160:7509),]
df0912 <- df0912[c(1:40, 52:683, 694:2424),]

fullSensor5m0 <- rbind(df0912, df1213)
fullSensor5m <- rbind(fullSensor5m0, df1014)
fullSensor5m$Date <- as.Date(fullSensor5m$DTH)

dayCounts <- table(fullSensor5m$Date)
#plot(fullSensor5m$DTH, fullSensor5m$Temp, xaxt = "n")
#axis(1, at = seq(as.POSIXct("2012-01-01"), as.POSIXct("2015-01-01"), by = "month"),
#     labels = as.Date(seq(as.POSIXct("2012-01-01"), as.POSIXct("2015-01-01"), 
#     by = "month")))

eventsTable <- table(fullSensor5m$Date)
lowDates <- as.Date(names(eventsTable[eventsTable < 24]))

fullSensor5m <- fullSensor5m[!(fullSensor5m$Date %in% lowDates),]

dailyAverages5m <- 
  do.call(data.frame, aggregate(Temp ~ Date, fullSensor5m, function(x) 
    c(mean = mean(x), sd = sd(x))))

#plot(dailyAverages5m$Date, dailyAverages5m$Temp.mean)
rm(list=setdiff(ls(), "dailyAverages5m"))
