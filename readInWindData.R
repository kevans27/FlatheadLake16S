windspeed <- 
  read.csv("~/FlatheadPublic/Weather/YBPoint_WindSpeed_Spring2013.csv")
winddirection <- 
  read.csv("~/FlatheadPublic/Weather/YBPoint_WindDirection_Spring2013.csv")

windspeed$id <- NULL
windspeed$parameterID <- NULL
windspeed$parameterName <- NULL
windspeed$parameterUnits <- NULL ##mph
winddirection$id <- NULL
winddirection$parameterID <- NULL
winddirection$parameterName <- NULL
winddirection$parameterUnits <- NULL ##degrees
colnames(winddirection) <- c("timestamp", "wdir")
colnames(windspeed) <- c("timestamp", "wspeed")

fulldf <- merge(windspeed, winddirection, by = "timestamp")

fulldf$DateTime <- as.POSIXct(fulldf$timestamp, format = "%m/%d/%Y %H:%M")
fulldf <- fulldf[order(fulldf$DateTime),]

#plot(fulldf$DateTime, fulldf$wspeed)
