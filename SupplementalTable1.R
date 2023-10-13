library(readxl)

fullset <- read.csv('~/FlatheadPublic/hydrolab.csv', na.strings="")
fullset <- fullset[fullset$Site == 'Flathead Lake, Midlake Deep',]
fullset <- fullset[14000:nrow(fullset),]
fullset$Date <- as.Date(fullset$Date, format = "%m/%d/%y")
units <- fullset[1,]
fullset <- fullset[-1,]
tempStart <- as.Date("2016-09-14")
tempEnd  <- as.Date("2018-12-01")

fullset <- fullset[fullset$Date <= tempEnd & fullset$Date >= tempStart,]
fullset$Month <- as.numeric(as.character(format(fullset$Date, "%m")))
fullset$D.O..1 <- as.numeric(as.character(fullset$D.O..1))
fullset <- fullset[!is.na(fullset$D.O..1),]
fullset$Depth <- round(as.numeric(as.character(fullset$Depth)))
DOMonthly <- do.call(data.frame, 
                      aggregate(D.O..1 ~ Month + Depth, fullset, 
                                function(x) c(mean = mean(x), sd = sd(x))))

##DCM
DCMDat <- read.csv("~/FlatheadPublic/DCMData.csv", header = TRUE, stringsAsFactors = FALSE)
DCMDat$Date <- as.Date(DCMDat$Date)
DCMDat$Date.1 <- as.Date(DCMDat$Date.1)
DCMDat$Date.1.1 <- as.Date(DCMDat$Date.1.1)
DCMDat <- DCMDat[order(DCMDat$Date),]
DCMDat <- aggregate(DCM~Date, DCMDat, FUN = mean)
DCMDat <- DCMDat[DCMDat$Date > tempStart & DCMDat$Date < tempEnd,]
DCMDat$Month <- format(DCMDat$Date, "%m")
DCMMonthly <- do.call(data.frame, 
                      aggregate(DCM ~ Month, DCMDat, 
                                function(x) c(mean = mean(x), sd = sd(x))))


##MLD
mldData <- read.csv("~/TMNT/MLDAll.csv", stringsAsFactors = FALSE)
mldData$Date <- as.Date(mldData$Date)
mldData <- mldData[order(mldData$Date),]
mldData <- mldData[mldData$Date > tempStart & mldData$Date < tempEnd,]
mldData$Month <- format(mldData$Date, "%m")
mldMonthly <- do.call(data.frame, 
                      aggregate(td05 ~ Month, mldData, 
                                function(x) c(mean = mean(x), sd = sd(x))))


##Light/DCM light
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

lightDat <- data.frame("Date" = rep(unique(kValsSub$Date), each = 181), 
                       "Depth" = seq(0, 90, by = 0.5), "Light" = NA)
lightDat$Light <- unlist(parDat[match(lightDat$Date, parDat$Date), "DailySum"]*
                           exp(-kValsSub[match(lightDat$Date, kValsSub$Date), 
                                         "kVal"]*lightDat$Depth))

lightDat <- lightDat[grep(as.Date("2018-08-13"), lightDat$Date, invert = TRUE),]
lightDat <- na.omit(lightDat)

DCMDat$Round <- round(DCMDat$DCM/0.5)*0.5
mergedLightDCM <- merge(lightDat, DCMDat, by.x = c("Date", "Depth"), 
                        by.y = c("Date", "Round"))
mergedLightDCM$Month <- format(mergedLightDCM$Date, "%m")

monthlyLight <- do.call(data.frame, aggregate(Light/1000~Month, mergedLightDCM, 
                                              function(x) c(mean = mean(x), 
                                                            sd = sd(x))))

##Chl
chlaDat <- read.csv("~/FlatheadPublic/FMPchla0-30.csv")

chlaDat$CollectDate <- as.Date(chlaDat$CollectDate, format = "%m/%d/%Y %H:%M:%S")
chlaDat <- chlaDat[chlaDat$Param == "Chl-a",]
chlaDat <- chlaDat[, c("CollectDate", "CorrectedReportedResult")]
subChlaDat <- chlaDat[chlaDat$CollectDate %in% unique(mldData$Date),]
subChlaDat$Month <- format(subChlaDat$CollectDate, "%m")
chlaDatmean <- do.call(data.frame, 
                       aggregate(CorrectedReportedResult~Month, subChlaDat, 
                                 function(x) c(mean = mean(x), sd = sd(x))))


##N+N
fullset <- read.csv('~/FlatheadPublic/CNPNutsFMP.csv', na.string = "")
fiveFull <- fullset[fullset$StartDepth == 5,]
fiveFull <- fiveFull[fiveFull$EndDepth == 5,]
fiveFull <- fiveFull[!is.na(fiveFull$Param),]
fiveNN <- fiveFull[fiveFull$Param %in% c("NO3/NO2"),]
fiveNN <- fiveNN[, c("CollectDate", "Param", "ReportedResult",
                     "CorrectedReportedResult")]

fiveNN$Date <- sapply(strsplit(as.character(fiveNN$CollectDate), split = " "), 
                      "[", 1)
fiveNN$Date <- as.Date(fiveNN$Date, format = "%m/%d/%Y")
fiveNN$CorrectedReportedResult <- as.numeric(as.character(fiveNN$CorrectedReportedResult))
subFiveNN <- fiveNN[fiveNN$Date %in% unique(mldData$Date),]

trimmedFiveNN <- subFiveNN[, c("Date", "CorrectedReportedResult", "Param")]
miniFiveNN <- aggregate(.~Date+Param, trimmedFiveNN, FUN = mean)

summedFiveNN <- aggregate(CorrectedReportedResult~Date, miniFiveNN, FUN = sum)
summedFiveNN$CorrectedReportedResult <- summedFiveNN$CorrectedReportedResult/14

summedFiveNN$Month <- format(summedFiveNN$Date, "%m")
NNMonthly5 <- do.call(data.frame, 
                       aggregate(CorrectedReportedResult~Month, summedFiveNN, 
                                 function(x) c(mean = mean(x), sd = sd(x))))


ninetyFull <- fullset[fullset$StartDepth == 90,]
ninetyFull <- ninetyFull[ninetyFull$EndDepth == 90,]
ninetyFull <- ninetyFull[!is.na(ninetyFull$Param),]
ninetyNN <- ninetyFull[ninetyFull$Param %in% c("NO3/NO2"),]
ninetyNN <- ninetyNN[, c("CollectDate", "Param", "ReportedResult", "CorrectedReportedResult")]

ninetyNN$Date <- sapply(strsplit(as.character(ninetyNN$CollectDate), split = " "), 
                        "[", 1)
ninetyNN$Date <- as.Date(ninetyNN$Date, format = "%m/%d/%Y")
ninetyNN$CorrectedReportedResult <- as.numeric(as.character(ninetyNN$CorrectedReportedResult))
subNinetyNN <- ninetyNN[ninetyNN$Date %in% unique(mldData$Date),]

trimmedNinetyNN <- subNinetyNN[, c("Date", "CorrectedReportedResult", "Param")]
miniNinetyNN <- aggregate(.~Date+Param, trimmedNinetyNN, FUN = mean)

summedNinetyNN <- aggregate(CorrectedReportedResult~Date, miniNinetyNN, FUN = sum)
summedNinetyNN$CorrectedReportedResult <- summedNinetyNN$CorrectedReportedResult/14
summedNinetyNN$Month <- format(summedNinetyNN$Date, "%m")
NNMonthly90 <- do.call(data.frame, 
                     aggregate(CorrectedReportedResult~Month, summedNinetyNN, 
                               function(x) c(mean = mean(x), sd = sd(x))))




##SRP
fullset <- read.csv('~/FlatheadPublic/FLBSPublicData_AllNuts_Jan2023.csv', 
                    na.string = "")
fullset <- fullset[fullset$Site == "Flathead Lake, Midlake Deep",]
fiveFull <- fullset[fullset$Start.Depth == 5,]
fiveFull <- fiveFull[fiveFull$End.Depth == 5,]
fiveFull <- fiveFull[!is.na(fiveFull$Parameter),]

fiveSRP <- fiveFull[fiveFull$Parameter %in% c("SRP"),]
fiveSRP <- fiveSRP[, c("Date", "Parameter", "Value")]


fiveSRP$Date <- as.Date(fiveSRP$Date, format = "%m/%d/%Y")
fiveSRP$Value <- as.numeric(as.character(fiveSRP$Value))
fiveSRP[is.na(fiveSRP$Value), "Value"] <- 0.8
subFiveSRP <- fiveSRP[fiveSRP$Date %in% unique(mldData$Date),]

miniFiveSRP <- aggregate(.~Date+Parameter, subFiveSRP, FUN = mean)

miniFiveSRP$Value <- miniFiveSRP$Value/31

miniFiveSRP$Month <- format(miniFiveSRP$Date, "%m")
SRPMonthly5 <- do.call(data.frame, 
                      aggregate(Value~Month, miniFiveSRP, 
                                function(x) c(mean = mean(x), sd = sd(x))))


ninetyFull <- fullset[fullset$Start.Depth == 90,]
ninetyFull <- ninetyFull[ninetyFull$End.Depth == 90,]
ninetyFull <- ninetyFull[!is.na(ninetyFull$Parameter),]
ninetySRP <- ninetyFull[ninetyFull$Parameter %in% c("SRP"),]
ninetySRP <- ninetySRP[, c("Date", "Parameter", "Value")]

ninetySRP$Date <- as.Date(ninetySRP$Date, format = "%m/%d/%Y")
ninetySRP$Value <- as.numeric(as.character(ninetySRP$Value))
ninetySRP[is.na(ninetySRP$Value), "Value"] <- 0.8
subNinetySRP <- ninetySRP[ninetySRP$Date %in% unique(mldData$Date),]

trimmedNinetySRP <- subNinetySRP[, c("Date", "Value", "Parameter")]
miniNinetySRP <- aggregate(.~Date+Parameter, trimmedNinetySRP, FUN = mean)

miniNinetySRP$Value <- miniNinetySRP$Value/31
miniNinetySRP$Month <- format(miniNinetySRP$Date, "%m")
SRPMonthly90 <- do.call(data.frame, 
                       aggregate(Value~Month, miniNinetySRP, 
                                 function(x) c(mean = mean(x), sd = sd(x))))
