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

#relTotals <- totals[dates <= endDate & dates >= startDate]
relTotals <- totals
allCounts <- NULL
for (i in seq(1, length(relTotals))){
  df <- relTotals[[i]]
  colnames(df)[grep("Sample", colnames(df))] <- "SampleID"
  colnames(df)[grep("Gate", colnames(df))] <- "GateName"
  df$SampleID <- gsub(" 14 2500", "", df$SampleID)
  df$Depth <- gsub('m','', gsub(' ', '', gsub('I', '', 
                                              gsub('\\s*\\([^\\)]+\\)','', df$SampleID))))
  date <- as.Date(unlist(lapply(strsplit(names(totals), " "), `[[`, 1)), format = "%m.%d.%y")[i]
  
  subdf <- df[grep("NA|Pop", df$GateName),]
  subdf$Concentration <- as.numeric(subdf$Concentration)
  
  agAll <- aggregate(Concentration~SampleID+Depth, subdf, FUN = sum)
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
allCounts[tolower(allCounts$Depth) == "20dcm", "Depth"] <- "DCM"
allCounts[tolower(allCounts$Depth) == "14dcm", "Depth"] <- "DCM"
allCounts[tolower(allCounts$Depth) == "14dc", "Depth"] <- "DCM"
allCounts[tolower(allCounts$Depth) == "12dcm", "Depth"] <- "DCM"
allCounts[tolower(allCounts$Depth) == "12dc", "Depth"] <- "DCM"
allCounts[tolower(allCounts$Depth) == "16mdc", "Depth"] <- "DCM"
allCounts[tolower(allCounts$Depth) == "16mdcm", "Depth"] <- "DCM"
allCounts[tolower(allCounts$Depth) == "15dcm", "Depth"] <- "DCM"

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
allCounts[allCounts$Date==as.Date("2018-12-10") & allCounts$Depth == 17, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2019-01-09") & allCounts$Depth == 14, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2019-04-15") & allCounts$Depth == 16, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2019-05-06") & allCounts$Depth == 12, "Depth"] <- "DCM"
allCounts[allCounts$Date==as.Date("2019-05-29") & allCounts$Depth == 14, "Depth"] <- "DCM"

goodDepthsCounts <- allCounts[which(allCounts$Depth %in% c(5, 10, "10DCM", "DCM", 50, 90)),]

#allDepths <- rbind(goodDepthsCounts, goodDepthsCounts[goodDepthsCounts$Depth == "10DCM",])
#allDepths[nrow(allDepths), "Depth"] <- "DCM"
#allDepths[allDepths$Depth == "10DCM", "Depth"] <- 10

allDepths <- aggregate(CellConc ~ Date+Depth, goodDepthsCounts, FUN = mean)

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

relTotals <- totals[dates <= endDate]
allEuks <- NULL
allCyas <- NULL
for (i in seq(1, length(relTotals))){
  df <- relTotals[[i]]
  colnames(df)[grep("Sample", colnames(df))] <- "SampleID"
  colnames(df)[grep("Gate", colnames(df))] <- "GateName"
  df$Depth <- gsub('m','', gsub(' ', '', gsub('I', '', 
                                              gsub('\\s*\\([^\\)]+\\)','', df$SampleID))))
  date <- as.Date(unlist(lapply(strsplit(names(totals), " "), `[[`, 1)), format = "%m.%d.%y")[i]
  
  print(date)
  
  subdfEuk <- df[grep("Euks", df$GateName),]
  subdfEuk$Concentration <- as.numeric(subdfEuk$Concentration)
  subdfCya <- df[grep("Cyan", df$GateName),]
  subdfCya$Concentration <- as.numeric(subdfCya$Concentration)
  
  agEuk <- aggregate(Concentration~SampleID+Depth, subdfEuk, FUN = sum)
  agEuk$Depth <- gsub("-", "", agEuk$Depth)
  
  meanEuk <- aggregate(Concentration~Depth, agEuk, FUN = mean)
  
  meanEuk$Date <- date
  colnames(meanEuk) <- c("Depth", "CellConc", "Date")
  
  if(!exists("allEuks") || is.null(allEuks)){
    allEuks <- meanEuk
  } else{
    allEuks <- rbind(allEuks, meanEuk)
  }
  
  agCya <- aggregate(Concentration~SampleID+Depth, subdfCya, FUN = sum)
  agCya$Depth <- gsub("-", "", agCya$Depth)
  
  meanCya <- aggregate(Concentration~Depth, agCya, FUN = mean)
  
  meanCya$Date <- date
  colnames(meanCya) <- c("Depth", "CellConc", "Date")
  
  if(!exists("allCyas") || is.null(allCyas)){
    allCyas <- meanCya
  } else{
    allCyas <- rbind(allCyas, meanCya)
  }
}
allEuks[tolower(allEuks$Depth) == "chlaax", "Depth"] <- "DCM"
allEuks[tolower(allEuks$Depth) == "chla", "Depth"] <- "DCM"
allEuks[tolower(allEuks$Depth) == "20dcm", "Depth"] <- "DCM"
allEuks[tolower(allEuks$Depth) == "14dcm", "Depth"] <- "DCM"
allEuks[tolower(allEuks$Depth) == "14dc", "Depth"] <- "DCM"
allEuks[tolower(allEuks$Depth) == "12dcm", "Depth"] <- "DCM"
allEuks[tolower(allEuks$Depth) == "12dc", "Depth"] <- "DCM"
allEuks[tolower(allEuks$Depth) == "16mdc", "Depth"] <- "DCM"
allEuks[tolower(allEuks$Depth) == "16mdcm", "Depth"] <- "DCM"
allEuks[tolower(allEuks$Depth) == "15dcm", "Depth"] <- "DCM"

allEuks[allEuks$Date==as.Date("2016-09-09") & allEuks$Depth == 18, "Depth"] <- "DCM"
allEuks[allEuks$Date==as.Date("2016-10-26") & allEuks$Depth == 17, "Depth"] <- "DCM"
allEuks[allEuks$Date==as.Date("2016-11-16") & allEuks$Depth == 14, "Depth"] <- "DCM"
allEuks[allEuks$Date==as.Date("2016-12-22") & allEuks$Depth == 17, "Depth"] <- "DCM"
allEuks[allEuks$Date==as.Date("2017-01-20") & allEuks$Depth == 15, "Depth"] <- "DCM"
allEuks[allEuks$Date==as.Date("2017-02-21") & allEuks$Depth == 11, "Depth"] <- "DCM"
allEuks[allEuks$Date==as.Date("2017-03-22") & allEuks$Depth == 15, "Depth"] <- "DCM"
allEuks[allEuks$Date==as.Date("2017-04-11") & allEuks$Depth == 14, "Depth"] <- "DCM"
allEuks[allEuks$Date==as.Date("2017-04-24") & allEuks$Depth == 14, "Depth"] <- "DCM"
allEuks[allEuks$Date==as.Date("2017-06-07") & allEuks$Depth == 14, "Depth"] <- "DCM"
allEuks[allEuks$Date==as.Date("2017-06-30") & allEuks$Depth == 13, "Depth"] <- "DCM"
allEuks[allEuks$Date==as.Date("2018-05-21") & allEuks$Depth == 14, "Depth"] <- "DCM"
allEuks[allEuks$Date==as.Date("2018-07-05") & allEuks$Depth == 14, "Depth"] <- "DCM"
allEuks[allEuks$Date==as.Date("2018-07-25") & allEuks$Depth == 14, "Depth"] <- "DCM"
allEuks[allEuks$Date==as.Date("2018-08-13") & allEuks$Depth == 20, "Depth"] <- "DCM"
allEuks[allEuks$Date==as.Date("2018-09-06") & allEuks$Depth == 26, "Depth"] <- "DCM"
allEuks[allEuks$Date==as.Date("2018-09-24") & allEuks$Depth == 20, "Depth"] <- "DCM"
allEuks[allEuks$Date==as.Date("2018-11-13") & allEuks$Depth == 15, "Depth"] <- "DCM"
allEuks[allEuks$Date==as.Date("2018-12-10") & allEuks$Depth == 17, "Depth"] <- "DCM"

allEuks <- aggregate(CellConc ~ Date+Depth, allEuks, FUN = mean)

allCyas[tolower(allCyas$Depth) == "chlaax", "Depth"] <- "DCM"
allCyas[tolower(allCyas$Depth) == "chla", "Depth"] <- "DCM"
allCyas[tolower(allCyas$Depth) == "20dcm", "Depth"] <- "DCM"
allCyas[tolower(allCyas$Depth) == "14dcm", "Depth"] <- "DCM"
allCyas[tolower(allCyas$Depth) == "14dc", "Depth"] <- "DCM"
allCyas[tolower(allCyas$Depth) == "12dcm", "Depth"] <- "DCM"
allCyas[tolower(allCyas$Depth) == "12dc", "Depth"] <- "DCM"
allCyas[tolower(allCyas$Depth) == "16mdc", "Depth"] <- "DCM"
allCyas[tolower(allCyas$Depth) == "16mdcm", "Depth"] <- "DCM"
allCyas[tolower(allCyas$Depth) == "15dcm", "Depth"] <- "DCM"

allCyas[allCyas$Date==as.Date("2016-09-09") & allCyas$Depth == 18, "Depth"] <- "DCM"
allCyas[allCyas$Date==as.Date("2016-10-26") & allCyas$Depth == 17, "Depth"] <- "DCM"
allCyas[allCyas$Date==as.Date("2016-11-16") & allCyas$Depth == 14, "Depth"] <- "DCM"
allCyas[allCyas$Date==as.Date("2016-12-22") & allCyas$Depth == 17, "Depth"] <- "DCM"
allCyas[allCyas$Date==as.Date("2017-01-20") & allCyas$Depth == 15, "Depth"] <- "DCM"
allCyas[allCyas$Date==as.Date("2017-02-21") & allCyas$Depth == 11, "Depth"] <- "DCM"
allCyas[allCyas$Date==as.Date("2017-03-22") & allCyas$Depth == 15, "Depth"] <- "DCM"
allCyas[allCyas$Date==as.Date("2017-04-11") & allCyas$Depth == 14, "Depth"] <- "DCM"
allCyas[allCyas$Date==as.Date("2017-04-24") & allCyas$Depth == 14, "Depth"] <- "DCM"
allCyas[allCyas$Date==as.Date("2017-06-07") & allCyas$Depth == 14, "Depth"] <- "DCM"
allCyas[allCyas$Date==as.Date("2017-06-30") & allCyas$Depth == 13, "Depth"] <- "DCM"
allCyas[allCyas$Date==as.Date("2018-05-21") & allCyas$Depth == 14, "Depth"] <- "DCM"
allCyas[allCyas$Date==as.Date("2018-07-05") & allCyas$Depth == 14, "Depth"] <- "DCM"
allCyas[allCyas$Date==as.Date("2018-07-25") & allCyas$Depth == 14, "Depth"] <- "DCM"
allCyas[allCyas$Date==as.Date("2018-08-13") & allCyas$Depth == 20, "Depth"] <- "DCM"
allCyas[allCyas$Date==as.Date("2018-09-06") & allCyas$Depth == 26, "Depth"] <- "DCM"
allCyas[allCyas$Date==as.Date("2018-09-24") & allCyas$Depth == 20, "Depth"] <- "DCM"
allCyas[allCyas$Date==as.Date("2018-11-13") & allCyas$Depth == 15, "Depth"] <- "DCM"
allCyas[allCyas$Date==as.Date("2018-12-10") & allCyas$Depth == 17, "Depth"] <- "DCM"

allCyas <- aggregate(CellConc ~ Date+Depth, allCyas, FUN = mean)

colnames(allCyas) <- c("Date", "Depth", "CyaConc")
colnames(allEuks) <- c("Date", "Depth", "EukConc")
colnames(allDepths) <- c("Date", "Depth", "TotalConc")

pigAlls <- merge(allCyas, allEuks, by = c("Date", "Depth"), all = TRUE)
cellAlls <- merge(pigAlls, allDepths, by = c("Date", "Depth"), all = TRUE)

goodDepthsCounts <- cellAlls[which(cellAlls$Depth %in% c(5, 10, "10DCM", "DCM", 50, 90)),]

allDepths <- rbind(goodDepthsCounts, goodDepthsCounts[goodDepthsCounts$Depth == "10DCM",])
allDepths[nrow(allDepths), "Depth"] <- "DCM"
allDepths[allDepths$Depth == "10DCM", "Depth"] <- 10

allDepths$TotalConc <- allDepths$TotalConc - (allDepths$CyaConc+allDepths$EukConc)
allDepths$DOY <- as.numeric(as.character(format(allDepths$Date, "%j")))

mixStart <- c(-10, 333)
mixEnd <- c(128, 380)

allDepths <- allDepths[order(allDepths$DOY),]

annCols <- data.frame("Year" = c("2016", "2017", "2018"), 
                      "Colors" = c("#D55E00", "#009E73", "#0072B2"))
allDepths$Col <- 
  as.character(annCols[match(format(allDepths$Date, "%Y"), annCols$Year), "Colors"])

monthDays <- 
  as.numeric(format(seq(as.Date("2011-1-1"), as.Date("2011-12-1"), by = "month"), "%j"))
thirdMonths <- monthDays[seq(3, 12, by = 3)]

allDepths <- allDepths[order(allDepths$Date),]

m5Depths <- allDepths[allDepths$Depth == "5",]
m10Depths <- allDepths[allDepths$Depth == "10",]
DCMDepths <- allDepths[allDepths$Depth == "DCM",]
m50Depths <- allDepths[allDepths$Depth == "50",]
m90Depths <- allDepths[allDepths$Depth == "90",]


epiCounts <- rbind(m5Depths, m10Depths)
hypoCounts <- rbind(m50Depths, m90Depths)
allCounts <- rbind(epiCounts, DCMDepths)
allCounts <- rbind(allCounts, hypoCounts)

#All
mean(hypoCounts$CyaConc)
sd(hypoCounts$CyaConc)
mean(epiCounts$CyaConc)
sd(epiCounts$CyaConc)

mean(hypoCounts$TotalConc)
sd(hypoCounts$TotalConc)
mean(epiCounts$TotalConc)
sd(epiCounts$TotalConc)

#Winter
mean(hypoCounts[hypoCounts$DOY < 125 | hypoCounts$DOY > 330, "CyaConc"])
sd(hypoCounts[hypoCounts$DOY < 125 | hypoCounts$DOY > 330, "CyaConc"])
mean(epiCounts[epiCounts$DOY < 125 | epiCounts$DOY > 330, "CyaConc"])
sd(epiCounts[epiCounts$DOY < 125 | epiCounts$DOY > 330, "CyaConc"])


mean(hypoCounts[hypoCounts$DOY < 125 | hypoCounts$DOY > 330, "TotalConc"])
sd(hypoCounts[hypoCounts$DOY < 125 | hypoCounts$DOY > 330, "TotalConc"])
mean(epiCounts[epiCounts$DOY < 125 | epiCounts$DOY > 330, "TotalConc"])
sd(epiCounts[epiCounts$DOY < 125 | epiCounts$DOY > 330, "TotalConc"])
mean(allCounts[allCounts$DOY < 125 | allCounts$DOY > 330, "TotalConc"])
sd(allCounts[allCounts$DOY < 125 | allCounts$DOY > 330, "TotalConc"])

#Mixed: <125 or >330

#Summer
mean(hypoCounts[hypoCounts$DOY > 125 & hypoCounts$DOY < 330, "CyaConc"])
sd(hypoCounts[hypoCounts$DOY > 125 & hypoCounts$DOY < 330, "CyaConc"])
mean(epiCounts[epiCounts$DOY > 125 & epiCounts$DOY < 330, "CyaConc"])
sd(epiCounts[epiCounts$DOY > 125 & epiCounts$DOY < 330, "CyaConc"])


mean(hypoCounts[hypoCounts$DOY > 125 & hypoCounts$DOY < 330, "TotalConc"])
sd(hypoCounts[hypoCounts$DOY > 125 & hypoCounts$DOY < 330, "TotalConc"])
mean(epiCounts[epiCounts$DOY > 125 & epiCounts$DOY < 330, "TotalConc"])
sd(epiCounts[epiCounts$DOY > 125 & epiCounts$DOY < 330, "TotalConc"])

lx1 <- 0.07
lx2 <- 0.46
rx1 <- 0.605
rx2 <- 0.995

earlyFirst <- as.Date("2016-12-01")
earlyLast <- as.Date("2017-05-01")
lateFirst <- as.Date("2017-12-01")
lateLast <- as.Date("2018-05-15")

allMonths <- seq.Date(as.Date("2016-03-01"), as.Date("2020-03-01"), by = "month")
allThirds <- allMonths[c(TRUE, FALSE, FALSE)]
allTabs <- allMonths[c(FALSE, TRUE, TRUE)]

xLabelDate5 <- as.Date("2016-10-10")
xLabelDateBig <- as.Date("2016-10-20")
xLabelChl <- as.Date("2016-12-01")
cexSmall <- 0.9

xlimDat <- as.Date(c("2016-09-01", "2018-12-01"))

letters <- c("5 m", "5 m", "10 m", "10 m", "Chl max", "Chl max", "50 m", 
             "50 m", "90 m", "90 m")

tiff("FigBinExtras/cellCountsDuo.tiff", width = 7, height = 6, 
    pointsize = 12, units = "in", res = 300)
plot.new()

##Total
###5m###
par(new = "TRUE",plt = c(lx1,lx2,0.82,0.98),las = 1, xpd = FALSE)
plot(m5Depths$Date, m5Depths$TotalConc, type = "l", pch = 18, yaxt = "n", 
     xaxt = "n", xlab = "", ylab = "", ylim = c(0, 1200), xlim = xlimDat)
rect(earlyFirst, -100, earlyLast, 1500, col = "lightgray", border = NA)
rect(lateFirst, -100, lateLast, 1500, col = "lightgray", border = NA)
lines(m5Depths$Date, m5Depths$TotalConc)
points(m5Depths$Date, m5Depths$TotalConc, pch = 18)
box(lwd = 1)
axis(1, tck = -0.05, padj = 1, labels = NA, at = allThirds)
axis(1, tck = -0.03, padj = 1, at = allTabs, labels = NA)
axis(2, tck = -0.03, padj = 1, at = seq(0, 1200, length.out = 5), labels = NA)
axis(2, tck = -0.05, labels = NA, at = seq(0, 1200, length.out = 3))
abline(v = seq.Date(as.Date("2016-01-01"), as.Date("2021-01-01"), by = "year"),
       lty = 2, col = "darkgrey")
axLab <- seq(0, 12, length.out = 3)
axAdj <- seq(4.1, -3.2, length.out = length(axLab))
mtext(axLab, side = 2, line = 0.5, padj = axAdj, cex = cexSmall)
text(xLabelDate5, 1100, letters[1], las = 0)

###10###
par(new = "TRUE",plt = c(lx1,lx2,0.64,0.8),las = 1)
plot(m10Depths$Date, m10Depths$TotalConc,
     type = "l", pch = 18, yaxt = "n", xaxt = "n", xlab = "", ylab = "",
     main = "", ylim = c(0, 1200), xlim = xlimDat)
rect(earlyFirst, -100, earlyLast, 1500, col = "lightgray", border = NA)
rect(lateFirst, -100, lateLast, 1500, col = "lightgray", border = NA)
lines(m10Depths$Date, m10Depths$TotalConc)
box(lwd = 1)
points(m10Depths$Date, m10Depths$TotalConc, pch = 18)
axis(1, tck = -0.05, padj = 1, labels = NA, at = allThirds)
axis(1, tck = -0.03, padj = 1, at = allTabs, labels = NA)
axis(2, tck = -0.03, padj = 1, at = seq(0, 1200, length.out = 5), labels = NA)
axis(2, tck = -0.05, labels = NA, at = seq(0, 1200, length.out = 3))
abline(v = seq.Date(as.Date("2016-01-01"), as.Date("2021-01-01"), by = "year"),
       lty = 2, col = "darkgrey")
axLab <- seq(0, 12, length.out = 3)
axAdj <- seq(4.1, -3.2, length.out = length(axLab))
mtext(axLab, side = 2, line = 0.5, padj = axAdj, cex = cexSmall)
text(xLabelDateBig, 1100, letters[3], las = 0)


#Chl Max
par(new = "TRUE",plt = c(lx1,lx2,0.46,0.62),las = 1)
plot(DCMDepths$Date, DCMDepths$TotalConc,
     type = "l", pch = 18, yaxt = "n", xaxt = "n", xlab = "", ylab = "",
     main = "", ylim = c(0, 1200), xlim = xlimDat)
rect(earlyFirst, -100, earlyLast, 1500, col = "lightgray", border = NA)
rect(lateFirst, -100, lateLast, 1500, col = "lightgray", border = NA)
lines(DCMDepths$Date, DCMDepths$TotalConc)
points(DCMDepths$Date, DCMDepths$TotalConc, pch = 18)
box(lwd = 1)
axis(1, tck = -0.05, padj = 1, labels = NA, at = allThirds)
axis(1, tck = -0.03, padj = 1, at = allTabs, labels = NA)
axis(2, tck = -0.03, padj = 1, at = seq(0, 1200, length.out = 5), labels = NA)
axis(2, tck = -0.05, labels = NA, at = seq(0, 1200, length.out = 3))
abline(v = seq.Date(as.Date("2016-01-01"), as.Date("2021-01-01"), by = "year"),
       lty = 2, col = "darkgrey")
axLab <- seq(0, 12, length.out = 3)
axAdj <- seq(4.1, -3.2, length.out = length(axLab))
mtext(axLab, side = 2, line = 0.5, padj = axAdj, cex = cexSmall)
text(xLabelChl, 1100, letters[5], las = 0)
mtext(bquote("Non-pigmented picoplankton  (10"^{"5"}~"Cells mL"^{"-1"}~")"), 
      side = 2, line = 1.2, las = 0)


##50 m
par(new = "TRUE",plt = c(lx1,lx2,0.28,0.44),las = 1)
plot(m50Depths$Date, m50Depths$TotalConc,
     type = "l", pch = 18, yaxt = "n", xaxt = "n", xlab = "", ylab = "",
     main = "", ylim = c(0, 1200), xlim = xlimDat)
rect(earlyFirst, -100, earlyLast, 1500, col = "lightgray", border = NA)
rect(lateFirst, -100, lateLast, 1500, col = "lightgray", border = NA)
lines(m50Depths$Date, m50Depths$TotalConc)
points(m50Depths$Date, m50Depths$TotalConc, pch = 18)
box(lwd = 1)
axis(1, tck = -0.05, padj = 1, labels = NA, at = allThirds)
axis(1, tck = -0.03, padj = 1, at = allTabs, labels = NA)
axis(2, tck = -0.03, padj = 1, at = seq(0, 1200, length.out = 5), labels = NA)
axis(2, tck = -0.05, labels = NA, at = seq(0, 1200, length.out = 3))
abline(v = seq.Date(as.Date("2016-01-01"), as.Date("2021-01-01"), by = "year"),
       lty = 2, col = "darkgrey")
axLab <- seq(0, 12, length.out = 3)
axAdj <- seq(4.1, -3.2, length.out = length(axLab))
mtext(axLab, side = 2, line = 0.5, padj = axAdj, cex = cexSmall)
text(xLabelDateBig, 1100, letters[7], las = 0)


#90
par(new = "TRUE",plt = c(lx1,lx2,0.1,0.26),las = 1)
plot(m90Depths$Date, m90Depths$TotalConc, type = "l", pch = 18, yaxt = "n", 
     xaxt = "n", xlab = "", ylab = "", main = "", ylim = c(0, 1200), xlim = xlimDat)
rect(earlyFirst, -100, earlyLast, 1500, col = "lightgray", border = NA)
rect(lateFirst, -100, lateLast, 1500, col = "lightgray", border = NA)
lines(m90Depths$Date, m90Depths$TotalConc)
points(m90Depths$Date, m90Depths$TotalConc, pch = 18)
box(lwd = 1)
axis(1, tck = -0.05, padj = -1.3, labels = format(allThirds, "%b"), at = allThirds,
     cex.axis = cexSmall)
axis(1, tck = -0.03, padj = 1, at = allTabs, labels = NA)
axis(2, tck = -0.03, padj = 1, at = seq(0, 1200, length.out = 5), labels = NA)
axis(2, tck = -0.05, labels = NA, at = seq(0, 1200, length.out = 3))
abline(v = seq.Date(as.Date("2016-01-01"), as.Date("2021-01-01"), by = "year"),
       lty = 2, col = "darkgrey")
axLab <- seq(0, 12, length.out = 3)
axAdj <- seq(4.1, -3.2, length.out = length(axLab))
mtext(axLab, side = 2, line = 0.5, padj = axAdj, cex = cexSmall)
text(xLabelDateBig, 1100, letters[9], las = 0)
text(as.Date("2016-10-01"), -600, "2016", xpd = NA, cex = cexSmall)
text(as.Date("2017-07-01"), -600, "2017", xpd = NA, cex = cexSmall)
text(as.Date("2018-07-01"), -600, "2018", xpd = NA, cex = cexSmall)
mtext("Date", 1, line = 2, adj = 1.24)



##Cya
###5m###
par(new = "TRUE",plt = c(rx1,rx2,0.82,0.98),las = 1, xpd = FALSE)
plot(m5Depths$Date, m5Depths$CyaConc,
     type = "l", pch = 18, yaxt = "n", xaxt = "n", xlab = "", ylab = "",
     main = "", ylim = c(0, 70), xlim = xlimDat)
rect(earlyFirst, -100, earlyLast, 1500, col = "lightgray", border = NA)
rect(lateFirst, -100, lateLast, 1500, col = "lightgray", border = NA)
lines(m5Depths$Date, m5Depths$CyaConc)
points(m5Depths$Date, m5Depths$CyaConc, pch = 18)
box(lwd = 1)
axis(1, tck = -0.05, padj = 1, labels = NA, at = allThirds)
axis(1, tck = -0.03, padj = 1, at = allTabs, labels = NA)
axis(2, tck = -0.03, padj = 1, at = seq(0, 70, length.out = 5), labels = NA)
axis(2,tck = -0.05, labels = NA, at = seq(0, 70, length.out = 3))
abline(v = seq.Date(as.Date("2016-01-01"), as.Date("2021-01-01"), by = "year"),
       lty = 2, col = "darkgrey")
axLab <- seq(0, 7, length.out = 3)
axAdj <- seq(4.1, -3.2, length.out = length(axLab))
mtext(axLab, side = 2, line = 0.5, padj = axAdj, cex = cexSmall)
text(xLabelDate5, 65, letters[2], las = 0)


###10###
par(new = "TRUE",plt = c(rx1,rx2,0.64,0.8),las = 1)
plot(m10Depths$Date, m10Depths$CyaConc,
     type = "l", pch = 18, yaxt = "n", xaxt = "n", xlab = "", ylab = "",
     main = "", ylim = c(0, 70), xlim = xlimDat)
rect(earlyFirst, -100, earlyLast, 1500, col = "lightgray", border = NA)
rect(lateFirst, -100, lateLast, 1500, col = "lightgray", border = NA)
lines(m10Depths$Date, m10Depths$CyaConc)
points(m10Depths$Date, m10Depths$CyaConc, pch = 18)
box(lwd = 1)
axis(1, tck = -0.05, padj = 1, labels = NA, at = allThirds)
axis(1, tck = -0.03, padj = 1, at = allTabs, labels = NA)
axis(2, tck = -0.03, padj = 1, at = seq(0, 70, length.out = 5), labels = NA)
axis(2,tck = -0.05, labels = NA, at = seq(0, 70, length.out = 3))
abline(v = seq.Date(as.Date("2016-01-01"), as.Date("2021-01-01"), by = "year"),
       lty = 2, col = "darkgrey")
axLab <- seq(0, 7, length.out = 3)
axAdj <- seq(4.1, -3.2, length.out = length(axLab))
mtext(axLab, side = 2, line = 0.5, padj = axAdj, cex = cexSmall)
text(xLabelDateBig, 65, letters[4], las = 0)


#Chl Max
par(new = "TRUE",plt = c(rx1,rx2,0.46,0.62),las = 1)
plot(DCMDepths$Date, DCMDepths$CyaConc,
     type = "l", pch = 18, yaxt = "n", xaxt = "n", xlab = "", ylab = "",
     main = "", ylim = c(0, 210), xlim = xlimDat)
rect(earlyFirst, -100, earlyLast, 1500, col = "lightgray", border = NA)
rect(lateFirst, -100, lateLast, 1500, col = "lightgray", border = NA)
lines(DCMDepths$Date, DCMDepths$CyaConc)
points(DCMDepths$Date, DCMDepths$CyaConc, pch = 18)
box(lwd = 1)
axis(1, tck = -0.05, padj = 1, labels = NA, at = allThirds)
axis(1, tck = -0.03, padj = 1, at = allTabs, labels = NA)
axis(2, tck = -0.03, padj = 1, at = seq(0, 210, length.out = 5), labels = NA)
axis(2, tck = -0.05, labels = NA, at = seq(0, 210, length.out = 3))
abline(v = seq.Date(as.Date("2016-01-01"), as.Date("2021-01-01"), by = "year"),
       lty = 2, col = "darkgrey")
mtext(bquote("Picoplanktonic cyanobacteria (10"^{"4"}~"Cells mL"^{"-1"}~")"), 
      side = 2, line = 2, las = 0)
axLab <- seq(0, 21, length.out = 3)
axAdj <- seq(4.1, -3.2, length.out = length(axLab))
mtext(axLab, side = 2, line = 0.5, padj = axAdj, cex = cexSmall)
text(xLabelChl, 195, letters[6], las = 0)


##50 m
par(new = "TRUE",plt = c(rx1,rx2,0.28,0.44),las = 1)
plot(m50Depths$Date, m50Depths$CyaConc,
     type = "l", pch = 18, yaxt = "n", xaxt = "n", xlab = "", ylab = "",
     main = "", ylim = c(0, 70), xlim = xlimDat)
rect(earlyFirst, -100, earlyLast, 1500, col = "lightgray", border = NA)
rect(lateFirst, -100, lateLast, 1500, col = "lightgray", border = NA)
lines(m50Depths$Date, m50Depths$CyaConc)
points(m50Depths$Date, m50Depths$CyaConc, pch = 18)
box(lwd = 1)
axis(1, tck = -0.05, padj = 1, labels = NA, at = allThirds)
axis(1, tck = -0.03, padj = 1, at = allTabs, labels = NA)
axis(2, tck = -0.03, padj = 1, at = seq(0, 70, length.out = 5), labels = NA)
axis(2, tck = -0.05, labels = NA, at = seq(0, 70, length.out = 3))
abline(v = seq.Date(as.Date("2016-01-01"), as.Date("2021-01-01"), by = "year"),
       lty = 2, col = "darkgrey")
axLab <- seq(0, 7, length.out = 3)
axAdj <- seq(4.1, -3.2, length.out = length(axLab))
mtext(axLab, side = 2, line = 0.5, padj = axAdj, cex = cexSmall)
text(xLabelDateBig, 65, letters[8], las = 0)


#90
par(new = "TRUE",plt = c(rx1,rx2,0.1,0.26),las = 1)
plot(m90Depths$Date, m90Depths$CyaConc,
     type = "l", pch = 18, yaxt = "n", xaxt = "n", xlab = "", ylab = "",
     main = "", ylim = c(0, 70), xlim = xlimDat)
rect(earlyFirst, -100, earlyLast, 1500, col = "lightgray", border = NA)
rect(lateFirst, -100, lateLast, 1500, col = "lightgray", border = NA)
lines(m90Depths$Date, m90Depths$CyaConc)
points(m90Depths$Date, m90Depths$CyaConc, pch = 18)
box(lwd = 1)
axis(1, tck = -0.05, padj = -1.3, cex.axis = cexSmall,
     labels = format(allThirds, "%b"), at = allThirds)
axis(1, tck = -0.03, padj = 1, at = allTabs, labels = NA)
axis(2, tck = -0.03, padj = 1, at = seq(0, 70, length.out = 5), labels = NA)
axis(2, tck = -0.05, hadj = 0.5, labels = NA, at = seq(0, 70, length.out = 3))
abline(v = seq.Date(as.Date("2016-01-01"), as.Date("2021-01-01"), by = "year"),
       lty = 2, col = "darkgrey")
axLab <- seq(0, 7, length.out = 3)
axAdj <- seq(4.1, -3.2, length.out = length(axLab))
mtext(axLab, side = 2, line = 0.5, padj = axAdj, cex = cexSmall)
text(xLabelDateBig, 65, letters[10], las = 0)
text(as.Date("2016-10-01"), -35, "2016", xpd = NA, cex = cexSmall)
text(as.Date("2017-07-01"), -35, "2017", xpd = NA, cex = cexSmall)
text(as.Date("2018-07-01"), -35, "2018", xpd = NA, cex = cexSmall)

dev.off()

churchDCM <- read.csv("~/FlatheadPublic/ChurchLab/ChurchDCM.csv", stringsAsFactors = FALSE)
churchDCM$Date <- as.Date(churchDCM$Date, "%m-%d-%Y")
allDepths$Col <- NULL
allDepths$DOY <- NULL
allDepths$EukConc <- NULL

trimCells <- merge(allDepths, churchDCM, all.x = TRUE, by = "Date")
trimCells <- trimCells[!is.na(trimCells$CyaConc),]
trimCells$isDCM <- FALSE
trimCells[trimCells$Depth.x == "DCM", "isDCM"] <- TRUE
trimCells[trimCells$isDCM == "TRUE", "Depth.x"] <- trimCells[trimCells$isDCM == "TRUE", "Depth.y"]
trimCells$Depth.y <- NULL
colnames(trimCells) <- c("Date", "Depth", "Cyano", "NonPigmented", "isDCM")
write.csv(trimCells, "~/FlatheadMicrobes/cellCountDat.csv", quote = FALSE, row.names = FALSE)
