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
allEuks <- NULL
allCyas <- NULL
for (i in seq(1, length(relTotals))){
  df <- relTotals[[i]]
  df$Depth <- gsub('m','', gsub(' ', '', gsub('I', '', 
                                              gsub('\\s*\\([^\\)]+\\)','', df$`Sample Name`))))
  date <- as.Date(unlist(lapply(strsplit(names(totals), " "), `[[`, 1)), format = "%m.%d.%y")[i]
  
  subdfEuk <- df[grep("Euks", df$`Gate Name`),]
  subdfEuk$Concentration <- as.numeric(subdfEuk$Concentration)
  subdfCya <- df[grep("Cyan", df$`Gate Name`),]
  subdfCya$Concentration <- as.numeric(subdfCya$Concentration)
  
  agEuk <- aggregate(Concentration~`Sample Name`+Depth, subdfEuk, FUN = sum)
  agEuk$Depth <- gsub("-", "", agEuk$Depth)
  
  meanEuk <- aggregate(Concentration~Depth, agEuk, FUN = mean)
  
  meanEuk$Date <- date
  colnames(meanEuk) <- c("Depth", "CellConc", "Date")
  
  if(!exists("allEuks") || is.null(allEuks)){
    allEuks <- meanEuk
  } else{
    allEuks <- rbind(allEuks, meanEuk)
  }
  
  agCya <- aggregate(Concentration~`Sample Name`+Depth, subdfCya, FUN = sum)
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


allCyas[tolower(allCyas$Depth) == "chlaax", "Depth"] <- "DCM"
allCyas[tolower(allCyas$Depth) == "chla", "Depth"] <- "DCM"

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

colnames(allCyas) <- c("Depth", "CyaConc", "Date")
colnames(allEuks) <- c("Depth", "EukConc", "Date")
colnames(allCounts) <- c("Depth", "TotalConc", "Date")

pigAlls <- merge(allCyas, allEuks, by = c("Date", "Depth"))
cellAlls <- merge(pigAlls, allCounts, by = c("Date", "Depth"))

goodDepthsCounts <- cellAlls[which(cellAlls$Depth %in% c(5, 10, "10DCM", "DCM", 50, 90)),]

allDepths <- rbind(goodDepthsCounts, goodDepthsCounts[goodDepthsCounts$Depth == "10DCM",])
allDepths[nrow(allDepths), "Depth"] <- "DCM"
allDepths[allDepths$Depth == "10DCM", "Depth"] <- 10

allDepths$TotalConc <- allDepths$TotalConc - (allDepths$CyaConc+allDepths$EukConc)
allDepths$DOY <- as.numeric(as.character(format(allDepths$Date, "%j")))

mixStart <- c(-10, 333)
mixEnd <- c(128, 380)

allDepths <- allDepths[order(allDepths$DOY),]
annCols <- data.frame("Year" = rev(unique(format(allDepths$Date, "%Y"))), 
                      "Colors" = c("#D55E00", "#009E73", "#0072B2"))
allDepths$Col <- 
  as.character(annCols[match(format(allDepths$Date, "%Y"), annCols$Year), "Colors"])

monthDays <- 
  as.numeric(format(seq(as.Date("2011-1-1"), as.Date("2011-12-1"), by = "month"), "%j"))
thirdMonths <- monthDays[seq(3, 12, by = 3)]
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

lx1 <- 0.14
lx2 <- 0.48
rx1 <- 0.65
rx2 <- 0.99

png("FigBinExtras/cellCountsDuo.png", width = 1200, height = 1200)
plot.new()

##Total
###5m###
par(new = "TRUE",plt = c(lx1,lx2,0.82,0.98),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(m5Depths$DOY, m5Depths$TotalConc,
     type = "l", pch = 18, yaxt = "n", xaxt = "n", xlab = "", ylab = "",
     cex.main = 2.5, cex = 2.5, lwd = 2.5, cex.lab = 2.5,
     main = "", ylim = c(0, 1200))
rect(mixStart[1], -100, mixEnd[1], 1500, col = "#B4B4B4")
rect(mixStart[2], -100, mixEnd[2], 1500, col = "#B4B4B4")
lines(m5Depths$DOY, m5Depths$TotalConc, cex = 2.5, lwd = 2.5)
points(m5Depths$DOY, m5Depths$TotalConc, pch = 18, cex = 2.5, lwd = 2.5, 
       col = m5Depths$Col)
box(lwd = 3)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, padj = 1,
     labels = NA, at = thirdMonths, font = 2)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1, 
     at = monthDays, labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1, 
     at = seq(0, 1200, length.out = 5), labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.05,
     labels = seq(0, 12, length.out = 3), at = seq(0, 1200, length.out = 3), 
     font = 2)
text(300, 1100, "5 m", cex = 2.5, font = 2, las = 0)

###10###
par(new = "TRUE",plt = c(lx1,lx2,0.64,0.8),las = 1, cex.axis = 2.5)
plot(m10Depths$DOY, m10Depths$TotalConc,
     type = "l", pch = 18, yaxt = "n", xaxt = "n", xlab = "", ylab = "",
     cex.main = 2.5, cex = 2.5, lwd = 2.5, cex.lab = 2.5,
     main = "", ylim = c(0, 1200))
rect(mixStart[1], -100, mixEnd[1], 1500, col = "#B4B4B4")
rect(mixStart[2], -100, mixEnd[2], 1500, col = "#B4B4B4")
lines(m10Depths$DOY, m10Depths$TotalConc, cex = 2.5, lwd = 2.5)
points(m10Depths$DOY, m10Depths$TotalConc, pch = 18, cex = 2.5, lwd = 2.5, 
       col = m10Depths$Col)
box(lwd = 3)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, padj = 1,
     labels = NA, at = thirdMonths, font = 2)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1, 
     at = monthDays, labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1, 
     at = seq(0, 1200, length.out = 5), labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.05,
     labels = seq(0, 12, length.out = 3), at = seq(0, 1200, length.out = 3), font = 2)
text(300, 1100, "10 m", cex = 2.5, font = 2, las = 0)


#Chl Max
par(new = "TRUE",plt = c(lx1,lx2,0.46,0.62),las = 1, cex.axis = 2.5)
plot(DCMDepths$DOY, DCMDepths$TotalConc,
     type = "l", pch = 18, yaxt = "n", xaxt = "n", xlab = "", ylab = "",
     cex.main = 2.5, cex = 2.5, lwd = 2.5, cex.lab = 2.5,
     main = "", ylim = c(0, 1200))
rect(mixStart[1], -100, mixEnd[1], 1500, col = "#B4B4B4")
rect(mixStart[2], -100, mixEnd[2], 1500, col = "#B4B4B4")
lines(DCMDepths$DOY, DCMDepths$TotalConc, cex = 2.5, lwd = 2.5)
points(DCMDepths$DOY, DCMDepths$TotalConc, pch = 18, cex = 2.5, lwd = 2.5, 
       col = DCMDepths$Col)
box(lwd = 3)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, padj = 1,
     labels = NA, at = thirdMonths, font = 2)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1, 
     at = monthDays, labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1, 
     at = seq(0, 1200, length.out = 5), labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.05,
     labels = seq(0, 12, length.out = 3), at = seq(0, 1200, length.out = 3), 
     font = 2)
text(275, 1100, "Chl max", cex = 2.5, font = 2, las = 0)
mtext("Non-pigmented Picoplankton Abundance", side = 2, line = 7, font = 2, 
      cex = 2.5, las = 0)
mtext(bquote(bold("10"^{"5"}~"Cells mL"^{"-1"})), 
      side = 2, cex = 2.5, font = 2, line = 4, las = 0)


##50 m
par(new = "TRUE",plt = c(lx1,lx2,0.28,0.44),las = 1, cex.axis = 2.5)
plot(m50Depths$DOY, m50Depths$TotalConc,
     type = "l", pch = 18, yaxt = "n", xaxt = "n", xlab = "", ylab = "",
     cex.main = 2.5, cex = 2.5, lwd = 2.5, cex.lab = 2.5,
     main = "", ylim = c(0, 1200))
rect(mixStart[1], -100, mixEnd[1], 1500, col = "#B4B4B4")
rect(mixStart[2], -100, mixEnd[2], 1500, col = "#B4B4B4")
lines(m50Depths$DOY, m50Depths$TotalConc, cex = 2.5, lwd = 2.5)
points(m50Depths$DOY, m50Depths$TotalConc, pch = 18, cex = 2.5, lwd = 2.5, 
       col = m50Depths$Col)
box(lwd = 3)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, padj = 1,
     labels = NA, at = thirdMonths, font = 2)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1, 
     at = monthDays, labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1, 
     at = seq(0, 1200, length.out = 5), labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.05,
     labels = seq(0, 12, length.out = 3), at = seq(0, 1200, length.out = 3), 
     font = 2)
text(300, 1100, "50 m", cex = 2.5, font = 2, las = 0)


#90
par(new = "TRUE",plt = c(lx1,lx2,0.1,0.26),las = 1, cex.axis = 2.5)
plot(m90Depths$DOY, m90Depths$TotalConc,
     type = "l", pch = 18, yaxt = "n", xaxt = "n", xlab = "", ylab = "",
     cex.main = 2.5, cex = 2.5, lwd = 2.5, cex.lab = 2.5,
     main = "", ylim = c(0, 1200))
rect(mixStart[1], -100, mixEnd[1], 1500, col = "#B4B4B4")
rect(mixStart[2], -100, mixEnd[2], 1500, col = "#B4B4B4")
lines(m90Depths$DOY, m90Depths$TotalConc, cex = 2.5, lwd = 2.5)
points(m90Depths$DOY, m90Depths$TotalConc, pch = 18, cex = 2.5, lwd = 2.5, 
       col = m90Depths$Col)
box(lwd = 3)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, padj = 1,
     labels = month.abb[seq(3, 12, by = 3)], at = thirdMonths, font = 2)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1, 
     at = monthDays, labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1, 
     at = seq(0, 1200, length.out = 5), labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.05,
     labels = seq(0, 12, length.out = 3), at = seq(0, 1200, length.out = 3), font = 2)
text(300, 1100, "90 m", cex = 2.5, font = 2, las = 0)



##Cya
###5m###
par(new = "TRUE",plt = c(rx1,rx2,0.82,0.98),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(m5Depths$DOY, m5Depths$CyaConc,
     type = "l", pch = 18, yaxt = "n", xaxt = "n", xlab = "", ylab = "",
     cex.main = 2.5, cex = 2.5, lwd = 2.5, cex.lab = 2.5,
     main = "", ylim = c(0, 70))
rect(mixStart[1], -100, mixEnd[1], 1500, col = "#B4B4B4")
rect(mixStart[2], -100, mixEnd[2], 1500, col = "#B4B4B4")
lines(m5Depths$DOY, m5Depths$CyaConc, cex = 2.5, lwd = 2.5)
points(m5Depths$DOY, m5Depths$CyaConc, pch = 18, cex = 2.5, lwd = 2.5, col = m5Depths$Col)
box(lwd = 3)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, padj = 1,
     labels = NA, at = thirdMonths, font = 2)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1, 
     at = monthDays, labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1, 
     at = seq(0, 70, length.out = 5), labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.05,
     labels = seq(0, 7, length.out = 3), at = seq(0, 70, length.out = 3), font = 2)
text(300, 65, "5 m", cex = 2.5, font = 2, las = 0)


###10###
par(new = "TRUE",plt = c(rx1,rx2,0.64,0.8),las = 1, cex.axis = 2.5)
plot(m10Depths$DOY, m10Depths$CyaConc,
     type = "l", pch = 18, yaxt = "n", xaxt = "n", xlab = "", ylab = "",
     cex.main = 2.5, cex = 2.5, lwd = 2.5, cex.lab = 2.5,
     main = "", ylim = c(0, 70))
rect(mixStart[1], -100, mixEnd[1], 1500, col = "#B4B4B4")
rect(mixStart[2], -100, mixEnd[2], 1500, col = "#B4B4B4")
lines(m10Depths$DOY, m10Depths$CyaConc, cex = 2.5, lwd = 2.5)
points(m10Depths$DOY, m10Depths$CyaConc, pch = 18, cex = 2.5, lwd = 2.5, col = m10Depths$Col)
box(lwd = 3)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, padj = 1,
     labels = NA, at = thirdMonths, font = 2)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1, 
     at = monthDays, labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1, 
     at = seq(0, 70, length.out = 5), labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.05,
     labels = seq(0, 7, length.out = 3), at = seq(0, 70, length.out = 3), font = 2)
text(300, 65, "10 m", cex = 2.5, font = 2, las = 0)


#Chl Max
par(new = "TRUE",plt = c(rx1,rx2,0.46,0.62),las = 1, cex.axis = 2.5)
plot(DCMDepths$DOY, DCMDepths$CyaConc,
     type = "l", pch = 18, yaxt = "n", xaxt = "n", xlab = "", ylab = "",
     cex.main = 2.5, cex = 2.5, lwd = 2.5, cex.lab = 2.5,
     main = "", ylim = c(0, 210))
rect(mixStart[1], -100, mixEnd[1], 1500, col = "#B4B4B4")
rect(mixStart[2], -100, mixEnd[2], 1500, col = "#B4B4B4")
lines(DCMDepths$DOY, DCMDepths$CyaConc, cex = 2.5, lwd = 2.5)
points(DCMDepths$DOY, DCMDepths$CyaConc, pch = 18, cex = 2.5, lwd = 2.5, 
       col = DCMDepths$Col)
box(lwd = 3)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, padj = 1,
     labels = NA, at = thirdMonths, font = 2)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1, 
     at = monthDays, labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1, 
     at = seq(0, 210, length.out = 5), labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.05,
     labels = seq(0, 21, length.out = 3), at = seq(0, 210, length.out = 3), font = 2)
mtext("Picoplanktonic Cyanobacteria Abundance", side = 2, font = 2, cex = 2.5, 
      line = 9, las = 0)
mtext(bquote(bold("10"^{"4"}~"Cells mL"^{"-1"})), 
      side = 2, cex = 2.5, font = 2, line = 5.5, las = 0)
text(275, 195, "Chl max", cex = 2.5, font = 2, las = 0)


##50 m
par(new = "TRUE",plt = c(rx1,rx2,0.28,0.44),las = 1, cex.axis = 2.5)
plot(m50Depths$DOY, m50Depths$CyaConc,
     type = "l", pch = 18, yaxt = "n", xaxt = "n", xlab = "", ylab = "",
     cex.main = 2.5, cex = 2.5, lwd = 2.5, cex.lab = 2.5,
     main = "", ylim = c(0, 70))
rect(mixStart[1], -100, mixEnd[1], 1500, col = "#B4B4B4")
rect(mixStart[2], -100, mixEnd[2], 1500, col = "#B4B4B4")
lines(m50Depths$DOY, m50Depths$CyaConc, cex = 2.5, lwd = 2.5)
points(m50Depths$DOY, m50Depths$CyaConc, pch = 18, cex = 2.5, lwd = 2.5, col = m50Depths$Col)
box(lwd = 3)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, padj = 1,
     labels = NA, at = thirdMonths, font = 2)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1, 
     at = monthDays, labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1, 
     at = seq(0, 70, length.out = 5), labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.05,
     labels = seq(0, 7, length.out = 3), at = seq(0, 70, length.out = 3), font = 2)
text(300, 65, "50 m", cex = 2.5, font = 2, las = 0)


#90
par(new = "TRUE",plt = c(rx1,rx2,0.1,0.26),las = 1, cex.axis = 2.5)
plot(m90Depths$DOY, m90Depths$CyaConc,
     type = "l", pch = 18, yaxt = "n", xaxt = "n", xlab = "", ylab = "",
     cex.main = 2.5, cex = 2.5, lwd = 2.5, cex.lab = 2.5,
     main = "", ylim = c(0, 70))
rect(mixStart[1], -100, mixEnd[1], 1500, col = "#B4B4B4")
rect(mixStart[2], -100, mixEnd[2], 1500, col = "#B4B4B4")
lines(m90Depths$DOY, m90Depths$CyaConc, cex = 2.5, lwd = 2.5)
points(m90Depths$DOY, m90Depths$CyaConc, pch = 18, cex = 2.5, lwd = 2.5, col = m90Depths$Col)
box(lwd = 3)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.05, padj = 1,
     labels = month.abb[seq(3, 12, by = 3)], at = thirdMonths, font = 2)
axis(1, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1, 
     at = monthDays, labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.03, padj = 1, 
     at = seq(0, 70, length.out = 5), labels = NA)
axis(2, cex.axis=2.5, lwd.ticks = 3, tck = -0.05,
     labels = seq(0, 7, length.out = 3), at = seq(0, 70, length.out = 3), font = 2)
text(300, 65, "90 m", cex = 2.5, font = 2, las = 0)



par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

legend("bottom", legend = annCols$Year, text.font = 2, bty = "n",
       fill = as.character(annCols$Colors), horiz = TRUE, cex = 2.5, box.lwd = 0)

dev.off()


mySheetsTrim <- totals[1:5]
rbind(mySheetsTrim)
df <- do.call("rbind", mySheetsTrim)
