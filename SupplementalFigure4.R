library(vegan)

fullBG <- read.csv("~/FlatheadPublic/disgustingPooledData_DepthFriendly.csv", 
                   stringsAsFactors = FALSE, row.names = 1)

smallAbund <- 
  read.csv("~/FlatheadMicrobes/SubsampledData/small_updatedNames.count_table",
           stringsAsFactors = FALSE, row.names = 1)
smallAbund$Representative_Sequence <- NULL
trimmedSmall <- smallAbund[seq(1:150),]

countMat <- t(as.matrix(trimmedSmall))
rownames(countMat) <- gsub("Sept", "Sep", rownames(countMat))
rownames(countMat) <- gsub("June", "Jun", rownames(countMat))
rownames(countMat) <- gsub("July", "Jul", rownames(countMat))
rownames(countMat) <- sub("[^_]*_(.*)", "\\1", rownames(countMat))
allSamps <- rownames(countMat)

for (i in 1:length(allSamps)){
  targSamp <- allSamps[i]
  twoParts1 <- strsplit(targSamp, "_")
  subDates1 <- sapply(twoParts1, `[[`, 1)
  subDates1 <- as.Date(subDates1, format = "%d%b%y")
  subDepth1 <- sapply(twoParts1, `[[`, 2)
  if (subDepth1 != "DCM"){
    subDepth1 <- as.numeric(gsub("m", "", subDepth1))
    targRow <- fullBG[fullBG$Date == subDates1 & fullBG$Depth == subDepth1,]
    if (nrow(targRow) < 1){
      miniFram <- fullBG[fullBG$Date == subDates1,]
      targRow <- miniFram[which.min(abs(unique(miniFram$Depth - subDepth1))),]
    }
  } else {
    miniFram <- fullBG[fullBG$Date == subDates1,]
    dcmDepth <- mean(miniFram[,"DCM"])
    targRow <- miniFram[which.min(abs(unique(miniFram$Depth - dcmDepth))),]
  }
  
  if(nrow(targRow) < 1){
    next
  }
  
  rownames(targRow) <- targSamp
  if (exists("enviroSmall")){
    enviroSmall <- rbind(enviroSmall, targRow)
  } else {
    enviroSmall <- targRow
  }
}

countMatTrim <- countMat[allSamps %in% rownames(enviroSmall),]
#Add column of randomly generated values
enviroSmall$Date <- NULL
#enviroSmall0 <- enviroSmall
enviroSmall$NO3 <- NULL
enviroSmall$DCM <- NULL
enviroSmall$COND <- NULL
enviroSmall$ORP <- NULL

trimmES <- na.omit(enviroSmall)
trimmCount <- countMatTrim[rownames(countMatTrim) %in% rownames(trimmES),]
library(tictoc)

selectES <- trimmES
selectES$NH3 <- NULL
selectES$pH <- NULL
selectES$D.O..1 <- NULL
selectES$TURB <- NULL
selectES$TOC <- NULL
selectES$TN <- NULL
selectES$TP <- NULL
selectES$DOC <- NULL
selectES$CHLa.1 <- NULL

colnames(selectES) <- c("Depth", "Temp", "Dissolved O2", "Chl a", "Surface temp", 
                        "Mixed layer depth", "N+N", "Incident PAR")

tic()
enviroControlsTrimmed <- bioenv(trimmCount ~ ., selectES, upto = 6)
toc()

#tic()
#enviroControlsFull <- bioenv(countMatTrim ~ ., enviroSmall0, upto = 6)
#toc()

#upto 2 = 210, 0.909 sec (0.004)
#upto 3 = 1350, 5.811 sec (0.004)
#upto 4 = 6195, 24.964 sec (0.004)
#upto 6 = 60459, 247.422 sec (0.004)
#upto 7 = 137979, 652.776 sec (0.005)
distMat <- bioenvdist(enviroControlsTrimmed, which = "best")


adonisResults <- 
  adonis2(trimmCount ~ ., data = selectES)
adonisResults <- 
  adonisResults[order(adonisResults$R2)[1:(nrow(adonisResults)-2)],]
pValuesAdonis <- adonisResults$`Pr(>F)`
pValuesDf <- data.frame("Params" = rownames(adonisResults), 
                        "PValues" = pValuesAdonis)[seq(1, 8),]
library(fuzzySim)
correctedP <- FDR(pvalues = pValuesDf)

y0 <- 0.06
yLength <- 0.9
y1 <- yLength+y0
xa0 <- 0.21
xLength <- 0.37
xDiff <- 0.04
xa1 <- xa0 + xLength
xb0 <- xa1 + xDiff
xb1 <- xb0 + xLength

png("FigBinExtras/adonisResults_trim.png", width = 1200, height = 1200)
plot.new()

par(new = "TRUE",plt = c(xa0, xa1, y0, y1),las = 1, cex.axis = 2.5)
barplot(adonisResults$R2, horiz = TRUE, 
        names.arg = gsub("`", "", rownames(adonisResults)), las = 1)

par(new = "TRUE",plt = c(xb0, xb1, y0, y1),las = 1, cex.axis = 2.5)
barplot(adonisResults$F, horiz = TRUE, las = 1)

dev.off()




esControl <- selectES
esControl$normalControl1 <- rnorm(nrow(esControl),100, 10)
esControl$normalControl2 <- rnorm(nrow(esControl),10, 10)
esControl$normalControl3 <- rnorm(nrow(esControl),1, 10)
esControl$betaControl1 <- rbeta(nrow(esControl),1,5)
adonisResults <- 
  adonis2(trimmCount ~ ., data = esControl)
adonisResults <- 
  adonisResults[order(adonisResults$R2)[1:(nrow(adonisResults)-2)],]


pVals <- data.frame("PVal" = adonisResults$`Pr(>F)`, 
                    "Names" = rownames(adonisResults),
                    "Sig" = NA)
pSigs <- data.frame("PVal" = c(0.001, 0.01, 0.05, 0.1), 
                    "Sigs" = c("***", "**", "*", "."))
for (i in 1:nrow(pVals)){
  if (pVals[i, "PVal"] <= pSigs[1, "PVal"]){
    pVals[i, "Sig"] <- pSigs[1, "Sigs"]
  } else if (pVals[i, "PVal"] < pSigs[2, "PVal"]){
    pVals[i, "Sig"] <- pSigs[2, "Sigs"]
  } else if (pVals[i, "PVal"] < pSigs[3, "PVal"]){
    pVals[i, "Sig"] <- pSigs[3, "Sigs"]
  } else if (pVals[i, "PVal"] < pSigs[4, "PVal"]){
    pVals[i, "Sig"] <- pSigs[4, "Sigs"]
  } else{
    pVals[i, "Sig"] <- ""
  }
}

y0 <- 0.06
ymax <- 0.98
xa0 <- 0.21
xb1 <- 0.95
xDiff <- 0.02
xLength <- (xb1 - xa0 - xDiff)/2
xa1 <- xa0 + xLength
xb0 <- xa1 + xDiff

rowLabels <- rownames(adonisResults)

png("FigBinExtras/adonisResultsControl_trim.png", width = 1200, height = 1200)
plot.new()

par(new = "TRUE",plt = c(xa0, xa1, y0, ymax),las = 1, cex.axis = 2.5)
barplot(adonisResults$R2, horiz = TRUE, names.arg = gsub("`", "", rowLabels), las = 1)
mtext("R2 Value", side = 3, line = -2, cex = 2.5, font = 2)

par(new = "TRUE",plt = c(xb0, xb1, y0, ymax),las = 1, cex.axis = 2.5)
p1 <- barplot(adonisResults$F, horiz = TRUE, las = 1)
mtext("F Value", side = 3, line = -2, cex = 2.5, font = 2)
text(x = 90, y = p1, pVals$Sig, xpd = NA, cex = 2.5, font = 2, adj = 0)

dev.off()

tiff("FigBinExtras/adonisResultsControl_trim.tiff", width = 7, height = 6, pointsize = 12, 
     units = "in", res = 300)

par(plt = c(xa0, xa1, y0, ymax),las = 1)
barplot(adonisResults$R2, horiz = TRUE, names.arg = gsub("`", "", rowLabels), las = 1)
mtext("R2 Value", side = 3, line = -0.5, font = 2)

par(new = "TRUE",plt = c(xb0, xb1, y0, ymax),las = 1)
p1 <- barplot(adonisResults$F, horiz = TRUE, las = 1)
mtext("F Value", side = 3, line = -0.5, font = 2)
text(x = 90, y = p1, pVals$Sig, xpd = NA, font = 2, adj = 0)

dev.off()

tic()
enviroControlsTrimmed <- bioenv(trimmCount ~ ., esControl, upto = 7,
                                parallel = 4, index = "jaccard")
toc()

#upto 6 = 60459, parallel = 4, 102.319 sec (0.001)
#upto 7 = 137979, parallel = 4, 244.129  sec (0.001)

##Best model - two parameters: temperature (of sample) and pH
##correlation 0.7252
##Interestingly, depth is not included in these. Nor is mixed layer depth

adonisResultsMargin <- 
  adonis2(trimmCount ~ ., data = esControl, by = "margin")
adonisResultsMargin <- adonisResultsMargin[order(
  adonisResultsMargin$R2)[1:(nrow(adonisResultsMargin)-2)],]

adonisResultsNull <- 
  adonis2(trimmCount ~ ., data = esControl, by = "NULL")
adonisResultsNull <- adonisResultsNull[order(
  adonisResultsNull$R2)[1:(nrow(adonisResultsNull)-2)],]


##ANCOVA

##Best model - two parameters: temperature (of sample) and pH
##correlation 0.7252
##Interestingly, depth is not included in these. Nor is mixed layer depth

data(varechem)
sol <- bioenv(wisconsin(varespec) ~ log(N) + P + K + Ca + pH + Al, varechem)
sol
summary(sol)

# adonis example in vegan, page 13
# library(vegan)
# data(dune)
# data(dune.env)
# adonis2(dune ~ Management*A1, data = dune.env)
# adonis2(dune ~ Management*A1, data = dune.env, by = NULL)
# dat <- expand.grid(rep=gl(2,1), NO3=factor(c(0,10)),field=gl(3,1) )
# dat
# Agropyron <- with(dat, as.numeric(field) + as.numeric(NO3)+2) +rnorm(12)/2
# Schizachyrium <- with(dat, as.numeric(field) - as.numeric(NO3)+2) +rnorm(12)/2
# total <- Agropyron + Schizachyrium
# dotplot(total ~ NO3, dat, jitter.x=TRUE, groups=field, type=c('p','a'), 
#         xlab="NO3", auto.key=list(columns=3, lines=TRUE) )
# Y <- data.frame(Agropyron, Schizachyrium)
# mod <- metaMDS(Y, trace = FALSE)
# plot(mod)
# with(dat, ordiellipse(mod, field, kind = "ehull", label = TRUE))
# with(dat, ordispider(mod, field, lty=3, col="red"))
# perm <- how(nperm = 199)
# adonis2(Y ~ NO3, data = dat, permutations = perm)
# setBlocks(perm) <- with(dat, field)
# adonis2(Y ~ NO3, data = dat, permutations = perm)
