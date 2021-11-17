library(plotly)
library(DESeq2)

countTableLarge <- read.csv("~/FlatheadMicrobes/SubsampledData/rrBothLarge.count_table", 
                            stringsAsFactors = FALSE, row.names = 1)
countTableLarge$Representative_Sequence <- NULL
countTableSmall <- read.csv("~/FlatheadMicrobes/SubsampledData/rrBothSmall.count_table", 
                            stringsAsFactors = FALSE, row.names = 1)
countTableSmall$Representative_Sequence <- NULL



#countTable <- merge(countTableLarge, countTableSmall, by = "row.names", all = TRUE)

for (i in 1){
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
  surfTemp <- tempDF[tempDF$Depth == 1,]
  
  DCMDat <- read.csv("~/FlatheadPublic/DCMData.csv", header = TRUE, stringsAsFactors = FALSE)
  DCMDat$Date <- as.Date(DCMDat$Date)
  DCMDat$Date.1 <- as.Date(DCMDat$Date.1)
  DCMDat$Date.1.1 <- as.Date(DCMDat$Date.1.1)
  DCMDat <- DCMDat[order(DCMDat$Date),]
  DCMDat <- aggregate(DCM~Date, DCMDat, FUN = mean)
}

colDataLarge <- data.frame("Depth" = rep(NA, ncol(countTableLarge)), "SurfTemp" = NA,
                           "Date" = as.Date(NA))
colDataSmall <- data.frame("Depth" = rep(NA, ncol(countTableSmall)), "SurfTemp" = NA, 
                           "Date" = as.Date(NA))

for (i in 1){
  colnames(countTableLarge) <- gsub("Sept", "Sep", colnames(countTableLarge))
  colnames(countTableLarge) <- gsub("June", "Jun", colnames(countTableLarge))
  colnames(countTableLarge) <- gsub("July", "Jul", colnames(countTableLarge))
  colnames(countTableSmall) <- gsub("Sept", "Sep", colnames(countTableSmall))
  colnames(countTableSmall) <- gsub("June", "Jun", colnames(countTableSmall))
  colnames(countTableSmall) <- gsub("July", "Jul", colnames(countTableSmall))
  
  twoPartsL <- strsplit(colnames(countTableLarge), "_")
  subDatesL <- sapply(twoPartsL, `[[`, 2)
  subDatesL <- as.Date(subDatesL, format = "%d%b%y")
  depthsL <- sapply(twoPartsL, `[[`, 3)
  
  twoPartsS <- strsplit(colnames(countTableSmall), "_")
  subDatesS <- sapply(twoPartsS, `[[`, 2)
  subDatesS <- as.Date(subDatesS, format = "%d%b%y")
  depthsS <- sapply(twoPartsS, `[[`, 3)
  
  colDataLarge$Depth <- depthsL
  colDataSmall$Depth <- depthsS
  
  colDataLarge$SurfTemp <- tempDF[match(subDatesL, tempDF$Date), "Temp"]
  colDataSmall$SurfTemp <- tempDF[match(subDatesS, tempDF$Date), "Temp"]
  
  colDataLarge$Date <- subDatesL
  colDataSmall$Date <- subDatesS
  
  colDataLarge[colDataLarge$Depth == "DCM", "Depth"] <- 
    DCMDat[match(colDataLarge[colDataLarge$Depth == "DCM", "Date"], DCMDat$Date), "DCM"]
  colDataSmall[colDataSmall$Depth == "DCM", "Depth"] <- 
    DCMDat[match(colDataSmall[colDataSmall$Depth == "DCM", "Date"], DCMDat$Date), "DCM"]
  
  colDataLarge$Depth <- gsub("m", "", colDataLarge$Depth)
  colDataLarge$Depth <- as.numeric(as.character(colDataLarge$Depth))
  colDataSmall$Depth <- gsub("m", "", colDataSmall$Depth)
  colDataSmall$Depth <- as.numeric(as.character(colDataSmall$Depth))
  
  removeLarge <- which(rowSums(is.na(colDataLarge)) > 0)
  colDataLarge <- colDataLarge[-removeLarge,]
  countTableLarge <- countTableLarge[, -removeLarge]
  removeSmall <- which(rowSums(is.na(colDataSmall)) > 0)
  colDataSmall <- colDataSmall[-removeSmall,]
  countTableSmall <- countTableSmall[, -removeSmall]
}

taxTable <- read.csv("RawData/prok.98.taxonomy")
taxTable[,"tclass"] <- gsub("\\s*\\([^\\)]+\\)","",as.character(taxTable$class))

for (i in unique(taxTable$tclass)){
  targs <- taxTable[taxTable$tclass == i, "seqID"]
  largeDF <- countTableLarge[rownames(countTableLarge) %in% targs,]
  smallDF <- countTableSmall[rownames(countTableSmall) %in% targs,]
  
  classVecLarge <- colSums(largeDF)
  classVecSmall <- colSums(smallDF)
  
  if(exists("largeClass")){
    largeClass <- rbind(largeClass, classVecLarge)
  } else{
    largeClass <- classVecLarge
  }
  
  if(exists("smallClass")){
    smallClass <- rbind(smallClass, classVecSmall)
  } else{
    smallClass <- classVecSmall
  }
}

rownames(largeClass) <- unique(taxTable$tclass)
rownames(smallClass) <- unique(taxTable$tclass)

###Adding euks
###
eukTax <- read.table("~/FlatheadMicrobes/RawData/euks.taxonomy", sep = "\t", 
                     header = FALSE, stringsAsFactors = FALSE)
processTax <- function(taxTable){
  taxDat <- data.frame(do.call("rbind", strsplit(as.character(taxTable$V2), ";", 
                                                 fixed = TRUE)))
  newTax <- cbind(taxTable$V1, taxDat)
  colnames(newTax) <- c("seqID", "domain", "supergroup", "phylum", "class", "subclass", 
                        "order", "suborder", "family", "genus", "species")
  return(newTax)
}
eukTax <- processTax(eukTax)

eukTax[,"tclass"] <- gsub("\\s*\\([^\\)]+\\)","",as.character(eukTax$class))

for (i in unique(eukTax$tclass)){
  targs <- eukTax[eukTax$tclass == i, "seqID"]
  largeDF <- countTableLarge[rownames(countTableLarge) %in% targs,]
  smallDF <- countTableSmall[rownames(countTableSmall) %in% targs,]
  
  classVecLarge <- colSums(largeDF)
  classVecSmall <- colSums(smallDF)
  
  if(exists("largeEuk")){
    largeEuk <- rbind(largeEuk, classVecLarge)
  } else{
    largeEuk <- classVecLarge
  }
  
  if(exists("smallEuk")){
    smallEuk <- rbind(smallEuk, classVecSmall)
  } else{
    smallEuk <- classVecSmall
  }
}

rownames(largeEuk) <- unique(eukTax$tclass)
rownames(smallEuk) <- unique(eukTax$tclass)


largeClass <- rbind(largeClass, largeEuk)
smallClass <- rbind(smallClass, smallEuk)

largeSums <- rowSums(largeClass)
largeClass <- largeClass[order(largeSums, decreasing = TRUE),]
smallSums <- rowSums(smallClass)
smallClass <- smallClass[order(smallSums, decreasing = TRUE),]

deseqDSLargeClassST <- DESeqDataSetFromMatrix(largeClass, colDataLarge, 
                                              ~ Depth + SurfTemp)
deseqDSSmallClassST <- DESeqDataSetFromMatrix(smallClass, colDataSmall, 
                                              ~ Depth + SurfTemp)

deLargeOutputClassST <- DESeq(deseqDSLargeClassST)
deSmallOutputClassST <- DESeq(deseqDSSmallClassST)

deseqDSLargeClassD <- DESeqDataSetFromMatrix(largeClass, colDataLarge, 
                                             ~ SurfTemp + Depth)
deseqDSSmallClassD <- DESeqDataSetFromMatrix(smallClass, colDataSmall, 
                                             ~ SurfTemp + Depth)

deLargeOutputClassD <- DESeq(deseqDSLargeClassD)
deSmallOutputClassD <- DESeq(deseqDSSmallClassD)



classLargeTop10Depth <- results(deLargeOutputClassD)
classLargeTop10ST <- results(deLargeOutputClassST)
classSmallTop10Depth <- results(deSmallOutputClassD)
classSmallTop10ST <- results(deSmallOutputClassST)





#####Ready to plot
plotReadyLarge <- data.frame("Organism" = rownames(classLargeTop10ST), 
                             "STDE" = classLargeTop10ST$log2FoldChange, 
                             "DepthDE" = classLargeTop10Depth$log2FoldChange,
                             "AvgAbundance" = classLargeTop10Depth$baseMean,
                             "OrgColor" = "#1A85FF", "Symbol" = "circle")
plotReadyLarge$OrgColor <- as.character(plotReadyLarge$OrgColor)
plotReadyLarge$Organism <- as.character(plotReadyLarge$Organism)
plotReadyLarge$Symbol <- as.character(plotReadyLarge$Symbol)
plotReadyLarge[plotReadyLarge$Organism %in% eukTax$tclass, "OrgColor"] = "#D41159"
plotReadyLarge[plotReadyLarge$AvgAbundance > 120, "Symbol"] <- "star"
plotReadySmall <- data.frame("Organism" = rownames(classSmallTop10ST), 
                             "STDE" = classSmallTop10ST$log2FoldChange, 
                             "DepthDE" = classSmallTop10Depth$log2FoldChange,
                             "AvgAbundance" = classSmallTop10Depth$baseMean,
                             "OrgColor" = "#1A85FF", "Symbol" = "circle")
plotReadySmall$OrgColor <- as.character(plotReadySmall$OrgColor)
plotReadySmall$Organism <- as.character(plotReadySmall$Organism)
plotReadySmall$Symbol <- as.character(plotReadySmall$Symbol)
plotReadySmall[plotReadySmall$Organism %in% eukTax$tclass, "OrgColor"] <- "#D41159"
plotReadySmall[plotReadySmall$AvgAbundance > 120, "Symbol"] <- "star"


fig <- plot_ly(type = 'scatter', mode = 'markers') 
smallFig <- fig %>%
  add_trace(
    x = plotReadySmall$STDE*20.53, 
    y = -plotReadySmall$DepthDE*85,
    text = 
      paste0(plotReadySmall$Organism, ", ", 
             signif(plotReadySmall$AvgAbundance/120,digits=3), "%"),
    hoverinfo = 'text',
    marker = list(color=plotReadySmall$OrgColor, symbol = plotReadySmall$Symbol),
    showlegend = F
  ) %>%
  layout(title = 'Small (0.2-3 μm)',
         xaxis = list(title = 'log2-Fold Change (Surface Temperature)'), 
         yaxis = list(title = 'log2-Fold Change (Depth)'), 
         legend = list(title=list(text='<b> Organism Type </b>')))

largeFig <- fig %>%
  add_trace(
    x = plotReadyLarge$STDE*20.53, 
    y = -plotReadyLarge$DepthDE*85,
    text = 
      paste0(plotReadyLarge$Organism, ", ", 
             signif(plotReadyLarge$AvgAbundance/120,digits=3), "%"),
    hoverinfo = 'text',
    marker = list(color=plotReadyLarge$OrgColor, symbol = plotReadyLarge$Symbol),
    showlegend = F
  ) %>%
  layout(title = 'Large (>3 μm)',
         xaxis = list(title = 'log2-Fold Change (Surface Temperature)'), 
         yaxis = list(title = 'log2-Fold Change (Depth)'), 
         legend = list(title=list(text='<b> Organism Type </b>')))

largeFig <- largeFig %>% add_annotations(
  x=0,
  y=2,
  xref = "x",
  yref = "y",
  text = "Center Anchor",
  xanchor = 'center',
  showarrow = F
)
fig <- subplot(smallFig, largeFig, nrows = 1)

largeFig
