smallAbund0 <- read.csv("SubsampledData/small_updatedNames.count_table", 
                        stringsAsFactors = FALSE, row.names = 1)
#Rename stuff
for (i in 1){
  colnames(smallAbund0) <- gsub("Sept", "Sep", colnames(smallAbund0))
  colnames(smallAbund0) <- gsub("June", "Jun", colnames(smallAbund0))
  colnames(smallAbund0) <- gsub("July", "Jul", colnames(smallAbund0))
  colnames(smallAbund0) <- sub("[^_]*_(.*)", "\\1", colnames(smallAbund0))
}

smallAbund <- smallAbund0
smallAbund$RelAbund <- rowSums(smallAbund0[,colnames(smallAbund0) != "Sequence"])/
  ((ncol(smallAbund0)-1)*100)
smallAbund$MinVal <- apply(smallAbund0[,colnames(smallAbund0) != "Sequence"]/100, 
                           1, min)
smallAbund$MaxVal <- apply(smallAbund0[,colnames(smallAbund0) != "Sequence"]/100, 
                           1, max)
smallAbund$Sd <- apply(smallAbund0[,colnames(smallAbund0) != "Sequence"]/100, 1, sd)

reorderedSmall <- smallAbund[order(smallAbund$RelAbund, decreasing = TRUE),]
easyReorderedSmall <- smallAbund0[order(smallAbund$RelAbund, decreasing = TRUE),]


top500Small <- reorderedSmall[seq(1,250),c("Sequence", "RelAbund", "MinVal", 
                                           "MaxVal", "Sd")]

top500Small[, colnames(top500Small) != "Sequence"] <- 
  round(top500Small[, colnames(top500Small) != "Sequence"], 2)

twoPartsS <- 
  strsplit(colnames(easyReorderedSmall[, 
                                       colnames(easyReorderedSmall) != "Sequence"]),
           "_")
subDatesS <- sapply(twoPartsS, `[[`, 1)
subDatesS <- as.Date(subDatesS, format = "%d%b%y")
subDepthsS <- sapply(twoPartsS, `[[`, 2)
subDepthsS[length(subDepthsS)+1] <- "Null"

library("lomb")
lspOTUS5 <- apply(easyReorderedSmall[seq(1,250),subDepthsS == "05m"], 1, 
                  function(x) lsp(unname(unlist(x)), 
                                  times = as.numeric(subDatesS)[subDepthsS == "05m"],
                                  type = "period", plot = FALSE))
lspOTUS10 <- apply(easyReorderedSmall[seq(1,250),subDepthsS == "10m"], 1, 
                   function(x) lsp(unname(unlist(x)), 
                                   times = 
                                     as.numeric(subDatesS)[subDepthsS == "10m"],
                                   type = "period", plot = FALSE))
lspOTUSDCM <- apply(easyReorderedSmall[seq(1,250),subDepthsS == "DCM"], 1, 
                    function(x) lsp(unname(unlist(x)), 
                                    times = 
                                      as.numeric(subDatesS)[subDepthsS == "DCM"],
                                    type = "period", plot = FALSE))
lspOTUS50 <- apply(easyReorderedSmall[seq(1,250),subDepthsS == "50m"], 1, 
                   function(x) lsp(unname(unlist(x)), 
                                   times = 
                                     as.numeric(subDatesS)[subDepthsS == "50m"],
                                   type = "period", plot = FALSE))
lspOTUS90 <- apply(easyReorderedSmall[seq(1,250),subDepthsS == "90m"], 1, 
                   function(x) lsp(unname(unlist(x)), 
                                   times = 
                                     as.numeric(subDatesS)[subDepthsS == "90m"],
                                   type = "period", plot = FALSE))

top500Small$Seasonality_5m <- unlist(lapply(lspOTUS5, `[`, 'p.value'))
top500Small$Period_5m <- sapply(lapply(lspOTUS5, `[`, 'peak.at'),`[[`, 1)[1,]

top500Small[top500Small$Period_5m > 400 | top500Small$Period_5m < 345, 
            "Seasonality_5m"] <- NA
top500Small[top500Small$Period_5m > 400 | top500Small$Period_5m < 345, 
            "Period_5m"] <- NA
top500Small[which(top500Small$Seasonality_5m > 0.05), "Period_5m"] <- NA
top500Small[which(top500Small$Seasonality_5m > 0.05), "Seasonality_5m"] <- NA
top500Small[which(top500Small$Seasonality_5m < 0.01), "Seasonality_5m"] <- "< 0.01"


top500Small$Seasonality_10m <- unlist(lapply(lspOTUS10, `[`, 'p.value'))
top500Small$Period_10m <- sapply(lapply(lspOTUS10, `[`, 'peak.at'),`[[`, 1)[1,]

top500Small[top500Small$Period_10m > 400 | top500Small$Period_10m < 345, 
            "Seasonality_10m"] <- NA
top500Small[top500Small$Period_10m > 400 | top500Small$Period_10m < 345, 
            "Period_10m"] <- NA
top500Small[which(top500Small$Seasonality_10m > 0.05), "Period_10m"] <- NA
top500Small[which(top500Small$Seasonality_10m > 0.05), "Seasonality_10m"] <- NA
top500Small[which(top500Small$Seasonality_10m < 0.01), "Seasonality_10m"] <- "< 0.01"

top500Small$Seasonality_DCM <- unlist(lapply(lspOTUSDCM, `[`, 'p.value'))
top500Small$Period_DCM <- sapply(lapply(lspOTUSDCM, `[`, 'peak.at'),`[[`, 1)[1,]

top500Small[top500Small$Period_DCM > 400 | top500Small$Period_DCM < 345, 
            "Seasonality_DCM"] <- NA
top500Small[top500Small$Period_DCM > 400 | top500Small$Period_DCM < 345, 
            "Period_DCM"] <- NA
top500Small[which(top500Small$Seasonality_DCM > 0.05), "Period_DCM"] <- NA
top500Small[which(top500Small$Seasonality_DCM > 0.05), "Seasonality_DCM"] <- NA
top500Small[which(top500Small$Seasonality_DCM < 0.01), "Seasonality_DCM"] <- "< 0.01"


top500Small$Seasonality_50m <- unlist(lapply(lspOTUS50, `[`, 'p.value'))
top500Small$Period_50m <- sapply(lapply(lspOTUS50, `[`, 'peak.at'),`[[`, 1)[1,]

top500Small[top500Small$Period_50m > 400 | top500Small$Period_50m < 345, 
            "Seasonality_50m"] <- NA
top500Small[top500Small$Period_50m > 400 | top500Small$Period_50m < 345, 
            "Period_50m"] <- NA
top500Small[which(top500Small$Seasonality_50m > 0.05), "Period_50m"] <- NA
top500Small[which(top500Small$Seasonality_50m > 0.05), "Seasonality_50m"] <- NA
top500Small[which(top500Small$Seasonality_50m < 0.01), "Seasonality_50m"] <- "< 0.01"



top500Small$Seasonality_90m <- unlist(lapply(lspOTUS90, `[`, 'p.value'))
top500Small$Period_90m <- sapply(lapply(lspOTUS90, `[`, 'peak.at'),`[[`, 1)[1,]

top500Small[top500Small$Period_90m > 400 | top500Small$Period_90m < 345, 
            "Seasonality_90m"] <- NA
top500Small[top500Small$Period_90m > 400 | top500Small$Period_90m < 345, 
            "Period_90m"] <- NA
top500Small[which(top500Small$Seasonality_90m > 0.05), "Period_90m"] <- NA
top500Small[which(top500Small$Seasonality_90m > 0.05), "Seasonality_90m"] <- NA
top500Small[which(top500Small$Seasonality_90m < 0.01), "Seasonality_90m"] <- "< 0.01"





taxTable <- read.csv("RawData/prok.trim.taxonomy")
allTop200 <- merge(top500Small, taxTable, by.x = "row.names", by.y = "seqID", 
                   all.x = TRUE)

allTop200 <- allTop200[order(as.numeric(gsub("Otu", "", allTop200$Sequence))),]
allTop200[, grep("Period", colnames(allTop200))] <- NULL
write.csv(allTop200, "seasonalityAndTax_250_Small_Jan2023.csv", quote = FALSE)






######Percent of phyla
taxTable <- read.csv("RawData/prok.trim.taxonomy")
smallAbund0 <- read.csv("SubsampledData/small_updatedNames.count_table", 
                        stringsAsFactors = FALSE, row.names = 1)
#Rename stuff
for (i in 1){
  colnames(smallAbund0) <- gsub("Sept", "Sep", colnames(smallAbund0))
  colnames(smallAbund0) <- gsub("June", "Jun", colnames(smallAbund0))
  colnames(smallAbund0) <- gsub("July", "Jul", colnames(smallAbund0))
  colnames(smallAbund0) <- sub("[^_]*_(.*)", "\\1", colnames(smallAbund0))
}

smallAbund <- smallAbund0
smallAbund$RelAbund <- rowSums(smallAbund0[,colnames(smallAbund0) != "Sequence"])/
  ((ncol(smallAbund0)-1)*100)
smallAbund$MinVal <- apply(smallAbund0[,colnames(smallAbund0) != "Sequence"]/100, 
                           1, min)
smallAbund$MaxVal <- apply(smallAbund0[,colnames(smallAbund0) != "Sequence"]/100, 
                           1, max)
smallAbund$Sd <- apply(smallAbund0[,colnames(smallAbund0) != "Sequence"]/100, 1, sd)

reorderedSmall <- smallAbund[order(smallAbund$RelAbund, decreasing = TRUE),]
easyReorderedSmall <- 
  smallAbund[order(smallAbund$MaxVal, decreasing = TRUE),
              c("Sequence", "RelAbund", "MinVal", "MaxVal", "Sd")]
trimmedSmall <- easyReorderedSmall[easyReorderedSmall$MaxVal > 0.01,]
trimmedSmall$Row <- 1:nrow(trimmedSmall)

allTop200 <- merge(trimmedSmall, taxTable, by.x = "row.names", by.y = "seqID", 
                   all.x = TRUE)
allTop200 <- allTop200[order(allTop200$Row), ]

chloroOTUs <- reorderedSmall[allTop200$phylum == "Chloroflexi",]
actinoOTUs <- reorderedSmall[allTop200$phylum == "Actinobacteria",]
proteoOTUs <- reorderedSmall[allTop200$phylum == "Proteobacteria",]
planctOTUs <- reorderedSmall[allTop200$phylum == "Planctomycetes",]
bacterOTUs <- reorderedSmall[allTop200$phylum == "Bacteroidetes",]


cyanos <- taxTable[grep("Cyanobium", taxTable$clade),]
cyanoFull <- merge(easyReorderedSmall, cyanos, by.x = "row.names", by.y = "seqID", 
                   all.x = TRUE)





###Assign seasonality
library("lomb")
lspOTUS5 <- apply(chloroOTUs, 1, 
                  function(x) lsp(unname(unlist(x)), 
                                  times = as.numeric(subDatesS),
                                  type = "period"))
lspOTUS10 <- apply(actinoOTUs, 1, 
                   function(x) lsp(unname(unlist(x)), 
                                   times = 
                                     as.numeric(subDatesS),
                                   type = "period"))
lspOTUSDCM <- apply(proteoOTUs, 1, 
                    function(x) lsp(unname(unlist(x)), 
                                    times = 
                                      as.numeric(subDatesS),
                                    type = "period"))
lspOTUS50 <- apply(planctOTUs, 1, 
                   function(x) lsp(unname(unlist(x)), 
                                   times = 
                                     as.numeric(subDatesS),
                                   type = "period"))
lspOTUS90 <- apply(bacterOTUs, 1, 
                   function(x) lsp(unname(unlist(x)), 
                                   times = 
                                     as.numeric(subDatesS),
                                   type = "period"))

top500Small$Seasonality_5m <- unlist(lapply(lspOTUS5, `[`, 'p.value'))
top500Small$Period_5m <- sapply(lapply(lspOTUS5, `[`, 'peak.at'),`[[`, 1)[1,]

top500Small[top500Small$Period_5m > 400 | top500Small$Period_5m < 345, 
            "Seasonality_5m"] <- NA
top500Small[top500Small$Period_5m > 400 | top500Small$Period_5m < 345, 
            "Period_5m"] <- NA
top500Small[which(top500Small$Seasonality_5m > 0.05), "Period_5m"] <- NA
top500Small[which(top500Small$Seasonality_5m > 0.05), "Seasonality_5m"] <- NA
top500Small[which(top500Small$Seasonality_5m < 0.01), "Seasonality_5m"] <- "< 0.01"


top500Small$Seasonality_10m <- unlist(lapply(lspOTUS10, `[`, 'p.value'))
top500Small$Period_10m <- sapply(lapply(lspOTUS10, `[`, 'peak.at'),`[[`, 1)[1,]

top500Small[top500Small$Period_10m > 400 | top500Small$Period_10m < 345, 
            "Seasonality_10m"] <- NA
top500Small[top500Small$Period_10m > 400 | top500Small$Period_10m < 345, 
            "Period_10m"] <- NA
top500Small[which(top500Small$Seasonality_10m > 0.05), "Period_10m"] <- NA
top500Small[which(top500Small$Seasonality_10m > 0.05), "Seasonality_10m"] <- NA
top500Small[which(top500Small$Seasonality_10m < 0.01), "Seasonality_10m"] <- "< 0.01"

top500Small$Seasonality_DCM <- unlist(lapply(lspOTUSDCM, `[`, 'p.value'))
top500Small$Period_DCM <- sapply(lapply(lspOTUSDCM, `[`, 'peak.at'),`[[`, 1)[1,]

top500Small[top500Small$Period_DCM > 400 | top500Small$Period_DCM < 345, 
            "Seasonality_DCM"] <- NA
top500Small[top500Small$Period_DCM > 400 | top500Small$Period_DCM < 345, 
            "Period_DCM"] <- NA
top500Small[which(top500Small$Seasonality_DCM > 0.05), "Period_DCM"] <- NA
top500Small[which(top500Small$Seasonality_DCM > 0.05), "Seasonality_DCM"] <- NA
top500Small[which(top500Small$Seasonality_DCM < 0.01), "Seasonality_DCM"] <- "< 0.01"


top500Small$Seasonality_50m <- unlist(lapply(lspOTUS50, `[`, 'p.value'))
top500Small$Period_50m <- sapply(lapply(lspOTUS50, `[`, 'peak.at'),`[[`, 1)[1,]

top500Small[top500Small$Period_50m > 400 | top500Small$Period_50m < 345, 
            "Seasonality_50m"] <- NA
top500Small[top500Small$Period_50m > 400 | top500Small$Period_50m < 345, 
            "Period_50m"] <- NA
top500Small[which(top500Small$Seasonality_50m > 0.05), "Period_50m"] <- NA
top500Small[which(top500Small$Seasonality_50m > 0.05), "Seasonality_50m"] <- NA
top500Small[which(top500Small$Seasonality_50m < 0.01), "Seasonality_50m"] <- "< 0.01"



top500Small$Seasonality_90m <- unlist(lapply(lspOTUS90, `[`, 'p.value'))
top500Small$Period_90m <- sapply(lapply(lspOTUS90, `[`, 'peak.at'),`[[`, 1)[1,]

top500Small[top500Small$Period_90m > 400 | top500Small$Period_90m < 345, 
            "Seasonality_90m"] <- NA
top500Small[top500Small$Period_90m > 400 | top500Small$Period_90m < 345, 
            "Period_90m"] <- NA
top500Small[which(top500Small$Seasonality_90m > 0.05), "Period_90m"] <- NA
top500Small[which(top500Small$Seasonality_90m > 0.05), "Seasonality_90m"] <- NA
top500Small[which(top500Small$Seasonality_90m < 0.01), "Seasonality_90m"] <- "< 0.01"
