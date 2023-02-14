library(coop)

countTableSmall <- 
  read.csv("~/FlatheadMicrobes/SubsampledData/small_updatedNames.count_table", 
           stringsAsFactors = FALSE, row.names = 1)

for (i in 1){
  countTableSmall$Representative_Sequence <- NULL
  
  colnames(countTableSmall) <- sub("[^_]*_(.*)", "\\1", colnames(countTableSmall))
  colnames(countTableSmall) <- gsub("_[^_]+$", "\\1", colnames(countTableSmall))
  colnames(countTableSmall) <- gsub("Sept", "Sep", colnames(countTableSmall))
  colnames(countTableSmall) <- gsub("June", "Jun", colnames(countTableSmall))
}
rm(i)

trimmedSmall <- countTableSmall[rowSums(countTableSmall)>1,]

taxTable[,"tclass"] <- gsub("\\s*\\([^\\)]+\\)","",as.character(taxTable$class))

smallClass5m0 <- trimmedSmall[, grep('05m', colnames(trimmedSmall))]
smallClass5m <- smallClass5m0[,sample(1:ncol(smallClass5m0),30)]
smallClass10m0 <- trimmedSmall[, grep('10m', colnames(trimmedSmall))]
smallClass10m <- smallClass10m0[,sample(1:ncol(smallClass10m0),30)]
smallClassDCM <- trimmedSmall[, grep('DCM', colnames(trimmedSmall))]
smallClass50m0 <- trimmedSmall[, grep('50m', colnames(trimmedSmall))]
smallClass50m <- smallClass50m0[,sample(1:ncol(smallClass50m0),30)]
smallClass90m0 <- trimmedSmall[, grep('90m', colnames(trimmedSmall))]
smallClass90m <- smallClass90m0[,sample(1:ncol(smallClass90m0),30)]

#Based largely on equation A2
lehrmenStab2000 <- function(df){
  bioVal <- mean(colSums(df))
  miniMatt <- t(df)
  
  covarMatt <- sum(covar(miniMatt))
  sumCV <- sum(covarMatt)
  
  tempStab <- bioVal/(sumCV**(1/2))
  
  return(tempStab)
}

library(vegan)
stability5m <- lehrmenStab2000(smallClass5m)
shann5m <- mean(diversity(smallClass5m, "shannon"))
stability10m <- lehrmenStab2000(smallClass10m)
shann10m <-  mean(diversity(smallClass10m, "shannon"))
stabilityDCM <- lehrmenStab2000(smallClassDCM)
shannDCM <-  mean(diversity(smallClassDCM, "shannon"))
stability50m <- lehrmenStab2000(smallClass50m)
shann50m <-  mean(diversity(smallClass50m, "shannon"))
stability90m <- lehrmenStab2000(smallClass90m)
shann90m <-  mean(diversity(smallClass90m, "shannon"))
