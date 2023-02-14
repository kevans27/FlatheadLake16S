taxTable <- read.csv("RawData/prok.98.taxonomy")
countTableSmall <- read.csv("SubsampledData/small_updatedNames.count_table",
                            stringsAsFactors = FALSE, 
                            row.names = 1)
transKey <- read.csv("RawData/prok.name.key")

#Rename/trim things
for (i in 1){
  countTableSmall$Representative_Sequence <- NULL
  
  countTableSmall$blastname <- rownames(countTableSmall)
  colnames(countTableSmall) <- sub("[^_]*_(.*)", "\\1", colnames(countTableSmall))
  colnames(countTableSmall) <- gsub("_[^_]+$", "\\1", colnames(countTableSmall))
  colnames(countTableSmall) <- gsub("Sept", "Sep", colnames(countTableSmall))
  colnames(countTableSmall) <- gsub("June", "Jun", colnames(countTableSmall))
}

`%notin%` <- Negate(`%in%`)
relabundTax <- function(counts, taxes, level){
  subTax <- taxes[,c("seqID", level)]
  subTax[,"replaced"] <- gsub("\\s*\\([^\\)]+\\)","",as.character(subTax[,level]))
  counts$level <- subTax[match(counts$blastname, subTax$seqID),
                         "replaced"]
  
  #Get rid of the sequence names when aggregating, then sum up every sample's thing
  targs <- 
    colnames(counts)[colnames(counts) %notin% c("blastname", "seqID", "level")]
  levelCounts <- aggregate(counts[targs], counts["level"], FUN=sum)
  
  levelCounts$Total <- rowSums(levelCounts[,-grep("level", colnames(levelCounts))])
  levelCounts <- levelCounts[order(-levelCounts$Total),]
  
  rownames(levelCounts) <- levelCounts[,"level"]
  levelCounts <- levelCounts[,-1]
  
  levelCounts$Total <- NULL
  mainVal <- colSums(levelCounts)
  miniMat <- sweep(levelCounts, 2, mainVal, FUN = '/')
  
  
  return(miniMat)
}

colMax <- function(data) sapply(data, max, na.rm = TRUE)
colSort <- function(data, ...) sapply(data, sort, ...)

kingTableSmall <- relabundTax(countTableSmall, taxTable, "kingdom")
phylTableSmall <- relabundTax(countTableSmall, taxTable, "phylum")
classTableSmall <- relabundTax(countTableSmall, taxTable, "class")
orderTableSmall <- relabundTax(countTableSmall, taxTable, "order")
linTableSmall <- relabundTax(countTableSmall, taxTable, "lineage")
cladTableSmall <- relabundTax(countTableSmall, taxTable, "clade")
tribTableSmall <- relabundTax(countTableSmall, taxTable, "tribe")

twoParts <- strsplit(colnames(phylTableSmall), "_")
subDepths <- sapply(twoParts, `[[`, 2)

targRow <- unname(unlist(phylTableSmall["Proteobacteria",]))
##Overall mean
mean(targRow)
for (i in 1:length(unique(subDepths))){
  targDep <- unique(subDepths)[i]
  depthMean <- mean(targRow[subDepths == targDep])
  print(c(targDep, depthMean))
}
