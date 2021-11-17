# library(vegan)
# library(cetcolor)
# 
# absSmall <- read.csv("SubsampledData/rrBothSmall.count_table", 
#                      header = TRUE, stringsAsFactors = FALSE, row.names = 1)
# absLarge <- read.csv("SubsampledData/rrBothLarge.count_table", 
#                      header = TRUE, stringsAsFactors = FALSE, row.names = 1)
# 
# remove <- c("Representative_Sequence")
# trimmed.small <- absSmall[,!colnames(absSmall) %in% remove]
# trimmed.large <- absLarge[,!colnames(absLarge) %in% remove]
# 
# relSmall <- trimmed.small/12000
# relLarge <- trimmed.large/12000
# 
# relSmallMat <- t(data.matrix(relSmall))
# relLargelMat <- t(data.matrix(relLarge))
# 
# vare.mds.S <- metaMDS(comm = relSmallMat, trace = FALSE, distance = "jaccard")
# vare.mds.L <- metaMDS(comm = relLargelMat, trace = FALSE, distance = "jaccard")
# 
# depthcolors <- data.frame("Depth" = c("05m", "10m", "DCM", "50m", "90m"), 
#                           "Color" = cet_pal(8, name = "cbl1")[seq(6,2)])
# depthpch <- data.frame("Depth" = c("05m", "10m", "DCM", "50m", "90m"), 
#                        "pch" = c(16, 17, 15, 0, 1))
# monthcolor <- data.frame("Month" = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
#                                      "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"),
#                          "Color" = cet_pal(12, name = "c1"))
# yearcolor <- data.frame("Year" = c("2016", "2017", "2018"), "Color" = 
#                           c("#117733", "#88CCEE", "#882255"))
# createMetaDat <- function(df){
#   samps <- colnames(df)
#   
#   metadat <- data.frame(lapply(1:4,function(i)sapply(strsplit(samps,"_"),"[",i)))
#   names(metadat) <- c("Num", "Date", "Depth", "Size")
#   rownames(metadat) <- samps
#   metadat$Date <- gsub("Sept", "Sep", metadat$Date)
#   metadat$Date <- as.Date(metadat$Date, "%d%b%y")
#   
#   metadat$Month <- format.Date(metadat$Date, format = "%m")
#   metadat$Year <- format.Date(metadat$Date, format = "%Y")
#   metadat$DepthColor <- depthcolors$Color[match(metadat$Depth, depthcolors$Depth)]
#   metadat$MonthColor <- monthcolor$Color[match(
#     format.Date(metadat$Date, "%b"), monthcolor$Month)]
#   metadat$YearColor <- yearcolor$Color[match(
#     format.Date(metadat$Date, "%Y"), yearcolor$Year)]
#   metadat$DepthPch <- depthpch$pch[match(metadat$Depth, depthpch$Depth)]
#   
#   return(metadat)
# }
# 
# metaDatS <- createMetaDat(relSmall)
# metaDatL <- createMetaDat(relLarge)
# 
# save.image(file="~/FlatheadMicrobes/SubsampledData/mdsAndMeta.RData")
load(file="~/FlatheadMicrobes/SubsampledData/mdsAndMeta.RData")

depthpch$Depth <- as.character(depthpch$Depth)
depthpch$Depth <- c("5 m", "10 m", "Chl. Max", "50 m", "90 m")

png("FigBin/mds_by_month.png", width = 500, height = 1000)
plot.new()

par(new = "TRUE",plt = c(0.1,0.95,0.6,0.98),las = 1, cex.axis = 1.5)
ordiplot(vare.mds.S, type = "none", xlab = "", ylab = "", main = "", ylim = c(-1.5, 1), 
         xaxt = "n", yaxt = "n")
points(vare.mds.S, col = as.character(metaDatS$MonthColor), pch = metaDatS$DepthPch, 
       cex = 1.5, lwd = 2)
text(-1.9, -2, "Small (0.2 - 3 μm)", cex = 1.5, font = 2, pos = 4)
text(-1.9, 1.25, "A", cex = 1.5, font = 2, pos = 4)  
axis(1, cex.axis = 1.5, tck = -0.02, at = seq(-2, 2, by = 1), lwd.ticks = 1.5)
axis(2, cex.axis = 1.5, tck = -0.02, at = seq(-2, 2, by = 1), lwd.ticks = 1.5, las = 2)
box(lwd = 1.5)

par(new = "TRUE",plt = c(0.1,0.95,0.18,0.56),las = 1, cex.axis = 1.5)
ordiplot(vare.mds.L, type = "none", main = "", xlab = "", ylab = "", ylim = c(-2, 1.5), 
         xaxt = "n", yaxt = "n")
points(vare.mds.L, col = as.character(metaDatL$MonthColor), pch = metaDatL$DepthPch, 
       cex = 1.5, lwd = 2)
text(-1.55, -1.9, "Large (>3 μm)", cex = 1.5, font = 2, pos = 4)
text(-1.55, 1.25, "B", cex = 1.5, font = 2, pos = 4)
box(lwd = 1.5)

legend("bottomleft", inset=c(0.1,-0.47), 
       legend = monthcolor$Month[c(1,5,9,2,6,10,3,7,11,4,8,12)], 
       fill = as.character(monthcolor$Color)[c(1,5,9,2,6,10,3,7,11,4,8,12)], 
       xpd = NA, cex = 1.5, ncol = 4,
       box.lwd = 0, text.font = 2, title = "")
mtext("Month", side = 1, line = 2, font = 2, cex = 1.5, adj = -0.1, padj = 7.8)

legend("bottomleft", inset=c(0.12,-0.27), legend = depthpch[c(1,4,2,5,3), "Depth"], 
       pch = depthpch[c(1,4,2,5,3), "pch"], pt.lwd = 2, ncol = 3, xpd = NA, cex = 1.5, 
       box.lwd = 0, text.font = 2, title = "")

mtext("Depth", side = 1, line = 1.2, font = 2, cex = 1.5, adj = -0.1, padj = 3.5)
axis(1, cex.axis = 1.5, tck = -0.02, at = seq(-2, 2, by = 1), lwd.ticks = 1.5)
axis(2, cex.axis = 1.5, tck = -0.02, at = seq(-2, 2, by = 1), lwd.ticks = 1.5, las = 2)

dev.off()
