fullFitDat <- readRDS("furhman_satAndSineRawFit.rds")
params_SS <- readRDS("furhman_permutedSatAndSineFitParams.rds")

sineAndSat <- function(alphPB, PBS, alph, omic, interc, days){
  jaccDist <- PBS*(1-exp(-alphPB*days/PBS))+alph*sin(6.283*days/365+omic)+interc
  return(jaccDist)
}
jacList <- readRDS("jacListData.rds")


xmin <- 0.12
xmax <- 0.998
ygap <- 0.02
ytop <- 0.99
ybottom <- 0.08
ywid <- (ytop - ybottom - 4*ygap)/5

textM <- c("5 m", "10 m", "Chl max", "50 m", "90 m")
cexSmall <- 0.9

tiff("FigBinExtras/jacModelFit_Quint.tiff", 
    width = 7, height = 8, pointsize = 12, units = "in", res = 1200)
plot.new()

for (i in 1:5){
  par(new = "TRUE",plt = c(xmin,xmax,ytop-i*ywid-(i-1)*ygap,
                           ytop-(i-1)*ywid-(i-1)*ygap), las = 1)
  plot(names(jacList[[i]]), jacList[[i]], yaxt = "n", xaxt = "n", ylab = "", 
       xlab = "", ylim = c(0.55, 0.95))
  abline(v=365, lty = 2, col = "gray75")
  abline(v=730, lty = 2, col = "gray75")
  lines(names(jacList[[i]]), 
        sineAndSat(mean(params_SS$alphPB[,i]), mean(params_SS$PBS[,i]), 
                   mean(params_SS$alph[,i]), mean(params_SS$omic[,i]), 
                   mean(params_SS$interc[,i]), as.numeric(names(jacList[[i]]))),
        col = "#E69F00", lwd = 3)
  axis(2, tck = -0.03, hadj = 0.5, cex = cexSmall)
  box(lwd = 1)
  text(720, 0.9, textM[i], adj = 1, cex = cexSmall)
  if (i < 5){
    axis(1, at = seq(100, 900, by = 200), tck = -0.03, hadj = 1.2, labels = NA)
    axis(1, tck = -0.05, hadj = 1.2, labels = NA)
  }else {
    axis(1, at = seq(100, 900, by = 200), tck = -0.03, hadj = 1.2, labels = NA)
    axis(1, tck = -0.05, padj = -1, cex = cexSmall)
    mtext("Days between samples", side = 1, line = 2, las = 0)
  }
  if (i == 3){
    mtext("Jaccard distance", side = 2, line = 2, las = 0)
  }
}

dev.off()

modelParams <- data.frame(alphPB <- colMeans(params_SS$alphPB), 
                          alphPB_sd <- apply(params_SS$alphPB, 2, sd),
                          PBS <- colMeans(params_SS$PBS), 
                          PBS_sd <- apply(params_SS$PBS, 2, sd),
                          alph <- colMeans(params_SS$alph), 
                          alph_sd <- apply(params_SS$alph, 2, sd),
                          omic <- colMeans(params_SS$omic), 
                          omic_sd <- apply(params_SS$omic, 2, sd),
                          interc <- colMeans(params_SS$interc),
                          interc_sd <- apply(params_SS$interc, 2, sd))

write.csv(modelParams, file = "~/FlatheadMicrobes/Drafts/modelFitParams.csv")
