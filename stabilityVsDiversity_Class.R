source("~/FlatheadMicrobes/stability_Class.R")
library(vegan)

largeDivInds <- data.frame("Org" = names(largeStabs), "Stabs" = largeStabs,
                           "Abund" = largeAbund,
                           "Shannon.Mean" = NA, "Shannon.Sd" = NA, 
                           "Rich.Mean" = NA, "Rich.Sd" = NA, 
                           "Even.Mean" = NA, "Even.Sd" = NA)
smallDivInds <- data.frame("Org" = names(smallStabs), "Stabs" = smallStabs,
                           "Abund" = smallAbund,
                           "Shannon.Mean" = NA, "Shannon.Sd" = NA, 
                           "Rich.Mean" = NA, "Rich.Sd" = NA, 
                           "Even.Mean" = NA, "Even.Sd" = NA)

for (i in 1:length(largeClass)){
  df <- largeClass[[i]]
  df <- t(df)
  shannon <- diversity(df, "shannon")
  largeDivInds[i, "Shannon.Sd"] <- sd(shannon)
  largeDivInds[i, "Shannon.Mean"] <- mean(shannon)
  
  rich <- specnumber(df)
  largeDivInds[i, "Rich.Sd"] <- sd(rich)
  largeDivInds[i, "Rich.Mean"] <- mean(rich)
  
  evenness <- shannon/log(rich)
  evenness <- evenness[!is.nan(evenness)]
  largeDivInds[i, "Even.Sd"] <- sd(evenness)
  largeDivInds[i, "Even.Mean"] <- mean(evenness)
}

for (i in 1:length(smallClass)){
  df <- smallClass[[i]]
  df <- t(df)
  shannon <- diversity(df, "shannon")
  smallDivInds[i, "Shannon.Sd"] <- sd(shannon)
  smallDivInds[i, "Shannon.Mean"] <- mean(shannon)
  
  rich <- specnumber(df)
  smallDivInds[i, "Rich.Sd"] <- sd(rich)
  smallDivInds[i, "Rich.Mean"] <- mean(rich)
  
  evenness <- shannon/log(rich)
  evenness <- evenness[!is.nan(evenness)]
  smallDivInds[i, "Even.Sd"] <- sd(evenness)
  smallDivInds[i, "Even.Mean"] <- mean(evenness)
}

largeDivInds$Org <- as.character(largeDivInds$Org)
largeDivInds[1:(nrow(largeDivInds)-numEuksLarge),"Org"] <- "#1A85FF"
largeDivInds[(nrow(largeDivInds)-numEuksLarge+1):nrow(largeDivInds),"Org"] <- "#D41159"
smallDivInds$Org <- as.character(smallDivInds$Org)
smallDivInds[1:(nrow(smallDivInds)-numEuksSmall),"Org"] <- "#1A85FF"
smallDivInds[(nrow(smallDivInds)-numEuksSmall+1):nrow(smallDivInds),"Org"] <- "#D41159"

largeDivInds <- largeDivInds[order(-largeDivInds$Abund),]
smallDivInds <- smallDivInds[order(-smallDivInds$Abund),]

largeDivInds <- largeDivInds[grep("unclassified", rownames(largeDivInds), invert = TRUE),]
smallDivInds <- smallDivInds[grep("unclassified", rownames(smallDivInds), invert = TRUE),]

largeDivIndTrimmed <- largeDivInds[1:10,]
smallDivIndTrimmed <- smallDivInds[1:10,]


png("FigBinExtras/ClassDS/allDivStabs.png", width = 2000, height = 2690)
plot.new()

par(new = "TRUE",plt = c(0.06,0.45,0.68,0.97),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(smallDivIndTrimmed$Stabs, smallDivIndTrimmed$Shannon.Mean, xlab = "", ylab = "",
     cex = smallDivIndTrimmed$Abund/50, col = smallDivIndTrimmed$Org,
     xlim = c(0.5, 3.5), xaxt = "n", lwd = 2, ylim = c(1, 3.5), yaxt = "n")
x_vec <- smallDivIndTrimmed$Stabs
x_vec[rownames(smallDivIndTrimmed) %in% c("Planctomycetacia", "Verrucomicrobiae")] <- 
  x_vec[rownames(smallDivIndTrimmed) %in% c("Planctomycetacia", "Verrucomicrobiae")]-0.45
x_vec[rownames(smallDivIndTrimmed) %in% c("Acidimicrobiia")] <- 
  x_vec[rownames(smallDivIndTrimmed) %in% c("Acidimicrobiia")]-0.35
x_vec[rownames(smallDivIndTrimmed) %in% c("Actinobacteria", "Gammaproteobacteria")] <- 
  x_vec[rownames(smallDivIndTrimmed) %in% c("Actinobacteria", "Gammaproteobacteria")]-0.15
x_vec[rownames(smallDivIndTrimmed) %in% c("Alphaproteobacteria")] <- 
  x_vec[rownames(smallDivIndTrimmed) %in% c("Alphaproteobacteria")]+0.5
y_vec <- smallDivIndTrimmed$Shannon.Mean
y_vec[rownames(smallDivIndTrimmed) %in% c("Gammaproteobacteria")] <- 
  y_vec[rownames(smallDivIndTrimmed) %in% c("Gammaproteobacteria")]-0.9
y_vec[rownames(smallDivIndTrimmed) %in% c("Bacteroidia")] <- 
  y_vec[rownames(smallDivIndTrimmed) %in% c("Bacteroidia")]-0.6
y_vec[rownames(smallDivIndTrimmed) %in% c("Acidimicrobiia")] <- 
  y_vec[rownames(smallDivIndTrimmed) %in% c("Acidimicrobiia")]-0.1
y_vec[rownames(smallDivIndTrimmed) %in% c("Phycisphaerae", "Cryptophyceae")] <- 
  y_vec[rownames(smallDivIndTrimmed) %in% c("Phycisphaerae", "Cryptophyceae")]+0.06
y_vec[rownames(smallDivIndTrimmed) %in% c("Actinobacteria")] <- 
  y_vec[rownames(smallDivIndTrimmed) %in% c("Actinobacteria")]+0.15
y_vec[rownames(smallDivIndTrimmed) %in% c("Oxyphotobacteria")] <- 
  y_vec[rownames(smallDivIndTrimmed) %in% c("Oxyphotobacteria")]+0.2
text(x_vec, y_vec, labels=rownames(smallDivIndTrimmed),
     cex= 1.5)
text(0.5, 3.5, "A", cex= 2.5, font = 2)
box(lwd = 2.5)
arrows(3.5, 2.12, 3.45, 2.82, lwd = 2, length = 0.15)
arrows(3.2, 2.3, 3.2, 2.65, lwd = 2, length = 0.15)
mtext(bquote(bold("Shannon's Diversity Index")), 2, font = 2, line = 5.5, 
      cex = 2.5, las = 0)
mtext("Small (0.2-3 μm)", 3, font = 2, line = 2, cex = 2.5)
axis(1, cex.axis=2.5, tck = -0.03, lwd.ticks = 2.5, padj = 1.2, font = 2, 
     at = seq(1, 6, by = 1), labels = NA)
axis(1, cex.axis=2.5, tck = -0.015, lwd.ticks = 2.5, padj = 1.2, font = 2, 
     at = seq(0.5, 6.5, by = 1), labels = NA)
axis(2, cex.axis=2.5, tck = -0.025, lwd.ticks = 2.5, hadj = 1.2, font = 2, 
     at = seq(1, 6, by = 1), labels = seq(1, 6, by = 1))
axis(2, cex.axis=2.5, tck = -0.015, lwd.ticks = 2.5, hadj = 1.2, font = 2, 
     at = seq(0.5, 6.5, by = 1), labels = NA)

par(new = "TRUE",plt = c(0.49,0.88,0.68,0.97),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(largeDivIndTrimmed$Stabs, largeDivIndTrimmed$Shannon.Mean, xlab = "", ylab = "",
     cex = largeDivIndTrimmed$Abund/50, col = largeDivIndTrimmed$Org,
     xlim = c(0.5, 3.5), xaxt = "n", lwd = 2, ylim = c(1, 3.5), yaxt = "n")
x_vec <- largeDivIndTrimmed$Stabs
x_vec[rownames(largeDivIndTrimmed) %in% c("Deltaproteobacteria")] <- 
  x_vec[rownames(largeDivIndTrimmed) %in% c("Deltaproteobacteria")]-0.6
x_vec[rownames(largeDivIndTrimmed) %in% c("Chrysophyceae")] <- 
  x_vec[rownames(largeDivIndTrimmed) %in% c("Chrysophyceae")]+0.15
x_vec[rownames(largeDivIndTrimmed) %in% c("Verrucomicrobiae")] <- 
  x_vec[rownames(largeDivIndTrimmed) %in% c("Verrucomicrobiae")]+0.45
y_vec <- largeDivIndTrimmed$Shannon.Mean
y_vec[rownames(largeDivIndTrimmed) %in% c("Bacillariophyta")] <- 
  y_vec[rownames(largeDivIndTrimmed) %in% c("Bacillariophyta")]-0.3
y_vec[rownames(largeDivIndTrimmed) %in% 
        c("Chrysophyceae", "Phycisphaerae", "Actinobacteria", 
          "Gammaproteobacteria", "Deltaproteobacteria")] <- 
  y_vec[rownames(largeDivIndTrimmed) %in% 
          c("Chrysophyceae", "Phycisphaerae", "Actinobacteria", 
            "Gammaproteobacteria", "Deltaproteobacteria")]+0.14
y_vec[rownames(largeDivIndTrimmed) %in% c("Planctomycetacia")] <- 
  y_vec[rownames(largeDivIndTrimmed) %in% c("Planctomycetacia")]+0.25
text(x_vec, y_vec, labels=rownames(largeDivIndTrimmed),
     cex= 1.5)
text(0.5, 3.5, "B", cex= 2.5, font = 2)
box(lwd = 2.5)
arrows(1.16, 2.85, 1.32, 2.75, lwd = 2, length = 0.15)
mtext("Large (>3 μm)", 3, font = 2, line = 2, cex = 2.5)
axis(1, cex.axis=2.5, tck = -0.03, lwd.ticks = 2.5, padj = 1.2, font = 2, 
     at = seq(1, 6, by = 1), labels = NA)
axis(1, cex.axis=2.5, tck = -0.015, lwd.ticks = 2.5, padj = 1.2, font = 2, 
     at = seq(0.5, 6.5, by = 1), labels = NA)
axis(2, cex.axis=2.5, tck = -0.025, lwd.ticks = 2.5, hadj = 1.2, font = 2, 
     at = seq(1, 6, by = 1), labels = seq(1, 6, by = 1))
axis(2, cex.axis=2.5, tck = -0.015, lwd.ticks = 2.5, hadj = 1.2, font = 2, 
     at = seq(0.5, 6.5, by = 1), labels = NA)


par(new = "TRUE",plt = c(0.06,0.45,0.36,0.65),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(smallDivIndTrimmed$Stabs, smallDivIndTrimmed$Rich.Mean, 
     cex = smallDivIndTrimmed$Abund/50, xlab = "", ylab = "", lwd = 2, 
     col = smallDivIndTrimmed$Org, xlim = c(0.5, 3.5), xaxt = "n", ylim = c(0, 440),
     yaxt = "n")
x_vec <- smallDivIndTrimmed$Stabs
x_vec[rownames(smallDivIndTrimmed) %in% c("Alphaproteobacteria")] <- 
  x_vec[rownames(smallDivIndTrimmed) %in% c("Alphaproteobacteria")]-0.3
x_vec[rownames(smallDivIndTrimmed) %in% c("Phycisphaerae")] <- 
  x_vec[rownames(smallDivIndTrimmed) %in% c("Phycisphaerae")]-0.2
x_vec[rownames(smallDivIndTrimmed) %in% c("Gammaproteobacteria")] <- 
  x_vec[rownames(smallDivIndTrimmed) %in% c("Gammaproteobacteria")]-0.1
x_vec[rownames(smallDivIndTrimmed) %in% c("Planctomycetacia", "Verrucomicrobiae")] <- 
  x_vec[rownames(smallDivIndTrimmed) %in% c("Planctomycetacia", "Verrucomicrobiae")]+0.45
y_vec <- smallDivIndTrimmed$Rich.Mean
y_vec[rownames(smallDivIndTrimmed) %in% c("Bacteroidia")] <- 
  y_vec[rownames(smallDivIndTrimmed) %in% c("Bacteroidia")]-40
y_vec[rownames(smallDivIndTrimmed) %in% 
        c("Cryptophyceae", "Phycisphaerae", "Planctomycetacia")] <- 
  y_vec[rownames(smallDivIndTrimmed) %in% 
          c("Cryptophyceae", "Phycisphaerae", "Planctomycetacia")]-10
y_vec[rownames(smallDivIndTrimmed) %in% c("Acidimicrobiia", "Alphaproteobacteria")] <- 
  y_vec[rownames(smallDivIndTrimmed) %in% c("Acidimicrobiia", "Alphaproteobacteria")]+30
y_vec[rownames(smallDivIndTrimmed) %in% c("Oxyphotobacteria")] <- 
  y_vec[rownames(smallDivIndTrimmed) %in% c("Oxyphotobacteria")]+35
y_vec[rownames(smallDivIndTrimmed) %in% c("Gammaproteobacteria")] <- 
  y_vec[rownames(smallDivIndTrimmed) %in% c("Gammaproteobacteria")]+40
text(x_vec, y_vec, labels=rownames(smallDivIndTrimmed),
     cex= 1.5)
text(0.5, 440, "C", cex= 2.5, font = 2)
arrows(1.72, 90, 1.9, 75, lwd = 2, length = 0.15)
box(lwd = 2.5)
mtext(bquote(bold("Richness")), 2, font = 2, line = 5.5, 
      cex = 2.5, las = 0)
axis(1, cex.axis=2.5, tck = -0.03, lwd.ticks = 2.5, padj = 1.2, font = 2, 
     at = seq(1, 6, by = 1), labels = NA)
axis(1, cex.axis=2.5, tck = -0.015, lwd.ticks = 2.5, padj = 1.2, font = 2, 
     at = seq(0.5, 6.5, by = 1), labels = NA)
axis(2, cex.axis=2.5, tck = -0.025, lwd.ticks = 2.5, hadj = 1.2, font = 2, 
     at = seq(0, 400, by = 200), labels = seq(0, 400, by = 200))
axis(2, cex.axis=2.5, tck = -0.015, lwd.ticks = 2.5, hadj = 1.2, font = 2, 
     at = seq(100, 500.8, by = 100), labels = NA)

par(new = "TRUE",plt = c(0.49,0.88,0.36,0.65),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(largeDivIndTrimmed$Stabs, largeDivIndTrimmed$Rich.Mean, 
     cex = largeDivIndTrimmed$Abund/50, xlab = "", ylab = "", lwd = 2, 
     col = largeDivIndTrimmed$Org, xlim = c(0.5, 3.5), xaxt = "n", ylim = c(0, 110),
     yaxt = "n")
x_vec <- largeDivIndTrimmed$Stabs
x_vec[rownames(largeDivIndTrimmed) %in% c("Oxyphotobacteria")] <- 
  x_vec[rownames(largeDivIndTrimmed) %in% c("Oxyphotobacteria")]-0.8
x_vec[rownames(largeDivIndTrimmed) %in% c("Phycisphaerae")] <- 
  x_vec[rownames(largeDivIndTrimmed) %in% c("Phycisphaerae")]-0.4
x_vec[rownames(largeDivIndTrimmed) %in% c("Chrysophyceae")] <- 
  x_vec[rownames(largeDivIndTrimmed) %in% c("Chrysophyceae")]+0.1
x_vec[rownames(largeDivIndTrimmed) %in% 
        c("Verrucomicrobiae", "Gammaproteobacteria", "Deltaproteobacteria")] <- 
  x_vec[rownames(largeDivIndTrimmed) %in% 
          c("Verrucomicrobiae", "Gammaproteobacteria", "Deltaproteobacteria")]+0.42
x_vec[rownames(largeDivIndTrimmed) %in% c("Planctomycetacia")] <- 
  x_vec[rownames(largeDivIndTrimmed) %in% c("Planctomycetacia")]+0.8
y_vec <- largeDivIndTrimmed$Rich.Mean
y_vec[rownames(largeDivIndTrimmed) %in% c("Deltaproteobacteria")] <- 
  y_vec[rownames(largeDivIndTrimmed) %in% c("Deltaproteobacteria")]-4
y_vec[rownames(largeDivIndTrimmed) %in% 
        c("Actinobacteria", "Gammaproteobacteria", "Verrucomicrobiae")] <- 
  y_vec[rownames(largeDivIndTrimmed) %in% 
          c("Actinobacteria", "Gammaproteobacteria", "Verrucomicrobiae")]+5
y_vec[rownames(largeDivIndTrimmed) %in% 
        c("Chrysophyceae")] <- 
  y_vec[rownames(largeDivIndTrimmed) %in% 
          c("Chrysophyceae")]+7
text(x_vec, y_vec, labels=rownames(largeDivIndTrimmed),
     cex= 1.5)
text(0.5, 110, "D", cex= 2.5, font = 2)
arrows(2.09, 74, 1.86, 74, lwd = 2, length = 0.15)
arrows(2, 92, 1.86, 90, lwd = 2, length = 0.15)
box(lwd = 2.5)
axis(1, cex.axis=2.5, tck = -0.03, lwd.ticks = 2.5, padj = 1.2, font = 2, 
     at = seq(1, 6, by = 1), labels = NA)
axis(1, cex.axis=2.5, tck = -0.015, lwd.ticks = 2.5, padj = 1.2, font = 2, 
     at = seq(0.5, 6.5, by = 1), labels = NA)
axis(2, cex.axis=2.5, tck = -0.025, lwd.ticks = 2.5, hadj = 1.2, font = 2, 
     at = seq(0, 100, by = 50), labels = seq(0, 100, by = 50))
axis(2, cex.axis=2.5, tck = -0.015, lwd.ticks = 2.5, hadj = 1.2, font = 2, 
     at = seq(25, 125, by = 50), labels = NA)

legend("topright", legend = c("Prokaryotes", "Ph. Euks"), fill = c("#1A85FF", "#D41159"), 
       inset = c(-0.28, 0), xpd = TRUE, box.lwd = 2.5, cex = 2)
par(xpd = TRUE)
points(x = rep(4.1, 3), y = c(65, 38, 20), pch = 1, lwd = 2, cex = c(1200, 600, 120)/50)
text(4.1, 90, "Average %", cex = 2, font = 2)
text(4.1, 85, "Abundance", cex = 2, font = 2)
text(4.1, 78, "10%", cex = 2)
text(4.1, 46, "5%", cex = 2)
text(4.1, 25, "1%", cex = 2)



par(new = "TRUE",plt = c(0.06,0.45,0.04,0.33),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(smallDivIndTrimmed$Stabs, smallDivIndTrimmed$Even.Mean, 
     cex = smallDivIndTrimmed$Abund/50, xlab = "", ylab = "", lwd = 2, 
     col = smallDivIndTrimmed$Org, xlim = c(0.5, 3.5), xaxt = "n", ylim = c(0.3, 0.8), 
     yaxt = "n")
x_vec <- smallDivIndTrimmed$Stabs
x_vec[rownames(smallDivIndTrimmed) %in% c("Chrysophyceae")] <- 
  x_vec[rownames(smallDivIndTrimmed) %in% c("Chrysophyceae")]+0.15
x_vec[rownames(smallDivIndTrimmed) %in% c("Alphaproteobacteria", "Planctomycetacia")] <- 
  x_vec[rownames(smallDivIndTrimmed) %in% c("Alphaproteobacteria", "Planctomycetacia")]-0.48
x_vec[rownames(smallDivIndTrimmed) %in% c("Bacteroidia", "Acidimicrobiia")] <- 
  x_vec[rownames(smallDivIndTrimmed) %in% c("Bacteroidia", "Acidimicrobiia")]-0.4
x_vec[rownames(smallDivIndTrimmed) %in% c("Gammaproteobacteria")] <- 
  x_vec[rownames(smallDivIndTrimmed) %in% c("Gammaproteobacteria")]-0.1
y_vec <- smallDivIndTrimmed$Even.Mean
y_vec[rownames(smallDivIndTrimmed) %in% 
        c("Thermoleophilia", "Deltaproteobacteria", "Ignavibacteria", "Holophagae", 
          "Cryptophyceae")] <- 
  y_vec[rownames(smallDivIndTrimmed) %in% 
          c("Thermoleophilia", "Deltaproteobacteria", "Ignavibacteria", "Holophagae",
            "Cryptophyceae")]+0.011
y_vec[rownames(smallDivIndTrimmed) %in% c("Phycisphaerae", "Bacteroidia")] <- 
  y_vec[rownames(smallDivIndTrimmed) %in% c("Phycisphaerae", "Bacteroidia")]+0.02
y_vec[rownames(smallDivIndTrimmed) %in% c("Verrucomicrobiae")] <- 
  y_vec[rownames(smallDivIndTrimmed) %in% c("Verrucomicrobiae")]+0.035
y_vec[rownames(smallDivIndTrimmed) %in% c("Oxyphotobacteria")] <- 
  y_vec[rownames(smallDivIndTrimmed) %in% c("Oxyphotobacteria")]+0.045
y_vec[rownames(smallDivIndTrimmed) %in% c("Gammaproteobacteria")] <- 
  y_vec[rownames(smallDivIndTrimmed) %in% c("Gammaproteobacteria")]+0.1
y_vec[rownames(smallDivIndTrimmed) %in% c("Chrysophyceae")] <- 
  y_vec[rownames(smallDivIndTrimmed) %in% c("Chrysophyceae")]-0.015
y_vec[rownames(smallDivIndTrimmed) %in% c("Alphaproteobacteria")] <- 
  y_vec[rownames(smallDivIndTrimmed) %in% c("Alphaproteobacteria")]-0.01
arrows(3.4, 0.745, 3.4, 0.69, lwd = 2, length = 0.15)
mtext(bquote(bold("S"["T"])), 1, font = 2, line = 6, 
      cex = 2.5, las = 0, adj = 1.08)
text(x_vec, y_vec, labels=rownames(smallDivIndTrimmed),
     cex= 1.5)
text(0.5, 0.8, "E", cex= 2.5, font = 2)
box(lwd = 2.5)
mtext(bquote(bold("Evenness")), 2, font = 2, line = 5.5, 
      cex = 2.5, las = 0)
axis(1, cex.axis=2.5, tck = -0.03, lwd.ticks = 2.5, padj = 1.2, font = 2, 
     at = seq(1, 6, by = 1), labels = seq(1, 6, by = 1))
axis(1, cex.axis=2.5, tck = -0.015, lwd.ticks = 2.5, padj = 1.2, font = 2, 
     at = seq(0.5, 6.5, by = 1), labels = NA)
axis(2, cex.axis=2.5, tck = -0.025, lwd.ticks = 2.5, hadj = 1.2, font = 2, 
     at = seq(0.1, 0.9, by = 0.2), labels = seq(0.1, 0.9, by = 0.2))
axis(2, cex.axis=2.5, tck = -0.015, lwd.ticks = 2.5, hadj = 1.2, font = 2, 
     at = seq(0.2, 0.8, by = 0.2), labels = NA)

par(new = "TRUE",plt = c(0.49,0.88,0.04,0.33),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(largeDivIndTrimmed$Stabs, largeDivIndTrimmed$Even.Mean, 
     cex = largeDivIndTrimmed$Abund/50, xlab = "", ylab = "", lwd = 2, 
     col = largeDivIndTrimmed$Org, xlim = c(0.5, 3.5), xaxt = "n", ylim = c(0.3, 0.8), 
     yaxt = "n")
x_vec <- largeDivIndTrimmed$Stabs
x_vec[rownames(largeDivIndTrimmed) %in% c("Prymnesiophyceae")] <- 
  x_vec[rownames(largeDivIndTrimmed) %in% c("Prymnesiophyceae")]-0.35
x_vec[rownames(largeDivIndTrimmed) %in% c("Cryptophyceae", "Acidimicrobiia")] <- 
  x_vec[rownames(largeDivIndTrimmed) %in% c("Cryptophyceae", "Acidimicrobiia")]+0.63
x_vec[rownames(largeDivIndTrimmed) %in% c("Phycisphaerae")] <- 
  x_vec[rownames(largeDivIndTrimmed) %in% c("Phycisphaerae")]+0.5
x_vec[rownames(largeDivIndTrimmed) %in% c("Deltaproteobacteria")] <- 
  x_vec[rownames(largeDivIndTrimmed) %in% c("Deltaproteobacteria")]-0.1
x_vec[rownames(largeDivIndTrimmed) %in% c("Chrysophyceae", "Verrucomicrobiae")] <- 
  x_vec[rownames(largeDivIndTrimmed) %in% c("Chrysophyceae", "Verrucomicrobiae")]+0.12
y_vec <- largeDivIndTrimmed$Even.Mean
y_vec[rownames(largeDivIndTrimmed) %in% c("Prymnesiophyceae")] <- 
  y_vec[rownames(largeDivIndTrimmed) %in% c("Prymnesiophyceae")]+0.011
y_vec[rownames(largeDivIndTrimmed) %in% c("Cryptophyceae")] <- 
  y_vec[rownames(largeDivIndTrimmed) %in% c("Cryptophyceae")]-0.045
y_vec[rownames(largeDivIndTrimmed) %in% c("Anaerolineae")] <- 
  y_vec[rownames(largeDivIndTrimmed) %in% c("Anaerolineae")]+0.02
y_vec[rownames(largeDivIndTrimmed) %in% 
        c("Alphaproteobacteria", "Actinobacteria")] <- 
  y_vec[rownames(largeDivIndTrimmed) %in% 
          c("Alphaproteobacteria", "Actinobacteria")]+0.03
y_vec[rownames(largeDivIndTrimmed) %in% c("Phycisphaerae")] <- 
  y_vec[rownames(largeDivIndTrimmed) %in% c("Phycisphaerae")]-0.03
y_vec[rownames(largeDivIndTrimmed) %in% 
        c("Chrysophyceae", "Gammaproteobacteria", "Verrucomicrobiae", "Acidimicrobiia")] <- 
  y_vec[rownames(largeDivIndTrimmed) %in% 
          c("Chrysophyceae", "Gammaproteobacteria", "Verrucomicrobiae", 
            "Acidimicrobiia")]+0.035
y_vec[rownames(largeDivIndTrimmed) %in% c("Deltaproteobacteria")] <- 
  y_vec[rownames(largeDivIndTrimmed) %in% c("Deltaproteobacteria")]+0.045
text(x_vec, y_vec, labels=rownames(largeDivIndTrimmed),
     cex= 1.5)
text(0.5, 0.8, "F", cex= 2.5, font = 2)
arrows(1.4, 0.76, 1.4, 0.745, lwd = 2, length = 0.15)
arrows(1.59, 0.6, 1.45, 0.615, lwd = 2, length = 0.15)
box(lwd = 2.5)
axis(1, cex.axis=2.5, tck = -0.03, lwd.ticks = 2.5, padj = 1.2, font = 2, 
     at = seq(1, 6, by = 1), labels = seq(1, 6, by = 1))
axis(1, cex.axis=2.5, tck = -0.015, lwd.ticks = 2.5, padj = 1.2, font = 2, 
     at = seq(0.5, 6.5, by = 1), labels = NA)
axis(2, cex.axis=2.5, tck = -0.025, lwd.ticks = 2.5, hadj = 1.2, font = 2, 
     at = seq(0.1, 0.9, by = 0.2), labels = seq(0.1, 0.9, by = 0.2))
axis(2, cex.axis=2.5, tck = -0.015, lwd.ticks = 2.5, hadj = 1.2, font = 2, 
     at = seq(0.2, 0.8, by = 0.2), labels = NA)


dev.off()




png("FigBinExtras/ClassDS/evenness.png", width = 1800, height = 945)
plot.new()
par(new = "TRUE",plt = c(0.06,0.48,0.1,0.9),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(smallDivIndTrimmed$Stabs, smallDivIndTrimmed$Even.Mean, 
     cex = smallDivIndTrimmed$Abund/50, xlab = "", ylab = "", lwd = 2, 
     col = smallDivIndTrimmed$Org, xlim = c(0.5, 3.5), xaxt = "n", ylim = c(0.3, 0.8), 
     yaxt = "n")
x_vec <- smallDivIndTrimmed$Stabs
x_vec[rownames(smallDivIndTrimmed) %in% c("Chrysophyceae")] <- 
  x_vec[rownames(smallDivIndTrimmed) %in% c("Chrysophyceae")]+0.15
x_vec[rownames(smallDivIndTrimmed) %in% c("Alphaproteobacteria", "Planctomycetacia")] <- 
  x_vec[rownames(smallDivIndTrimmed) %in% c("Alphaproteobacteria", "Planctomycetacia")]-0.48
x_vec[rownames(smallDivIndTrimmed) %in% c("Bacteroidia", "Acidimicrobiia")] <- 
  x_vec[rownames(smallDivIndTrimmed) %in% c("Bacteroidia", "Acidimicrobiia")]-0.4
x_vec[rownames(smallDivIndTrimmed) %in% c("Gammaproteobacteria")] <- 
  x_vec[rownames(smallDivIndTrimmed) %in% c("Gammaproteobacteria")]-0.1
y_vec <- smallDivIndTrimmed$Even.Mean
y_vec[rownames(smallDivIndTrimmed) %in% 
        c("Thermoleophilia", "Deltaproteobacteria", "Ignavibacteria", "Holophagae", 
          "Cryptophyceae")] <- 
  y_vec[rownames(smallDivIndTrimmed) %in% 
          c("Thermoleophilia", "Deltaproteobacteria", "Ignavibacteria", "Holophagae",
            "Cryptophyceae")]+0.011
y_vec[rownames(smallDivIndTrimmed) %in% c("Phycisphaerae", "Bacteroidia")] <- 
  y_vec[rownames(smallDivIndTrimmed) %in% c("Phycisphaerae", "Bacteroidia")]+0.02
y_vec[rownames(smallDivIndTrimmed) %in% c("Verrucomicrobiae")] <- 
  y_vec[rownames(smallDivIndTrimmed) %in% c("Verrucomicrobiae")]+0.035
y_vec[rownames(smallDivIndTrimmed) %in% c("Oxyphotobacteria")] <- 
  y_vec[rownames(smallDivIndTrimmed) %in% c("Oxyphotobacteria")]+0.045
y_vec[rownames(smallDivIndTrimmed) %in% c("Gammaproteobacteria")] <- 
  y_vec[rownames(smallDivIndTrimmed) %in% c("Gammaproteobacteria")]+0.1
y_vec[rownames(smallDivIndTrimmed) %in% c("Chrysophyceae")] <- 
  y_vec[rownames(smallDivIndTrimmed) %in% c("Chrysophyceae")]-0.015
y_vec[rownames(smallDivIndTrimmed) %in% c("Alphaproteobacteria")] <- 
  y_vec[rownames(smallDivIndTrimmed) %in% c("Alphaproteobacteria")]-0.01
arrows(3.4, 0.745, 3.4, 0.7, lwd = 2, length = 0.15)
text(x_vec, y_vec, labels=rownames(smallDivIndTrimmed),
     cex= 1.5)
box(lwd = 2.5)
mtext(bquote(bold("S"["T"]~"(Lehman & Tilman, 2000)")), 1, font = 2, line = 5, 
      cex = 2.5, las = 0, adj = 2)
mtext(bquote(bold("Evenness")), 2, font = 2, line = 5, 
      cex = 2.5, las = 0)
mtext("Small (0.2-3 μm)", 3, font = 2, line = 2, cex = 2.5)
axis(1, cex.axis=2.5, tck = -0.03, lwd.ticks = 2.5, padj = 1.2, font = 2, 
     at = seq(1, 6, by = 1), labels = seq(1, 6, by = 1))
axis(1, cex.axis=2.5, tck = -0.015, lwd.ticks = 2.5, padj = 1.2, font = 2, 
     at = seq(0.5, 6.5, by = 1), labels = NA)
axis(2, cex.axis=2.5, tck = -0.025, lwd.ticks = 2.5, hadj = 1.2, font = 2, 
     at = seq(0.1, 0.9, by = 0.2), labels = seq(0.1, 0.9, by = 0.2))
axis(2, cex.axis=2.5, tck = -0.015, lwd.ticks = 2.5, hadj = 1.2, font = 2, 
     at = seq(0.2, 0.8, by = 0.2), labels = NA)

par(new = "TRUE",plt = c(0.56,0.98,0.1,0.9),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(largeDivIndTrimmed$Stabs, largeDivIndTrimmed$Even.Mean, 
     cex = largeDivIndTrimmed$Abund/50, xlab = "", ylab = "", lwd = 2, 
     col = largeDivIndTrimmed$Org, xlim = c(0.5, 3.5), xaxt = "n", ylim = c(0.3, 0.8), 
     yaxt = "n")
x_vec <- largeDivIndTrimmed$Stabs
x_vec[rownames(largeDivIndTrimmed) %in% c("Prymnesiophyceae")] <- 
  x_vec[rownames(largeDivIndTrimmed) %in% c("Prymnesiophyceae")]-0.35
x_vec[rownames(largeDivIndTrimmed) %in% c("Cryptophyceae", "Acidimicrobiia")] <- 
  x_vec[rownames(largeDivIndTrimmed) %in% c("Cryptophyceae", "Acidimicrobiia")]+0.63
x_vec[rownames(largeDivIndTrimmed) %in% c("Phycisphaerae")] <- 
  x_vec[rownames(largeDivIndTrimmed) %in% c("Phycisphaerae")]+0.5
x_vec[rownames(largeDivIndTrimmed) %in% c("Deltaproteobacteria")] <- 
  x_vec[rownames(largeDivIndTrimmed) %in% c("Deltaproteobacteria")]-0.1
x_vec[rownames(largeDivIndTrimmed) %in% c("Chrysophyceae", "Verrucomicrobiae")] <- 
  x_vec[rownames(largeDivIndTrimmed) %in% c("Chrysophyceae", "Verrucomicrobiae")]+0.12
y_vec <- largeDivIndTrimmed$Even.Mean
y_vec[rownames(largeDivIndTrimmed) %in% c("Prymnesiophyceae")] <- 
  y_vec[rownames(largeDivIndTrimmed) %in% c("Prymnesiophyceae")]+0.011
y_vec[rownames(largeDivIndTrimmed) %in% c("Cryptophyceae")] <- 
  y_vec[rownames(largeDivIndTrimmed) %in% c("Cryptophyceae")]-0.045
y_vec[rownames(largeDivIndTrimmed) %in% c("Anaerolineae")] <- 
  y_vec[rownames(largeDivIndTrimmed) %in% c("Anaerolineae")]+0.02
y_vec[rownames(largeDivIndTrimmed) %in% 
        c("Alphaproteobacteria", "Actinobacteria")] <- 
  y_vec[rownames(largeDivIndTrimmed) %in% 
          c("Alphaproteobacteria", "Actinobacteria")]+0.03
y_vec[rownames(largeDivIndTrimmed) %in% c("Phycisphaerae")] <- 
  y_vec[rownames(largeDivIndTrimmed) %in% c("Phycisphaerae")]-0.03
y_vec[rownames(largeDivIndTrimmed) %in% 
        c("Chrysophyceae", "Gammaproteobacteria", "Verrucomicrobiae", "Acidimicrobiia")] <- 
  y_vec[rownames(largeDivIndTrimmed) %in% 
          c("Chrysophyceae", "Gammaproteobacteria", "Verrucomicrobiae", 
            "Acidimicrobiia")]+0.035
y_vec[rownames(largeDivIndTrimmed) %in% c("Deltaproteobacteria")] <- 
  y_vec[rownames(largeDivIndTrimmed) %in% c("Deltaproteobacteria")]+0.045
text(x_vec, y_vec, labels=rownames(largeDivIndTrimmed),
     cex= 1.5)
mtext("Large (>3 μm)", 3, font = 2, line = 2, cex = 2.5)
arrows(1.4, 0.76, 1.4, 0.745, lwd = 2, length = 0.15)
arrows(1.59, 0.6, 1.45, 0.615, lwd = 2, length = 0.15)
box(lwd = 2.5)
axis(1, cex.axis=2.5, tck = -0.03, lwd.ticks = 2.5, padj = 1.2, font = 2, 
     at = seq(1, 6, by = 1), labels = seq(1, 6, by = 1))
axis(1, cex.axis=2.5, tck = -0.015, lwd.ticks = 2.5, padj = 1.2, font = 2, 
     at = seq(0.5, 6.5, by = 1), labels = NA)
axis(2, cex.axis=2.5, tck = -0.025, lwd.ticks = 2.5, hadj = 1.2, font = 2, 
     at = seq(0.1, 0.9, by = 0.2), labels = seq(0.1, 0.9, by = 0.2))
axis(2, cex.axis=2.5, tck = -0.015, lwd.ticks = 2.5, hadj = 1.2, font = 2, 
     at = seq(0.2, 0.8, by = 0.2), labels = NA)

dev.off()


png("FigBinExtras/ClassDS/richness.png", width = 1800, height = 945)
plot.new()
par(new = "TRUE",plt = c(0.06,0.48,0.1,0.9),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(smallDivIndTrimmed$Stabs, smallDivIndTrimmed$Rich.Mean, 
     cex = smallDivIndTrimmed$Abund/50, xlab = "", ylab = "", lwd = 2, 
     col = smallDivIndTrimmed$Org, xlim = c(0.5, 3.5), xaxt = "n", ylim = c(0, 440),
     yaxt = "n")
x_vec <- smallDivIndTrimmed$Stabs
x_vec[rownames(smallDivIndTrimmed) %in% c("Alphaproteobacteria")] <- 
  x_vec[rownames(smallDivIndTrimmed) %in% c("Alphaproteobacteria")]-0.3
x_vec[rownames(smallDivIndTrimmed) %in% c("Phycisphaerae")] <- 
  x_vec[rownames(smallDivIndTrimmed) %in% c("Phycisphaerae")]-0.2
x_vec[rownames(smallDivIndTrimmed) %in% c("Gammaproteobacteria")] <- 
  x_vec[rownames(smallDivIndTrimmed) %in% c("Gammaproteobacteria")]-0.1
x_vec[rownames(smallDivIndTrimmed) %in% c("Planctomycetacia", "Verrucomicrobiae")] <- 
  x_vec[rownames(smallDivIndTrimmed) %in% c("Planctomycetacia", "Verrucomicrobiae")]+0.45
y_vec <- smallDivIndTrimmed$Rich.Mean
y_vec[rownames(smallDivIndTrimmed) %in% c("Bacteroidia")] <- 
  y_vec[rownames(smallDivIndTrimmed) %in% c("Bacteroidia")]-40
y_vec[rownames(smallDivIndTrimmed) %in% 
        c("Cryptophyceae", "Phycisphaerae", "Planctomycetacia")] <- 
  y_vec[rownames(smallDivIndTrimmed) %in% 
          c("Cryptophyceae", "Phycisphaerae", "Planctomycetacia")]-10
y_vec[rownames(smallDivIndTrimmed) %in% c("Acidimicrobiia", "Alphaproteobacteria")] <- 
  y_vec[rownames(smallDivIndTrimmed) %in% c("Acidimicrobiia", "Alphaproteobacteria")]+30
y_vec[rownames(smallDivIndTrimmed) %in% c("Oxyphotobacteria")] <- 
  y_vec[rownames(smallDivIndTrimmed) %in% c("Oxyphotobacteria")]+35
y_vec[rownames(smallDivIndTrimmed) %in% c("Gammaproteobacteria")] <- 
  y_vec[rownames(smallDivIndTrimmed) %in% c("Gammaproteobacteria")]+40
text(x_vec, y_vec, labels=rownames(smallDivIndTrimmed),
     cex= 1.5)
arrows(1.72, 90, 1.9, 75, lwd = 2, length = 0.15)
box(lwd = 2.5)
mtext(bquote(bold("S"["T"]~"(Lehman & Tilman, 2000)")), 1, font = 2, line = 5, 
      cex = 2.5, las = 0, adj = 2)
mtext(bquote(bold("Richness")), 2, font = 2, line = 5, 
      cex = 2.5, las = 0)
mtext("Small (0.2-3 μm)", 3, font = 2, line = 2, cex = 2.5)
axis(1, cex.axis=2.5, tck = -0.03, lwd.ticks = 2.5, padj = 1.2, font = 2, 
     at = seq(1, 6, by = 1), labels = seq(1, 6, by = 1))
axis(1, cex.axis=2.5, tck = -0.015, lwd.ticks = 2.5, padj = 1.2, font = 2, 
     at = seq(0.5, 6.5, by = 1), labels = NA)
axis(2, cex.axis=2.5, tck = -0.025, lwd.ticks = 2.5, hadj = 1.2, font = 2, 
     at = seq(0, 400, by = 200), labels = seq(0, 400, by = 200))
axis(2, cex.axis=2.5, tck = -0.015, lwd.ticks = 2.5, hadj = 1.2, font = 2, 
     at = seq(100, 500.8, by = 100), labels = NA)

par(new = "TRUE",plt = c(0.56,0.98,0.1,0.9),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(largeDivIndTrimmed$Stabs, largeDivIndTrimmed$Rich.Mean, 
     cex = largeDivIndTrimmed$Abund/50, xlab = "", ylab = "", lwd = 2, 
     col = largeDivIndTrimmed$Org, xlim = c(0.5, 3.5), xaxt = "n", ylim = c(0, 110),
     yaxt = "n")
x_vec <- largeDivIndTrimmed$Stabs
x_vec[rownames(largeDivIndTrimmed) %in% c("Oxyphotobacteria")] <- 
  x_vec[rownames(largeDivIndTrimmed) %in% c("Oxyphotobacteria")]-0.8
x_vec[rownames(largeDivIndTrimmed) %in% c("Phycisphaerae")] <- 
  x_vec[rownames(largeDivIndTrimmed) %in% c("Phycisphaerae")]-0.4
x_vec[rownames(largeDivIndTrimmed) %in% c("Chrysophyceae")] <- 
  x_vec[rownames(largeDivIndTrimmed) %in% c("Chrysophyceae")]+0.1
x_vec[rownames(largeDivIndTrimmed) %in% 
        c("Verrucomicrobiae", "Gammaproteobacteria", "Deltaproteobacteria")] <- 
  x_vec[rownames(largeDivIndTrimmed) %in% 
          c("Verrucomicrobiae", "Gammaproteobacteria", "Deltaproteobacteria")]+0.42
x_vec[rownames(largeDivIndTrimmed) %in% c("Planctomycetacia")] <- 
  x_vec[rownames(largeDivIndTrimmed) %in% c("Planctomycetacia")]+0.8
y_vec <- largeDivIndTrimmed$Rich.Mean
y_vec[rownames(largeDivIndTrimmed) %in% c("Deltaproteobacteria")] <- 
  y_vec[rownames(largeDivIndTrimmed) %in% c("Deltaproteobacteria")]-4
y_vec[rownames(largeDivIndTrimmed) %in% 
        c("Actinobacteria", "Gammaproteobacteria", "Verrucomicrobiae")] <- 
  y_vec[rownames(largeDivIndTrimmed) %in% 
          c("Actinobacteria", "Gammaproteobacteria", "Verrucomicrobiae")]+5
y_vec[rownames(largeDivIndTrimmed) %in% 
        c("Chrysophyceae")] <- 
  y_vec[rownames(largeDivIndTrimmed) %in% 
          c("Chrysophyceae")]+7
text(x_vec, y_vec, labels=rownames(largeDivIndTrimmed),
     cex= 1.5)
arrows(2.09, 74, 1.86, 74, lwd = 2, length = 0.15)
arrows(2, 92, 1.86, 90, lwd = 2, length = 0.15)
mtext("Large (>3 μm)", 3, font = 2, line = 2, cex = 2.5)
box(lwd = 2.5)
axis(1, cex.axis=2.5, tck = -0.03, lwd.ticks = 2.5, padj = 1.2, font = 2, 
     at = seq(1, 6, by = 1), labels = seq(1, 6, by = 1))
axis(1, cex.axis=2.5, tck = -0.015, lwd.ticks = 2.5, padj = 1.2, font = 2, 
     at = seq(0.5, 6.5, by = 1), labels = NA)
axis(2, cex.axis=2.5, tck = -0.025, lwd.ticks = 2.5, hadj = 1.2, font = 2, 
     at = seq(0, 100, by = 50), labels = seq(0, 100, by = 50))
axis(2, cex.axis=2.5, tck = -0.015, lwd.ticks = 2.5, hadj = 1.2, font = 2, 
     at = seq(25, 125, by = 50), labels = NA)

dev.off()




png("FigBinExtras/ClassDS/diversity.png", width = 1800, height = 945)
plot.new()
par(new = "TRUE",plt = c(0.06,0.48,0.1,0.9),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(smallDivIndTrimmed$Stabs, smallDivIndTrimmed$Shannon.Mean, xlab = "", ylab = "",
     cex = smallDivIndTrimmed$Abund/50, col = smallDivIndTrimmed$Org,
     xlim = c(0.5, 3.5), xaxt = "n", lwd = 2, ylim = c(1, 3.5), yaxt = "n")
x_vec <- smallDivIndTrimmed$Stabs
x_vec[rownames(smallDivIndTrimmed) %in% c("Planctomycetacia", "Verrucomicrobiae")] <- 
  x_vec[rownames(smallDivIndTrimmed) %in% c("Planctomycetacia", "Verrucomicrobiae")]-0.45
x_vec[rownames(smallDivIndTrimmed) %in% c("Acidimicrobiia")] <- 
  x_vec[rownames(smallDivIndTrimmed) %in% c("Acidimicrobiia")]-0.35
x_vec[rownames(smallDivIndTrimmed) %in% c("Actinobacteria", "Gammaproteobacteria")] <- 
  x_vec[rownames(smallDivIndTrimmed) %in% c("Actinobacteria", "Gammaproteobacteria")]-0.15
x_vec[rownames(smallDivIndTrimmed) %in% c("Alphaproteobacteria")] <- 
  x_vec[rownames(smallDivIndTrimmed) %in% c("Alphaproteobacteria")]+0.5
y_vec <- smallDivIndTrimmed$Shannon.Mean
y_vec[rownames(smallDivIndTrimmed) %in% c("Gammaproteobacteria")] <- 
  y_vec[rownames(smallDivIndTrimmed) %in% c("Gammaproteobacteria")]-0.9
y_vec[rownames(smallDivIndTrimmed) %in% c("Bacteroidia")] <- 
  y_vec[rownames(smallDivIndTrimmed) %in% c("Bacteroidia")]-0.6
y_vec[rownames(smallDivIndTrimmed) %in% c("Acidimicrobiia")] <- 
  y_vec[rownames(smallDivIndTrimmed) %in% c("Acidimicrobiia")]-0.1
y_vec[rownames(smallDivIndTrimmed) %in% c("Phycisphaerae", "Cryptophyceae")] <- 
  y_vec[rownames(smallDivIndTrimmed) %in% c("Phycisphaerae", "Cryptophyceae")]+0.06
y_vec[rownames(smallDivIndTrimmed) %in% c("Actinobacteria")] <- 
  y_vec[rownames(smallDivIndTrimmed) %in% c("Actinobacteria")]+0.15
y_vec[rownames(smallDivIndTrimmed) %in% c("Oxyphotobacteria")] <- 
  y_vec[rownames(smallDivIndTrimmed) %in% c("Oxyphotobacteria")]+0.2
text(x_vec, y_vec, labels=rownames(smallDivIndTrimmed),
     cex= 1.5)
box(lwd = 2.5)
arrows(3.5, 2.12, 3.45, 2.82, lwd = 2, length = 0.15)
arrows(3.2, 2.3, 3.2, 2.65, lwd = 2, length = 0.15)
mtext(bquote(bold("S"["T"]~"(Lehman & Tilman, 2000)")), 1, font = 2, line = 5, 
      cex = 2.5, las = 0, adj = 2)
mtext(bquote(bold("Shannon's Diversity Index")), 2, font = 2, line = 5, 
      cex = 2.5, las = 0)
mtext("Small (0.2-3 μm)", 3, font = 2, line = 2, cex = 2.5)
axis(1, cex.axis=2.5, tck = -0.03, lwd.ticks = 2.5, padj = 1.2, font = 2, 
     at = seq(1, 6, by = 1), labels = seq(1, 6, by = 1))
axis(1, cex.axis=2.5, tck = -0.015, lwd.ticks = 2.5, padj = 1.2, font = 2, 
     at = seq(0.5, 6.5, by = 1), labels = NA)
axis(2, cex.axis=2.5, tck = -0.025, lwd.ticks = 2.5, hadj = 1.2, font = 2, 
     at = seq(1, 6, by = 1), labels = seq(1, 6, by = 1))
axis(2, cex.axis=2.5, tck = -0.015, lwd.ticks = 2.5, hadj = 1.2, font = 2, 
     at = seq(0.5, 6.5, by = 1), labels = NA)

par(new = "TRUE",plt = c(0.56,0.98,0.1,0.9),las = 1, cex.axis = 2.5, xpd = FALSE)
plot(largeDivIndTrimmed$Stabs, largeDivIndTrimmed$Shannon.Mean, xlab = "", ylab = "",
     cex = largeDivIndTrimmed$Abund/50, col = largeDivIndTrimmed$Org,
     xlim = c(0.5, 3.5), xaxt = "n", lwd = 2, ylim = c(1, 3.5), yaxt = "n")
x_vec <- largeDivIndTrimmed$Stabs
x_vec[rownames(largeDivIndTrimmed) %in% c("Deltaproteobacteria")] <- 
  x_vec[rownames(largeDivIndTrimmed) %in% c("Deltaproteobacteria")]-0.6
x_vec[rownames(largeDivIndTrimmed) %in% c("Chrysophyceae")] <- 
  x_vec[rownames(largeDivIndTrimmed) %in% c("Chrysophyceae")]+0.15
x_vec[rownames(largeDivIndTrimmed) %in% c("Verrucomicrobiae")] <- 
  x_vec[rownames(largeDivIndTrimmed) %in% c("Verrucomicrobiae")]+0.45
y_vec <- largeDivIndTrimmed$Shannon.Mean
y_vec[rownames(largeDivIndTrimmed) %in% c("Bacillariophyta")] <- 
  y_vec[rownames(largeDivIndTrimmed) %in% c("Bacillariophyta")]-0.3
y_vec[rownames(largeDivIndTrimmed) %in% 
        c("Chrysophyceae", "Phycisphaerae", "Actinobacteria", 
          "Gammaproteobacteria", "Deltaproteobacteria")] <- 
  y_vec[rownames(largeDivIndTrimmed) %in% 
          c("Chrysophyceae", "Phycisphaerae", "Actinobacteria", 
            "Gammaproteobacteria", "Deltaproteobacteria")]+0.14
y_vec[rownames(largeDivIndTrimmed) %in% c("Planctomycetacia")] <- 
  y_vec[rownames(largeDivIndTrimmed) %in% c("Planctomycetacia")]+0.25
text(x_vec, y_vec, labels=rownames(largeDivIndTrimmed),
     cex= 1.5)
mtext("Large (>3 μm)", 3, font = 2, line = 2, cex = 2.5)
box(lwd = 2.5)
arrows(1.16, 2.85, 1.32, 2.75, lwd = 2, length = 0.15)
axis(1, cex.axis=2.5, tck = -0.03, lwd.ticks = 2.5, padj = 1.2, font = 2, 
     at = seq(1, 6, by = 1), labels = seq(1, 6, by = 1))
axis(1, cex.axis=2.5, tck = -0.015, lwd.ticks = 2.5, padj = 1.2, font = 2, 
     at = seq(0.5, 6.5, by = 1), labels = NA)
axis(2, cex.axis=2.5, tck = -0.025, lwd.ticks = 2.5, hadj = 1.2, font = 2, 
     at = seq(1, 6, by = 1), labels = seq(1, 6, by = 1))
axis(2, cex.axis=2.5, tck = -0.015, lwd.ticks = 2.5, hadj = 1.2, font = 2, 
     at = seq(0.5, 6.5, by = 1), labels = NA)

dev.off()