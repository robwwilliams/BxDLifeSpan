#
# vita_zoom_figures.R - A tight zoom of Vita1a, and overviews of proximal Chr 1 and distal Chr 2: BXD and UM-HET3 maps
# beside the founder-specific variants
#
# Rob Williams's fork, 24 November 2026. Exploratory. Rob: "We definitely want a very tight zoom for the bioRxiv paper and
# the next grant application of this region with the LTR prominently displayed, as well as any other strain specific
# oddities that might help us make sense of what is going on across the proximal 40 Mb on Chr 1 and much of Chr 2."
#   - Vita1a zoom, 4.75 to 5.60 Mb (mm10):
#       a  BXD chow, T770: -log10 P of strain means alone and with the F1s, their 1.5-LOD intervals, and each map's 5%
#          threshold; the block of identical genotypes (4.858-5.316 Mb) shaded;
#       b  the same map with the F1s at every T-age (signed: warm D, cool B);
#       c  the BXD strains whose genotypes change in the window, and the F1 parents, as B and D segments;
#       d  RefSeq genes, and the 11 kb LTR in Atp6v1h (only the D haplotype carries it, to Rob's knowledge; its position is
#          a placeholder until David Ashbrook and Flavia Villani's long reads place it);
#       e  founder-specific variants (Supplementary Table 17 of Arends et al. 2026) in 5-kb bins, by pattern, with the
#          coding sites marked.
#     UM-HET3 has no markers inside the window (its nearest are at 4.3 and 6.2 Mb); its paternal and maternal values there
#     are given in the margin.
#   - Overviews, Chr 1 0-45 Mb and Chr 2 100-182 Mb: the strongest signed -log10 P over T-ages at each marker, for BXD (strain
#     and F1 means, chow and high fat) and UM-HET3 females (paternal and maternal); Rob's Vita loci (from his spreadsheet
#     QTL_Lifespan_Master_16Oct2026_rob.xlsx) and Mrps5 (Houtkooper et al. 2013); and founder-specific variants in 250-kb
#     bins.
# Inputs: this fork's output/ (maps); ../bxd_mortality_rates_2026/output/ (founder variants) and its UMHET3/vita1a_refseq_mm10.txt
# Writes output/maps/Vita1a_zoom.svg, output/maps/Chr1_Chr2_overview.svg
#   Rscript vita_zoom_figures.R
#

library(svglite)
source("shared.R")
main <- "../bxd_mortality_rates_2026"
cB <- "#7d4598"; cD <- "#ed8f22"; cF1 <- "#008300"; cCD <- "#a86a28"; cHF <- "#2a78d6"   # notes/figure_style.md
short <- c("D-private" = "D only", "C3H-private" = "C3H only", "B6-private" = "B6 only", "BALB/c-private" = "BALB/c only", "H = D, against C = B" = "H = D | C = B",
           "B = H, against C = D" = "B = H | C = D", "B = D, against C = H" = "B = D | C = H", "three or more alleles" = "3 or more alleles")
pcol <- c("D-private" = cD, "C3H-private" = "black", "B6-private" = cB, "BALB/c-private" = "grey55",
          "H = D, against C = B" = "#f6c48a", "B = H, against C = D" = "#b996c8", "B = D, against C = H" = "#6e6e6e", "three or more alleles" = "#c42e2f")
rd <- function(f) read.table(f, sep = "\t", check.names = FALSE)
lp1 <- function(x, cols) apply(as.matrix(x[, cols]), 1, function(v){ v[is.na(v)] <- 0; v[which.max(abs(v))] })   # the strongest signed value over T-ages
vita <- data.frame(name = c("Vita1a", "Vita1b", "Vita1x", "Vita1y", "Vita2b"), chr = c("1", "1", "1", "1", "2"),
                   from = c(4.858, 22.361, 25.535, 32.325, 153.171), to = c(10.499, 25.257, 31.408, 43.007, 158.948),
                   core.from = c(5.019, 22.852, NA, NA, 154.141), core.to = c(5.522, 23.515, NA, NA, 157.312))   # Rob's spreadsheet
ltr <- c(from = 5.110, to = 5.121)                                                    # PLACEHOLDER: 11 kb inside Atp6v1h, position not yet known
bxdT <- as.character(seq(20, 770, 30))

### Vita1a zoom
X <- c(4.75, 5.60); block <- c(4.858261, 5.316455)
sF <- rd("output/qtl2_CD_meansF1_add_Progressive.txt"); sS <- rd("output/qtl2_CD_means_equal_Progressive.txt")
k <- rownames(sF)[sF$Chr == "1" & sF$Mb >= X[1] - 0.2 & sF$Mb <= X[2] + 0.2]; k <- k[order(sF[k, "Mb"])]
iv <- read.table("output/chr1_zoom_intervals.txt", skip = 5, header = TRUE, nrows = 104, fill = TRUE)
ivT <- function(m) iv[iv$diet == "CD" & iv$map == m & iv$T == 770, ]
pm <- read.delim("output/qtl2_permutation_maxima.txt"); pf <- read.delim("output/qtl2_strainF1_perm_maxima.txt", check.names = FALSE)
th <- c(quantile(pm[pm$diet == "CD", "T770"], 0.95), quantile(pf[pf$diet == "CD" & pf$map == "add" & pf$perm > 0, "T770"], 0.95))
gm <- read.delim(file.path(main, "UMHET3/vita1a_refseq_mm10.txt"), comment.char = "#")
sites <- read.delim(file.path(main, "output/Vita1a_founder_variants_sites.txt"))
sites <- sites[sites$pos / 1e6 >= X[1] & sites$pos / 1e6 <= X[2], ]
sites$label[sites$label %in% c("three alleles", "three or more alleles")] <- "three or more alleles"
sites$label <- sub("H-private \\(C3H\\)", "C3H-private", sub("B-private \\(B6\\)", "B6-private", sub("C-private \\(BALB/c\\)", "BALB/c-private", sites$label)))
cds <- read.delim(file.path(main, "UMHET3/vita1a_refseq_cds_mm10.txt"), comment.char = "#")
sites$coding <- sapply(sites$pos, function(p) any(p >= cds$cds_from & p <= cds$cds_to))
showStrains <- c("BXD65", "BXD65b", "BXD102", "BXD43", "BXD79", "BXD52", "BXD11", "BXD77", "BXD73b", "BXD9", "BXD170", "BXD172", "BXD197")
gk <- rownames(map)[map$Chr == "1" & map$Mb >= X[1] - 0.3 & map$Mb <= X[2] + 0.3]; gk <- gk[order(map[gk, "Mb"])]
um <- list(paternal = rd("output/umhet3_paternal_female_Progressive.txt"), maternal = rd("output/umhet3_maternal_female_Progressive.txt"))
umAt <- function(m, mb){ x <- um[[m]]; i <- which(x$Chr == "1" & abs(x$Mb - mb) < 0.05)[1]; v <- unlist(x[i, as.character(seq(20, 860, 30))]); median(v) }

svglite("output/maps/Vita1a_zoom.svg", width = 7.2, height = 9.6, bg = "white", system_fonts = list(sans = "Arial"), fix_text_size = FALSE)
layout(matrix(1:5), heights = c(2.0, 1.9, 1.7, 1.5, 2.5))
par(oma = c(3.6, 0, 0.3, 0), cex = 0.7, mgp = c(2.1, 0.5, 0), tcl = -0.25, las = 1)
xax <- function(lab = FALSE){ axis(1, at = seq(4.75, 5.6, 0.05), labels = FALSE, tcl = -0.15); axis(1, at = seq(4.8, 5.6, 0.1), labels = lab) }
shade <- function(){ rect(block[1], par("usr")[3], block[2], par("usr")[4], col = "#f3efe6", border = NA) }
letter <- function(l, txt){ mtext(l, side = 3, line = 0.2, at = grconvertX(0.06, "inches", "user"), adj = 0, font = 2, cex = 0.95, xpd = NA); mtext(txt, side = 3, line = 0.2, adj = 0, cex = 0.8, font = 2) }
# a: profiles at T770
par(mai = c(0.15, 0.75, 0.3, 0.2))
plot(NA, xlim = X, ylim = c(-1.2, 7), xaxs = "i", xaxt = "n", bty = "l", xlab = "", ylab = "−log10 P"); shade()
lines(sS[k, "Mb"], abs(sS[k, "770"]), col = "grey35", lwd = 1.3, type = "s"); lines(sF[k, "Mb"], abs(sF[k, "770"]), col = cF1, lwd = 1.3, type = "s")
abline(h = th, lty = 2, lwd = 0.8, col = c("grey35", cF1))
for(m in c("strains", "strains+F1")){ s <- ivT(m); y <- if(m == "strains") -0.4 else -0.9; segments(s$lo, y, s$hi, y, lwd = 4, lend = 1, col = if(m == "strains") "grey35" else cF1) }
legend(X[2] - 0.005, 3.7, bty = "n", cex = 0.78, xjust = 1, lwd = c(1.3, 1.3, 0.8, 4), lty = c(1, 1, 2, 1), col = c("grey35", cF1, "grey50", "grey60"), yjust = 0.5,
       legend = c("BXD strain means (46)", "strains + F1 types (46 + 13)", "genome-wide 5%, T770", "1.5-LOD intervals"))
text(mean(block), 6.85, "block of identical BXD genotypes", cex = 0.75, col = "grey40")
letter("a", "BXD females on chow, T770: the D allele lives longer"); xax()
# b: the strain + F1 map at every T-age
par(mai = c(0.15, 0.75, 0.3, 0.2))
plot(NA, xlim = X, ylim = c(0, 26), xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "T-age (days)")
mb <- sF[k, "Mb"]; xs <- c(X[1] - 1, (mb[-1] + mb[-length(mb)]) / 2); xe <- c((mb[-1] + mb[-length(mb)]) / 2, X[2] + 1)
br <- c(0, 1, 2, 3, 4, 5, 100); cw <- c("#fceea7", "#fdb863", "#e66101", "#d01c8b", "#c42e2f"); cc <- c("#bfe6d8", "#66c2a4", "#2b8cbe", "#253494", "#081d58")
for(j in seq_along(bxdT)){ v <- sF[k, bxdT[j]]; v[is.na(v)] <- 0; bin <- findInterval(abs(v), br, left.open = TRUE)
  for(i in which(abs(v) > 1)) rect(max(xs[i], X[1]), j - 1, min(xe[i], X[2]), j, col = if(v[i] > 0) cw[bin[i] - 1] else cc[bin[i] - 1], border = NA) }
abline(h = 0:26, col = "grey85", lwd = 0.3); box(lwd = 0.6)
axis(2, at = seq(1, 26, 4) - 0.5, labels = bxdT[seq(1, 26, 4)])
rect(block[1], 0, block[2], 26, border = "grey30", lwd = 0.6, lty = 3)
letter("b", "BXD strain + F1 means on chow at every T-age (warm: D longer; white: −log10 P below 1)"); xax()
# c: crossovers
par(mai = c(0.15, 0.75, 0.3, 0.2))
plot(NA, xlim = X, ylim = c(length(showStrains) + 0.5, 0.5), xaxs = "i", xaxt = "n", yaxt = "n", bty = "l", xlab = "", ylab = ""); shade()
axis(2, at = seq_along(showStrains), labels = showStrains, cex.axis = 0.75, tick = FALSE, line = -0.5)
gmb <- map[gk, "Mb"]; gx <- c(X[1] - 1, (gmb[-1] + gmb[-length(gmb)]) / 2, X[2] + 1)
for(j in seq_along(showStrains)){ g <- as.character(geno[gk, showStrains[j]])
  for(i in seq_along(gk)) if(!is.na(g[i]) && g[i] %in% c("B", "D")) rect(max(gx[i], X[1]), j - 0.32, min(gx[i + 1], X[2]), j + 0.32, col = if(g[i] == "B") cB else cD, border = NA) }
legend(X[2], 0.5, bty = "n", cex = 0.8, fill = c(cB, cD), border = NA, legend = c("B (C57BL/6J)", "D (DBA/2J)"), horiz = TRUE, xjust = 1, yjust = 0, xpd = NA)
letter("c", "BXD strains with crossovers here, and the F1 parents"); xax()
# d: genes and the LTR
par(mai = c(0.15, 0.75, 0.3, 0.2))
tx <- gm[!duplicated(paste(gm$gene, gm$tx_from, gm$tx_to)), ]; tx <- tx[tx$tx_to / 1e6 > X[1] & tx$tx_from / 1e6 < X[2], ]
tx <- tx[order(tx$tx_from), ]
plot(NA, xlim = X, ylim = c(nrow(tx) + 1.2, 0.2), xaxs = "i", xaxt = "n", yaxt = "n", bty = "l", xlab = "", ylab = ""); shade()
for(j in seq_len(nrow(tx))){ r <- tx[j, ]; segments(r$tx_from / 1e6, j, r$tx_to / 1e6, j, col = "grey30", lwd = 0.8)
  es <- as.numeric(strsplit(r$exon_from, ",")[[1]]) / 1e6; ee <- as.numeric(strsplit(r$exon_to, ",")[[1]]) / 1e6
  rect(es, j - 0.28, ee, j + 0.28, col = "grey20", border = NA)
  x0 <- min(max(X[1] + 0.002, r$tx_from / 1e6), X[2] - 0.17)
  text(x0, j - 0.36, r$gene, adj = c(0, 0), cex = 0.62, font = 3, col = "grey15")
  text(x0 + strwidth(r$gene, cex = 0.62, font = 3) + 0.004, j - 0.36, sprintf("%s (%s)", r$transcript, r$strand), adj = c(0, 0), cex = 0.55, col = "grey35") }
rect(ltr[1], 0.35, ltr[2], nrow(tx) + 1.1, col = adjustcolor("#c42e2f", 0.18), border = "#c42e2f", lwd = 1.2, lty = 2)
text(ltr[2] + 0.004, nrow(tx) + 0.75, "11 kb LTR, D haplotype only (position to come from long reads)", adj = c(0, 0.5), cex = 0.68, col = "#c42e2f", font = 2)
letter("d", "Genes (RefSeq, mm10) and the Atp6v1h LTR"); xax()
# e: founder variants
par(mai = c(0.15, 0.75, 0.3, 0.2))
pats <- names(pcol); bw <- 0.005; bins <- seq(X[1], X[2], bw)
plot(NA, xlim = X, ylim = c(length(pats) + 1.5, 0.5), xaxs = "i", xaxt = "n", yaxt = "n", bty = "l", xlab = "", ylab = ""); shade()
axis(2, at = seq_along(pats), labels = short[pats], cex.axis = 0.7, tick = FALSE, line = -0.5)
mx <- 0
for(j in seq_along(pats)){ h <- hist(sites$pos[sites$label == pats[j]] / 1e6, breaks = bins, plot = FALSE)$counts; mx <- max(mx, h)
  for(i in which(h > 0)) rect(bins[i], j - 0.4, bins[i + 1], j + 0.4, col = adjustcolor(pcol[j], min(1, 0.15 + h[i] / 25)), border = NA)
  cs <- sites[sites$label == pats[j] & sites$coding, ]; if(nrow(cs)) points(cs$pos / 1e6, rep(j, nrow(cs)), pch = 25, bg = "white", col = "black", cex = 0.9) }
text(X[1] + 0.005, length(pats) + 1.05, "Shade: variants per 5 kb (darker = more). Triangles: coding sites. Founders: C BALB/cByJ, B C57BL/6NJ, H C3H/HeJ, D DBA/2J.", adj = 0, cex = 0.6, col = "grey30")
letter("e", "Founder-specific variants (SNPs and short indels; Arends et al. 2026, Supplementary Table 17)"); xax(TRUE)
mtext("Chr 1 position (Mb, mm10)", side = 1, line = 1.6, cex = 0.75)
mtext(sprintf("UM-HET3 has no markers in this window. At its nearest (4.3 and 6.2 Mb), females' median signed −log10 P over T20–T860 is %+.1f and %+.1f paternal (H longer),",
              umAt("paternal", 4.3), umAt("paternal", 6.2)), side = 1, outer = TRUE, line = 1.9, cex = 0.6, col = "grey25", adj = 0, at = 0.02)
mtext(sprintf("and %+.1f and %+.1f maternal (B against C): a flat maternal map.", umAt("maternal", 4.3), umAt("maternal", 6.2)), side = 1, outer = TRUE, line = 2.7, cex = 0.6, col = "grey25", adj = 0, at = 0.02)
invisible(dev.off())

### Overviews: proximal Chr 1 and distal Chr 2
fp <- list("1" = read.delim(file.path(main, "output/founder_patterns_chr1_0-45Mb.txt"), check.names = FALSE), "2" = read.delim(file.path(main, "output/founder_patterns_chr2_100-182Mb.txt"), check.names = FALSE))
bx <- list(CD = rd("output/qtl2_CD_meansF1_add_Progressive.txt"), HF = rd("output/qtl2_HF_meansF1_add_Progressive.txt"))
um2 <- list(paternal = um$paternal, maternal = um$maternal)
svglite("output/maps/Chr1_Chr2_overview.svg", width = 7.2, height = 8.4, bg = "white", system_fonts = list(sans = "Arial"), fix_text_size = FALSE)
layout(matrix(1:6, ncol = 2), heights = c(2.2, 2.2, 3.0))
par(oma = c(2.2, 0, 0.3, 0.1), cex = 0.62, mgp = c(2, 0.5, 0), tcl = -0.25, las = 1)
for(ch in c("1", "2")){
  X <- if(ch == "1") c(0, 45) else c(100, 182); L <- if(ch == "1") c("a", "b", "c") else c("d", "e", "f")
  vl <- vita[vita$chr == ch, ]
  sh <- function(){ for(i in seq_len(nrow(vl))){ rect(vl$from[i], par("usr")[3], vl$to[i], par("usr")[4], col = "#f3efe6", border = NA)
    if(!is.na(vl$core.from[i])) rect(vl$core.from[i], par("usr")[3], vl$core.to[i], par("usr")[4], col = "#e8dcc3", border = NA) } }
  # BXD
  par(mai = c(0.15, 0.85, 0.35, 0.1))
  plot(NA, xlim = X, ylim = c(-6.5, 6.5), xaxs = "i", xaxt = "n", bty = "l", xlab = "", ylab = "signed −log10 P"); sh(); abline(h = 0, col = "grey70", lwd = 0.5)
  for(d in c("CD", "HF")){ x <- bx[[d]]; i <- which(x$Chr == ch & x$Mb >= X[1] & x$Mb <= X[2]); i <- i[order(x$Mb[i])]; lines(x$Mb[i], lp1(x[i, ], bxdT), col = if(d == "CD") cCD else cHF, lwd = 0.9) }
  text(vl$from + (vl$to - vl$from) / 2, 6.1, vl$name, cex = 0.8, font = 3)
  if(ch == "2"){ abline(v = 127.59, lty = 3, col = "grey30"); text(127.59, 5.2, "Mrps5", cex = 0.75, font = 3, adj = c(-0.1, 0.5)); rect(124, 6.2, 129, 6.5, col = "grey40", border = NA) }
  legend("bottomleft", bty = "n", cex = 0.8, lwd = 0.9, col = c(cCD, cHF), legend = c("BXD + F1 means, chow", "high fat"), inset = c(0, 0))
  mtext(L[1], side = 3, line = 0.3, at = grconvertX(0.01, "nfc", "user"), adj = 0, font = 2, cex = 0.95, xpd = NA)
  mtext(sprintf("Chr %s: BXD + F1 means (positive: D longer)", ch), side = 3, line = 0.3, adj = 0, cex = 0.75, font = 2)
  axis(1, labels = FALSE)
  # UM-HET3
  plot(NA, xlim = X, ylim = c(-9.5, 9.5), xaxs = "i", xaxt = "n", bty = "l", xlab = "", ylab = "signed −log10 P"); sh(); abline(h = 0, col = "grey70", lwd = 0.5)
  for(m in c("paternal", "maternal")){ x <- um2[[m]]; i <- which(x$Chr == ch & x$Mb >= X[1] & x$Mb <= X[2]); i <- i[order(x$Mb[i])]
    lines(x$Mb[i], lp1(x[i, ], as.character(seq(20, 1100, 30))), col = if(m == "paternal") cD else cB, lwd = 0.9, type = "o", pch = 16, cex = 0.35) }
  legend("bottomleft", bty = "n", cex = 0.8, lwd = 0.9, col = c(cD, cB), legend = c("UM-HET3 females, paternal D − H", "maternal B − C"))
  mtext(L[2], side = 3, line = 0.3, at = grconvertX(0.01, "nfc", "user"), adj = 0, font = 2, cex = 0.95, xpd = NA)
  mtext(sprintf("Chr %s: UM-HET3 females", ch), side = 3, line = 0.3, adj = 0, cex = 0.75, font = 2)
  axis(1, labels = FALSE)
  # founder variants
  par(mai = c(0.15, 0.85, 0.35, 0.1))
  f <- fp[[ch]]; pats <- names(pcol)
  plot(NA, xlim = X, ylim = c(length(pats) + 0.5, 0.5), xaxs = "i", xaxt = "n", yaxt = "n", bty = "l", xlab = "", ylab = ""); sh()
  axis(2, at = seq_along(pats), labels = short[pats], cex.axis = 0.72, tick = FALSE, line = -0.6)
  for(j in seq_along(pats)){ h <- f[[pats[j]]]; q <- quantile(h[h > 0], 0.95)
    for(i in which(h > 0)) rect(f$bin.from.Mb[i], j - 0.4, f$bin.from.Mb[i] + 0.25, j + 0.4, col = adjustcolor(pcol[j], min(1, 0.1 + 0.9 * h[i] / q)), border = NA) }
  axis(1); mtext(sprintf("Chr %s position (Mb, mm10)", ch), side = 1, line = 1.6, cex = 0.7)
  mtext(L[3], side = 3, line = 0.3, at = grconvertX(0.01, "nfc", "user"), adj = 0, font = 2, cex = 0.95, xpd = NA)
  mtext("Founder-specific variants per 250 kb", side = 3, line = 0.3, adj = 0, cex = 0.75, font = 2)
}
invisible(dev.off())
cat("Wrote output/maps/Vita1a_zoom.svg and output/maps/Chr1_Chr2_overview.svg\n")
