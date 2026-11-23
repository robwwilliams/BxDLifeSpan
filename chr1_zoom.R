#
# chr1_zoom.R - Proximal Chr 1 (0 to 50 Mb) mapped with and without the BXD x BXD F1s: zoomed actuarial maps, profiles,
# and the 1.5-LOD support interval at every T-age
#
# Rob Williams's fork, 23 November 2026. Exploratory. Rob: "What I do not see at the top is our fun and games that include
# the BXD-by-BXD F1 data. That is of greatest interest to me in improving confidence of loci and narrowing 1.5 LOD (or
# –logP) of the support interval. ... the maps need to be compared to excluding the F1s. Need zoom-in versions of the
# proximal 50 Mb of Chr 1."
# No new scans: reads the maps of qtl2_strainmeans.R (strains only, equal weights) and qtl2_strainF1.R (strains and F1
# types, additive and dominance), on chow (CD) and high fat (HF).
#   - Support interval: at each T-age, the run of neighbouring markers around the strongest marker in 0-50 Mb whose LOD is
#     within 1.5 of the peak's. LOD comes back from -log10 P on 1 d.f. A marker outside that run that also comes within
#     1.5 LOD is flagged, since the interval then assumes one locus where there may be two.
#   - Intervals are drawn on the maps where the peak reaches -log10 P 3 (about the suggestive level of the strain-means
#     maps); the table lists every T-age.
# Writes output/maps/Chr1_zoom_maps.svg (six zoomed maps), output/maps/Chr1_zoom_profiles.svg (profiles at the peak
# T-ages, with the F1 types' heterozygous segments) and output/chr1_zoom_intervals.txt
#   Rscript chr1_zoom.R
#

library(svglite)
source("shared.R")
source("actuary_map.R")

zoomMb <- 50
tps <- seq(20, 770, 30)
short <- c("BXD9", "BXD170", "BXD172", "BXD197"); long <- c("BXD65", "BXD73b", "BXD77", "BXD102")
col.strains <- "grey35"; col.F1 <- "#008300"                                      # F1 hybrids green (notes/figure_style.md)

k <- rownames(map)[map[, "Chr"] == "1" & map[, "Mb"] <= zoomMb]
k <- k[order(map[k, "Mb"])]
mb <- as.numeric(map[k, "Mb"])
readMap <- function(f){ x <- read.table(file.path("output", f), sep = "\t"); as.matrix(x[k, paste0("X", tps)]) }
res <- list(CD = list(strains = readMap("qtl2_CD_means_equal_Progressive.txt"), add = readMap("qtl2_CD_meansF1_add_Progressive.txt"), dom = readMap("qtl2_CD_meansF1_dom_Progressive.txt")),
            HF = list(strains = readMap("qtl2_HF_means_equal_Progressive.txt"), add = readMap("qtl2_HF_meansF1_add_Progressive.txt"), dom = readMap("qtl2_HF_meansF1_dom_Progressive.txt")))
nS <- read.table("output/qtl2_means_N.txt", sep = "\t")
nF <- read.table("output/qtl2_meansF1_N.txt", sep = "\t", header = TRUE)
nUnits <- function(d, m){ if(m == "strains") nS[as.character(tps), d] else { n <- nF[nF[, "diet"] == d, ]; (n[, "strains"] + n[, "F1"])[match(tps, n[, "T"])] } }
nF1 <- function(d){ n <- nF[nF[, "diet"] == d, ]; n[match(tps, n[, "T"]), "F1"] }
lod <- function(p) qchisq(pmin(10^-abs(p), 1), 1, lower.tail = FALSE) / (2 * log(10))

### The support interval at each T-age
supportInterval <- function(p){                                                    # p: signed -log10 P over the markers k
  p[is.na(p)] <- 0                                                                 # markers not scanned count as no signal
  L <- lod(p); i <- which.max(abs(p)); cut <- L[i] - 1.5
  lo <- i; while(lo > 1 && L[lo - 1] >= cut) lo <- lo - 1
  hi <- i; while(hi < length(L) && L[hi + 1] >= cut) hi <- hi + 1
  c(peak = mb[i], p = unname(p[i]), lod = unname(L[i]), lo = mb[lo], hi = mb[hi], width = mb[hi] - mb[lo], other = any(L[-(lo:hi)] >= cut))
}
iv <- do.call(rbind, lapply(c("CD", "HF"), function(d) do.call(rbind, lapply(c("strains", "add"), function(m){
  s <- t(apply(res[[d]][[m]], 2, supportInterval))
  data.frame(diet = d, map = c(strains = "strains", add = "strains+F1")[m], T = tps, units = nUnits(d, m), F1 = if(m == "add") nF1(d) else 0, s, row.names = NULL)
}))))
out <- iv; out[, c("peak", "lo", "hi", "width")] <- round(out[, c("peak", "lo", "hi", "width")], 2); out[, c("p", "lod")] <- round(out[, c("p", "lod")], 2)
out[, "other"] <- ifelse(out[, "other"] == 1, "yes", "no")
sink("output/chr1_zoom_intervals.txt", split = TRUE)
cat("Proximal Chr 1 (0 to", zoomMb, "Mb): strongest marker and its 1.5-LOD support interval at each T-age, strain means alone\n")
cat("(qtl2_strainmeans.R, equal weights) against strain and F1 means together (qtl2_strainF1.R, additive). p: signed -log10 P,\n")
cat("positive when the D allele lives longer. lo-hi: the run of markers within 1.5 LOD of the peak (Mb). other: another\n")
cat("marker in 0-50 Mb, outside the run, also within 1.5 LOD. units: strains (+ F1 types) with at least 3 mice alive at T.\n\n")
print(out, row.names = FALSE)
cat("\nWhere both maps reach -log10 P 3 at a T-age: median width of the interval, strains alone against strains + F1s\n")
for(d in c("CD", "HF")){
  a <- iv[iv$diet == d & iv$map == "strains", ]; b <- iv[iv$diet == d & iv$map == "strains+F1", ]
  ok <- abs(a$p) >= 3 & abs(b$p) >= 3
  cat(sprintf("  %s: %d T-ages (%s); %.1f Mb against %.1f Mb\n", d, sum(ok), paste0("T", tps[ok], collapse = " "), median(a$width[ok]), median(b$width[ok])))
}
cat("\nDominance (qtl2_strainF1.R): the largest -log10 P in 0-50 Mb at any T-age\n")
for(d in c("CD", "HF")){
  x <- res[[d]]$dom; w <- which(abs(x) == max(abs(x), na.rm = TRUE), arr.ind = TRUE)[1, ]
  cat(sprintf("  %s: %+.2f at %.1f Mb, T%d (positive: BD lives longer than the mean of BB and DD)\n", d, x[w[1], w[2]], mb[w[1]], tps[w[2]]))
}
sink()

### Six zoomed actuarial maps, with the support interval of each T-age as a black bar across its row
map0 <- map[k, ]; map0[, "Mb"] <- 1e6 * map0[, "Mb"]
lay <- mapLayout(map0, chrs = "1")
lay$length["1"] <- zoomMb * 1e6; lay$l.x <- lay$start["1"] + zoomMb * 1e6 + 1
panels <- list(
  list(d = "CD", m = "strains", main = "Chow: strain means alone"),
  list(d = "CD", m = "add", main = "Chow: strain and F1 means, additive effect"),
  list(d = "CD", m = "dom", main = "Chow: strain and F1 means, dominance deviation"),
  list(d = "HF", m = "strains", main = "High fat: strain means alone"),
  list(d = "HF", m = "add", main = "High fat: strain and F1 means, additive effect"),
  list(d = "HF", m = "dom", main = "High fat: strain and F1 means, dominance deviation"))
svglite("output/maps/Chr1_zoom_maps.svg", width = 14, height = 5.5 * length(panels), bg = "white", system_fonts = list(sans = "Arial"), fix_text_size = FALSE)
par(mfrow = c(length(panels), 1), mai = map.mai, omi = c(0, 0, 0, 0), cex = 1)
for(i in seq_along(panels)){
  p <- panels[[i]]; x <- res[[p$d]][[p$m]]; rownames(x) <- k
  signs <- if(p$m == "dom") c("*BD* > midpoint", "*BD* < midpoint") else c("*D*: longer life", "*B*: longer life")
  actuaryMap(byTPS(x), map0, lay, N = nUnits(p$d, if(p$m == "strains") "strains" else "add")[match(plotTPS, tps)],
             N.title = if(p$m == "strains") "Strains" else "Strains + F1s", main = p$main, sign.labels = signs, chr.numbers = FALSE, letter = letters[i])
  x0 <- lay$start["1"]
  axis(1, at = x0 + seq(0, zoomMb, 5) * 1e6, labels = seq(0, zoomMb, 5), cex.axis = map.cex, tcl = -0.4, mgp = c(3, 0.7, 0))
  mtext("Chr 1 position (Mb)", side = 1, line = 2.6, cex = map.cex.title)
  if(p$m != "dom"){
    s <- iv[iv$diet == p$d & iv$map == c(strains = "strains", add = "strains+F1")[p$m], ]
    for(j in which(abs(s$p) >= 3 & s$T %in% plotTPS)){
      r <- match(s$T[j], plotTPS)
      segments(x0 + s$lo[j] * 1e6, r + 0.5, x0 + s$hi[j] * 1e6, r + 0.5, lwd = 3, col = "black", lend = 1)
      points(x0 + s$peak[j] * 1e6, r + 0.5, pch = 21, bg = "white", cex = 1.1)
    }
  }
}
invisible(dev.off())

### Profiles at the peak T-ages, strains alone against strains + F1s, over the F1 types' heterozygous segments
f1types <- c("BXD170x172", "BXD172x197", "BXD197x170", "BXD197x9", "BXD9x170", "BXD9x172",
             "BXD102x73b", "BXD102x77", "BXD65x102", "BXD65x77", "BXD73bx65", "BXD73bx77", "BXD77x65")
parentsOf <- function(t) c(sub("x.*", "", t), paste0("BXD", sub(".*x", "", t)))
hetSeg <- lapply(f1types, function(t){ g <- geno[k, parentsOf(t)]; ok <- !is.na(g[, 1]) & !is.na(g[, 2]) & g[, 1] %in% c("B", "D") & g[, 2] %in% c("B", "D"); ok & g[, 1] != g[, 2] })
names(hetSeg) <- f1types
views <- list(list(d = "CD", T = 770), list(d = "HF", T = 20), list(d = "HF", T = 260), list(d = "HF", T = 500))
svglite("output/maps/Chr1_zoom_profiles.svg", width = 7, height = 2.1 * length(views) + 2.6, bg = "white", system_fonts = list(sans = "Arial"), fix_text_size = FALSE)
layout(matrix(1:(length(views) + 1)), heights = c(rep(2.1, length(views)), 2.6))
par(mai = c(0.15, 0.75, 0.25, 0.45), cex = 0.75, mgp = c(2.2, 0.6, 0), tcl = -0.3)
ymax <- ceiling(max(sapply(views, function(v) max(abs(c(res[[v$d]]$strains[, paste0("X", v$T)], res[[v$d]]$add[, paste0("X", v$T)]))))))
for(i in seq_along(views)){
  v <- views[[i]]; cl <- paste0("X", v$T)
  plot(NA, xlim = c(0, zoomMb), ylim = c(-1.2, ymax), xaxs = "i", xaxt = "n", yaxt = "n", bty = "l", xlab = "", ylab = "−log10 P")
  axis(1, at = seq(0, zoomMb, 5), labels = FALSE); axis(2, at = 0:ymax, las = 1)
  mtext(sprintf("%s, T%d", c(CD = "Chow", HF = "High fat")[v$d], v$T), side = 3, line = 0.1, adj = 0, font = 2, cex = 0.8)
  mtext(letters[i], side = 3, line = 0.1, at = grconvertX(0.05, "inches", "user"), adj = 0, font = 2, cex = 1, xpd = NA)
  abline(h = 4.4, lty = 2, col = "grey55")
  text(zoomMb, 4.4, "strain means, genome-wide 5%", adj = c(1, 1.4), cex = 0.8, col = "grey40")
  for(m in c("strains", "add")){
    y <- abs(res[[v$d]][[m]][, cl]); co <- if(m == "strains") col.strains else col.F1
    lines(mb, y, col = co, lwd = 1.3)
    s <- iv[iv$diet == v$d & iv$map == c(strains = "strains", add = "strains+F1")[m] & iv$T == v$T, ]
    yb <- if(m == "strains") -0.4 else -0.9
    segments(s$lo, yb, s$hi, yb, col = co, lwd = 4, lend = 1)
    lab <- sprintf("%.1f–%.1f Mb", s$lo, s$hi)                                   # beside the bar, on whichever side has room
    if(s$hi + 0.6 + strwidth(lab, cex = 0.75) < zoomMb) text(s$hi + 0.6, yb, lab, adj = c(0, 0.5), cex = 0.75, col = co) else text(s$lo - 0.6, yb, lab, adj = c(1, 0.5), cex = 0.75, col = co)
  }
  n <- nF[nF$diet == v$d & nF$T == v$T, ]
  legend("topleft", bty = "n", cex = 0.8, lwd = 1.3, col = c(col.strains, col.F1), inset = c(0.01, 0),
         legend = c(sprintf("Strains alone (%d)", n$strains), sprintf("Strains + F1 types (%d + %d)", n$strains, n$F1)))
}
par(mai = c(0.55, 0.75, 0.25, 0.45))
plot(NA, xlim = c(0, zoomMb), ylim = c(length(f1types) + 0.5, 0.5), xaxs = "i", yaxt = "n", bty = "l", xlab = "Chr 1 position (Mb)", ylab = "")
mtext("F1 types: where the parents differ", side = 3, line = 0.1, adj = 0, font = 2, cex = 0.8)
mtext(letters[length(views) + 1], side = 3, line = 0.1, at = grconvertX(0.05, "inches", "user"), adj = 0, font = 2, cex = 1, xpd = NA)
axis(2, at = seq_along(f1types), labels = f1types, las = 1, cex.axis = 0.6, tick = FALSE, line = -0.6)
for(j in seq_along(f1types)){
  h <- hetSeg[[j]]; r <- rle(h); e <- cumsum(r$lengths); s0 <- e - r$lengths + 1
  for(q in which(r$values)) rect(mb[s0[q]], j - 0.35, mb[e[q]], j + 0.35, col = col.F1, border = NA)
}
abline(h = 6.5, col = "grey60", lwd = 0.5)
for(g in list(c(3.5, "short × short"), c(10, "long × long"))) text(grconvertX(grconvertX(zoomMb, "user", "inches") + 0.15, "inches", "user"), as.numeric(g[1]), g[2], srt = -90, cex = 0.7, col = "grey30", xpd = NA)
invisible(dev.off())
cat("Wrote output/maps/Chr1_zoom_maps.svg, output/maps/Chr1_zoom_profiles.svg and output/chr1_zoom_intervals.txt\n")
