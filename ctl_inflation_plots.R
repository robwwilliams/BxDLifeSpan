#
# ctl_inflation_plots.R - Q-Q (inflation) plots of the CTL test under three assumptions, in maps with genotypes scrambled among strains
# (no CTL can exist) and in the real data
#
# Rob Williams's fork, 26 November 2026. Exploratory, asked for by Danny Arends and Rob ("inflation plots on the CTL ... the mice being
# the N, the strains being the N, and then the standard"). Written by Claude (Opus 5.5). BXD mice, Danny's design (weight at T).
# Run from the fork's root after ctl_inputs.R: Rscript ctl_inflation_plots.R
#
# The three tests, at every marker and T-age:
#   mice as n       Danny's CTL.R: weighted Spearman r of the mice in B and in D, Fisher z weighted n - 1, n = mice
#   strains as n    the same formula on strain means (strains with at least 2 mice), n = strains, weights n - 1
#   standard        strain means, Fisher z weighted (n - 3) / 1.06 (the variance of a Spearman z, Fieller et al. 1957), and at least
#                   20 strains in each group
# A valid test lies on the diagonal in the scrambled maps: observed -log10 P equals expected. lambda = median chi-square / 0.455.
# Writes output/maps/CTL_inflation.svg (and .png) and output/ctl_inflation.txt
#

source("ctl_functions.R"); library(svglite)
DIR <- "output/ctl_scan"
N.SCRAMBLED <- 10                                                          # scrambled maps, seeds 20261124 + 1 to 10 (as in ctl_scan.R)
SEED <- 20261124
TEST.COLOR <- c(`mice as n` = "black", `strains as n` = "#d55e00", standard = "#009e73")    # Okabe-Ito
DIET.NAME <- c(CD = "Chow", HF = "High fat")

inp <- readInputs(DIR)
Tages <- as.numeric(sub("^w", "", grep("^w[0-9]+$", colnames(inp$mice), value = TRUE)))
scans <- buildScans(inp$mice, inp$strains, c(actuarial = NA), Tages, 2)

# All three tests of one map (genotype order o), for one diet: a list of vectors of -log10 P (distinct patterns x scans)
threeTests <- function(o, diet){
  G <- inp$G[o, , drop = FALSE]; diet.scans <- scans[sapply(scans, `[[`, "diet") == diet]
  strain.maps <- scanAllStrains(G, diet.scans)
  list(`mice as n` = as.vector(abs(scanAll(G, diet.scans))), `strains as n` = as.vector(abs(strain.maps$coded)), standard = as.vector(abs(strain.maps$standard)))
}

real <- lapply(c(CD = "CD", HF = "HF"), function(diet) threeTests(seq_along(inp$strains), diet))
scrambled <- lapply(c(CD = "CD", HF = "HF"), function(diet){
  per.map <- lapply(1:N.SCRAMBLED, function(i) threeTests(permuteStrains(inp$stratum, SEED + i), diet))
  lapply(names(TEST.COLOR), function(t) unlist(lapply(per.map, `[[`, t))) |> setNames(names(TEST.COLOR)) })

lambda <- function(log10P){ p <- 10^-log10P[!is.na(log10P)]; median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1) }
qqPoints <- function(log10P){                                              # observed against expected -log10 P, every test
  obs <- sort(log10P[!is.na(log10P)], decreasing = TRUE); list(x = -log10(seq_along(obs) / (length(obs) + 1)), y = obs) }
thinToDrawing <- function(p, xlim, ylim, cells = 2000){                    # one point per cell of a cells x cells grid over the panel
  cell <- paste(round(p$x / diff(xlim) * cells), round(p$y / diff(ylim) * cells))   # (0.1 pt here): closer points overlap on the page,
  keep <- !duplicated(cell); list(x = p$x[keep], y = p$y[keep]) }                   # so the drawing is the same and the SVG stays small

# ---- The figure: rows scrambled (no CTL) and real; columns chow and high fat ----------------------------------------------------------
dir.create("output/maps", showWarnings = FALSE)
svglite("output/maps/CTL_inflation.svg", width = 7.2, height = 6.4, bg = "white", system_fonts = list(sans = "Arial"), fix_text_size = FALSE)
layout(matrix(1:4, 2, byrow = TRUE)); par(mai = c(0.55, 0.6, 0.35, 0.15), mgp = c(2, 0.45, 0), tcl = -0.2, las = 1, cex = 0.75)
panel <- 0
for(which.maps in c("scrambled", "real")) for(diet in c("CD", "HF")){
  sets <- if(which.maps == "scrambled") scrambled[[diet]] else real[[diet]]
  pts <- lapply(sets, qqPoints); top <- max(sapply(pts, function(p) max(p$y)), sapply(pts, function(p) max(p$x)))
  xlim <- c(0, max(sapply(pts, function(p) max(p$x)))); ylim <- c(0, top)
  plot(NA, xlim = xlim, ylim = ylim, xlab = "Expected -log10 P", ylab = "Observed -log10 P", bty = "l")
  abline(0, 1, col = "grey60", lty = 2)
  for(t in names(sets)){ p <- thinToDrawing(pts[[t]], xlim, ylim); points(p$x, p$y, col = TEST.COLOR[t], pch = 16, cex = 0.25) }
  legend("topleft", sprintf("%s (lambda %.2f)", names(sets), sapply(sets, lambda)), col = TEST.COLOR, pch = 16, bty = "n", cex = 0.85)
  panel <- panel + 1; mtext(letters[panel], 3, line = 0.5, at = grconvertX(0.01, "nfc", "user"), adj = 0, font = 2, cex = 1.15)
  mtext(sprintf("%s, %s", DIET.NAME[diet], if(which.maps == "scrambled") sprintf("%d maps with genotypes scrambled among strains", N.SCRAMBLED) else "real data"), 3, line = 0.5, cex = 0.8)
}
invisible(dev.off())

sink("output/ctl_inflation.txt", split = TRUE)
cat("CTL test under three assumptions, BXD mice, Danny's design. Scrambled: genotypes permuted among strains (", N.SCRAMBLED, "maps)\n")
cat("Rates are over distinct strain patterns x T-ages; a valid test gives lambda 1, P < 0.001 in 0.1% and P < 1e-4 in 0.01%.\n\n")
for(which.maps in c("scrambled", "real")) for(diet in c("CD", "HF")) for(t in names(TEST.COLOR)){
  x <- (if(which.maps == "scrambled") scrambled else real)[[diet]][[t]]; x <- x[!is.na(x)]
  cat(sprintf("%-9s %-2s %-13s tests %8d  lambda %.2f  P<0.001 %6.3f%%  P<1e-4 %6.3f%%  largest %.2f\n", which.maps, diet, t, length(x), lambda(x),
              100 * mean(x > 3), 100 * mean(x > 4), max(x)))
}
sink()
