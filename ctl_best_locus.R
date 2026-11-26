#
# ctl_best_locus.R - The strongest CTL locus of each diet, displayed as in Fig. 4i-p of Arends et al. 2026: the weight-lifespan
# correlation of B and D mice at every T-age, and the mice and strains behind it at the peak T-age
#
# Rob Williams's fork, 26 November 2026. Exploratory, asked for by Danny Arends ("displays of the CTLs like we have in Figure 4 of the
# Nature paper ... to actually see what the slope differences are"). Written by Claude (Opus 5.5). BXD mice only, Danny's design
# (weight at T). Run from the fork's root after ctl_inputs.R: Rscript ctl_best_locus.R
#
# Steps: 1 the real CTL map (ctl_functions.R); 2 the strongest marker and T-age per diet; 3 its correlations at every T-age, with 50% CI
# from Fisher's z; 4 the figure; 5 a table.
#   The 50% CI band uses the number of strains as n: the strain is the unit of inference, as mice of a strain share their genome. The
#   table also gives the mice (Danny's n) and the Kish effective number of the weighted mice, (sum w)^2 / sum w^2 with w = 1 / (mice of
#   the strain), which allows for the unequal weights but not for strain mates sharing a genome.
# Writes output/maps/CTL_best_locus.svg (and .png) and output/ctl_best_locus.txt
#

source("ctl_functions.R"); library(svglite)
DIR <- "output/ctl_scan"
COLOR <- c(B = "#7d4598", D = "#ed8f22")                                   # haplotype colors (notes/figure_style.md)
DIET.NAME <- c(CD = "Chow", HF = "High fat")

# ---- 1. The real map -------------------------------------------------------------------------------------------------------------------
inp <- readInputs(DIR)
Tages <- as.numeric(sub("^w", "", grep("^w[0-9]+$", colnames(inp$mice), value = TRUE)))
scans <- buildScans(inp$mice, inp$strains, c(actuarial = NA), Tages, 2)
real.map <- scanAll(inp$G, scans)                                           # distinct patterns x scans, signed -log10 P

# ---- 2-3. The strongest marker and T-age of each diet, and its correlations at every T-age ------------------------------------------
kish <- function(w) sum(w)^2 / sum(w^2)
best <- list()
for(diet in c("CD", "HF")){
  columns <- which(sapply(scans, `[[`, "diet") == diet)
  peak <- which(abs(real.map[, columns]) == max(abs(real.map[, columns]), na.rm = TRUE), arr.ind = TRUE)[1, ]   # (pattern, scan)
  pattern <- peak[1]; peak.scan <- scans[[columns[peak[2]]]]
  markers <- which(inp$pat == pattern)                                      # every marker with this strain pattern
  by.T <- do.call(rbind, lapply(scans[columns], function(d){
    r <- ctlCor(inp$G[, pattern, drop = FALSE], d$s, d$wt, d$ow, d$tw, d$oy, d$ty, d$nall)[1, ]            # rB, rD, nB, nD
    g <- inp$G[d$s + 1L, pattern]                                           # genotype of each mouse: 1 B, 2 D
    data.frame(Tage = d$Tage, rB = r[1], rD = r[2], nB = r[3], nD = r[4], effB = kish(d$wt[g == 1]), effD = kish(d$wt[g == 2]),
               strainsB = length(unique(d$s[g == 1])), strainsD = length(unique(d$s[g == 2])), log10P = ctlTest(r[1], r[2], r[3], r[4])) }))
  best[[diet]] <- list(pattern = pattern, scan = peak.scan, markers = markers, by.T = by.T)
}

ciBand <- function(r, n, level = 0.5){ z <- 0.5 * log((1 + r) / (1 - r)); h <- qnorm(0.5 + level / 2) / sqrt(n - 3); cbind(lo = tanh(z - h), hi = tanh(z + h)) }

# ---- 4. The figure: per diet, a the correlations by T-age, b the mice at the peak T-age, c the strain means at the peak T-age ----------
dir.create("output/maps", showWarnings = FALSE)
svglite("output/maps/CTL_best_locus.svg", width = 7.2, height = 5.0, bg = "white", system_fonts = list(sans = "Arial"), fix_text_size = FALSE)
layout(matrix(1:6, 2, byrow = TRUE), widths = c(1.25, 1, 1)); par(mai = c(0.5, 0.62, 0.3, 0.5), mgp = c(2.2, 0.45, 0), tcl = -0.2, las = 1, cex = 0.7)
panel <- 0
for(diet in c("CD", "HF")){
  b <- best[[diet]]; x <- b$by.T; d <- b$scan; mk <- inp$map[b$markers, ]
  where <- if(nrow(mk) == 1) sprintf("Chr %s, %.1f Mb", mk$Chr[1], mk$Mb) else sprintf("Chr %s, %.1f-%.1f Mb", mk$Chr[1], min(mk$Mb), max(mk$Mb))
  # a: correlations of B and D mice at every T-age, with 50% CI bands (n = strains)
  plot(NA, xlim = range(Tages), ylim = c(-0.8, 0.8), xlab = "Survivorship, T-age (days)", ylab = "Weight-lifespan correlation (Spearman)", bty = "l")
  abline(h = 0, col = "grey70", lwd = 0.5)
  for(k in c("B", "D")){ r <- x[[paste0("r", k)]]
    band <- ciBand(r, x[[paste0("strains", k)]]); polygon(c(x$Tage, rev(x$Tage)), c(band[, "lo"], rev(band[, "hi"])), col = adjustcolor(COLOR[k], 0.25), border = NA)
    lines(x$Tage, r, col = COLOR[k], lwd = 1.5) }
  abline(v = d$Tage, col = "grey50", lty = 3)
  par(new = TRUE); plot(x$Tage, abs(x$log10P), type = "l", lty = 2, col = "black", lwd = 0.8, axes = FALSE, xlab = "", ylab = "", xlim = range(Tages), ylim = c(0, max(12, abs(x$log10P), na.rm = TRUE)))
  axis(4); mtext("CTL -log10 P (mice as n)", 4, line = 1.7, las = 0, cex = 0.7)
  legend("topleft", c("B mice", "D mice", "-log10 P"), col = c(COLOR, "black"), lty = c(1, 1, 2), lwd = c(1.5, 1.5, 0.8), bty = "n", cex = 0.85)
  panel <- panel + 1; mtext(letters[panel], 3, line = 0.4, at = grconvertX(0.01, "nfc", "user"), adj = 0, font = 2, cex = 1.15)
  mtext(paste(DIET.NAME[diet], "-", where), 3, line = 0.4, cex = 0.75)
  # b: the mice at the peak T-age, adjusted weight against adjusted lifespan
  g <- inp$G[d$s + 1L, b$pattern]; keep <- g %in% 1:2
  plot(d$w[keep], d$y[keep], pch = 16, cex = 0.45, col = adjustcolor(COLOR[g[keep]], 0.55), xlab = sprintf("Adjusted weight at %d days (g)", d$Tage), ylab = "Adjusted lifespan (days)", bty = "l", ylim = range(d$y[keep]) + c(0, 0.25) * diff(range(d$y[keep])))
  for(k in 1:2){ i <- g == k; abline(lm(d$y[i] ~ d$w[i], weights = d$wt[i]), col = COLOR[k], lwd = 1.5) }
  label <- sprintf("%s: r = %+.2f (%d mice, %d strains)", c("B", "D"), c(x$rB[x$Tage == d$Tage], x$rD[x$Tage == d$Tage]), c(sum(g == 1), sum(g == 2)),
                   c(x$strainsB[x$Tage == d$Tage], x$strainsD[x$Tage == d$Tage]))
  text(grconvertX(0.02, "npc"), grconvertY(c(0.98, 0.91), "npc"), label, col = COLOR, adj = c(0, 1), cex = 0.68)
  panel <- panel + 1; mtext(letters[panel], 3, line = 0.4, at = grconvertX(0.01, "nfc", "user"), adj = 0, font = 2, cex = 1.15)
  mtext(sprintf("T%d, each mouse", d$Tage), 3, line = 0.4, cex = 0.75)
  # c: the strain means at the peak T-age, point area by mice
  strain.weight <- tapply(d$w, d$s, mean); strain.life <- tapply(d$y, d$s, mean); strain.mice <- tapply(d$w, d$s, length)
  strain.geno <- inp$G[as.integer(names(strain.weight)) + 1L, b$pattern]; k2 <- strain.geno %in% 1:2
  plot(strain.weight[k2], strain.life[k2], pch = 21, bg = adjustcolor(COLOR[strain.geno[k2]], 0.6), col = "white", lwd = 0.3, cex = 0.5 + sqrt(strain.mice[k2]) / 2.5,
       xlab = sprintf("Strain mean weight at %d days (g)", d$Tage), ylab = "Strain mean lifespan (days)", bty = "l", ylim = range(strain.life[k2]) + c(0, 0.25) * diff(range(strain.life[k2])))
  rs <- sapply(1:2, function(k) cor(strain.weight[strain.geno == k], strain.life[strain.geno == k], method = "spearman"))
  text(grconvertX(0.03, "npc"), grconvertY(c(0.97, 0.90), "npc"), sprintf("%s: r = %+.2f, %d strains", c("B", "D"), rs, c(sum(strain.geno == 1), sum(strain.geno == 2))), col = COLOR, adj = c(0, 1), cex = 0.8)
  panel <- panel + 1; mtext(letters[panel], 3, line = 0.4, at = grconvertX(0.01, "nfc", "user"), adj = 0, font = 2, cex = 1.15)
  mtext(sprintf("T%d, strain means", d$Tage), 3, line = 0.4, cex = 0.75)
  best[[diet]]$strain.r <- rs
}
invisible(dev.off())

# ---- 5. The table -------------------------------------------------------------------------------------------------------------------
sink("output/ctl_best_locus.txt", split = TRUE)
cat("Strongest CTL marker of each diet (BXD mice, weight at T, Danny's test), and its correlations at every T-age\n")
cat("rB, rD: weighted Spearman correlations of weight and lifespan; nB, nD: mice; effB, effD: Kish effective number; strainsB, strainsD.\n")
for(diet in c("CD", "HF")){ b <- best[[diet]]; mk <- inp$map[b$markers, ]
  cat(sprintf("\n%s: peak at T%d, -log10 P %.2f; %d markers share the strain pattern (Chr %s, %.1f-%.1f Mb; first %s). Strain-mean Spearman at the peak: B %+.2f, D %+.2f\n",
              DIET.NAME[diet], b$scan$Tage, max(abs(b$by.T$log10P), na.rm = TRUE), nrow(mk), mk$Chr[1], min(mk$Mb), max(mk$Mb), mk$Locus[1], b$strain.r[1], b$strain.r[2]))
  print(cbind(b$by.T[c("Tage", "nB", "nD", "strainsB", "strainsD")], round(b$by.T[c("rB", "rD", "effB", "effD", "log10P")], 2)), row.names = FALSE)
}
sink()
