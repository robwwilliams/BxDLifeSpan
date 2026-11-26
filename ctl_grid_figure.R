#
# SUPERSEDED on 26 November 2026 by ctl_scan.R (BXD only, one path, tested by ctl_tests.R). Kept so the committed outputs it made can be traced.
# ctl_grid_figure.R - One heat map of the CTL scans: every weighing age against every later survivorship
#
# Rob Williams's fork, 25 November 2026. Exploratory. Written by Claude (Opus 5.5). Each cell is one genome scan, colored by the genome-wide P
# of its strongest marker from 1,000 hierarchical permutations (ctl_permutations.R, fixed W, with F1s). A null family of scans colors
# about 5% of its cells at P <= 0.05. Writes output/maps/CTL_grid.svg
#   Rscript ctl_grid_figure.R
#

library(svglite)
Ws <- seq(80, 740, 30); Ts <- seq(80, 770, 30); P <- list()
for(d in c("CD", "HF")){ m <- matrix(NA, length(Ws), length(Ts), dimnames = list(Ws, Ts))
  for(W in Ws){ x <- read.delim(sprintf("output/ctl_permutations_fixed%d_F1_maxima.txt", W), check.names = FALSE)
    for(k in grep(paste0("^", d, "\\.T"), colnames(x), value = TRUE)) m[as.character(W), sub(".*\\.T", "", k)] <- (1 + sum(x[-1, k] >= x[1, k])) / nrow(x) }
  P[[d]] <- m }
tiers <- c(0, 0.05, 0.10, 0.63, 1.01); lab <- c("≤ 0.05", "0.05–0.10", "0.10–0.63", "> 0.63")
col <- list(CD = c("#6b3f10", "#a86a28", "#dcb98f", "#f4ece1"), HF = c("#123f7a", "#2a78d6", "#a9c9f0", "#e8f0fb"))
name <- c(CD = "Chow", HF = "High fat")
svglite("output/maps/CTL_grid.svg", width = 7.2, height = 4.1, bg = "white", system_fonts = list(sans = "Arial"), fix_text_size = FALSE)
layout(matrix(1:2, 1)); par(mai = c(0.6, 0.6, 0.4, 0.1), mgp = c(2, 0.5, 0), tcl = -0.2, las = 1, cex = 0.75)
for(d in c("CD", "HF")){
  plot(NA, xlim = c(65, 785), ylim = c(65, 785), xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "l", xlab = "Age at weighing (days)", ylab = "Survivorship, T-age (days)")
  axis(1, seq(80, 740, 120)); axis(2, seq(80, 770, 120))
  for(i in seq_along(Ws)) for(j in seq_along(Ts)) if(!is.na(P[[d]][i, j])) rect(Ws[i] - 15, Ts[j] - 15, Ws[i] + 15, Ts[j] + 15, border = NA, col = col[[d]][findInterval(P[[d]][i, j], tiers, left.open = TRUE)])
  mtext(letters[match(d, c("CD", "HF"))], 3, line = 0.6, at = grconvertX(0.02, "nfc", "user"), adj = 0, font = 2, cex = 1)
  mtext(name[d], 3, line = 0.6, font = 2, cex = 0.85)
  text(480, 300 - 45 * seq_along(lab), lab, adj = 0, cex = 0.9); points(rep(450, 4), 300 - 45 * seq_along(lab), pch = 22, bg = col[[d]], col = "grey60", cex = 1.6)
  text(440, 305, "Genome-wide P", adj = 0, cex = 0.9)
  cat(d, "cells", sum(!is.na(P[[d]])), "<=0.05:", sum(P[[d]] <= 0.05, na.rm = TRUE), "\n")
}
invisible(dev.off())
