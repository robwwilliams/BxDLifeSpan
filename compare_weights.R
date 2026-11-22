#
# compare_weights.R - The three weightings of the strain-means maps (qtl2_strainmeans.R) against each other
#
# Rob Williams's fork, 21 November 2026. Exploratory. At each diet and T-age: the correlation over markers (signed
# -log10 P) of equal weights with N weights and with variance weights, and the strongest marker of each map. Also the
# strongest -log10 P on each chromosome over all T-ages, for each weighting.
#   Rscript compare_weights.R
# Writes output/compare_weights.txt
#

source("shared.R")
tps <- seq(20, 770, 30)
schemes <- c("equal", "N", "var")
peak <- function(v){ i <- which.max(abs(v)); sprintf("Chr %s %.1f Mb %s %.2f", map[i, "Chr"], map[i, "Mb"], if(v[i] > 0) "D" else "B", abs(v[i])) }

sink("output/compare_weights.txt", split = TRUE)
cat("Strain-means maps with LOCO kinship (qtl2_strainmeans.R): equal weights against N weights (Danny's) and variance weights.\n")
cat("r: correlation over markers of signed -log10 P. Peaks: strongest marker, with the allele that raises lifespan.\n")
for(d in c("CD", "HF")){
  m <- lapply(schemes, function(s){ x <- read.table(paste0("output/qtl2_", d, "_means_", s, "_Progressive.txt"), sep = "\t"); as.matrix(x[rownames(map), paste0("X", tps)]) })
  names(m) <- schemes
  cat("\n=====", d, "\n")
  print(data.frame(T = tps,
                   r.equal.N = sapply(seq_along(tps), function(i) round(cor(m$equal[, i], m$N[, i], use = "complete"), 3)),
                   r.equal.var = sapply(seq_along(tps), function(i) round(cor(m$equal[, i], m$var[, i], use = "complete"), 3)),
                   peak.equal = sapply(seq_along(tps), function(i) peak(m$equal[, i])),
                   peak.N = sapply(seq_along(tps), function(i) peak(m$N[, i])),
                   peak.var = sapply(seq_along(tps), function(i) peak(m$var[, i]))), row.names = FALSE)
  top <- function(L) sapply(c(1:19, "X"), function(ch){ k <- which(map[, "Chr"] == ch); a <- abs(L[k, ]); w <- which(a == max(a, na.rm = TRUE), arr.ind = TRUE)[1, ]
                                                       sprintf("%.2f (%.0f Mb, T%d, %s)", a[w[1], w[2]], map[k[w[1]], "Mb"], tps[w[2]], if(L[k[w[1]], w[2]] > 0) "D" else "B") })
  cat("\nStrongest -log10 P on each chromosome over all T-ages:\n")
  print(data.frame(chr = c(1:19, "X"), equal = top(m$equal), N = top(m$N), var = top(m$var)), row.names = FALSE)
}
sink()
