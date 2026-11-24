#
# significant_loci_F1.R - The loci of the maps of strain and F1 means together (qtl2_strainF1.R) that pass the genome-wide
# thresholds of qtl2_strainF1_perm.R, additive and dominance
#
# Rob Williams's fork, 23 November 2026. Exploratory. The same rules as significant_loci.R for the strain-means maps:
#   - levels: significant (genome-wide P = 0.05), P = 0.10, and suggestive (P = 0.63; Lander & Kruglyak 1995);
#   - per scan (one T-age) and any T-age; a stable range of T-ages for each diet and map, from T20 to the last T-age
#     before the 5% per-scan threshold first moves more than 0.3 from its median over T20-T500;
#   - for each chromosome, the strongest marker against the per-scan threshold at each T-age in the stable range.
# Beyond the stable range, chromosomes whose strongest marker passes the any-T-age 5% threshold of the stable range (the
# strictest of them) are listed as well, with the per-scan threshold at their own T-age.
#   Rscript significant_loci_F1.R
# Writes output/significant_loci_F1.txt
#

source("shared.R")
tps <- seq(20, 770, 30)
levels <- c(significant = 0.95, `P 0.10` = 0.90, suggestive = 0.37)
r <- read.delim("output/qtl2_strainF1_perm_maxima.txt", check.names = FALSE)
sink("output/significant_loci_F1.txt", split = TRUE)
cat("Strain and F1 means together (qtl2_strainF1.R): genome-wide thresholds (-log10 P) from", sum(r$diet == "CD" & r$map == "add" & r$perm > 0),
    "permutations (qtl2_strainF1_perm.R)\n")
cat("Levels: significant (genome-wide P = 0.05), P = 0.10, and suggestive (P = 0.63; Lander & Kruglyak 1995)\n")
for(d in c("CD", "HF")) for(m in c("add", "dom")){
  rr <- as.matrix(r[r$diet == d & r$map == m & r$perm > 0, paste0("T", tps)])
  s5 <- apply(rr, 2, quantile, 0.95); ref <- median(s5[tps <= 500])
  off <- which(abs(s5 - ref) > 0.3 & tps > 500)
  tmax <- if(length(off) == 0) max(tps) else tps[min(off) - 1]
  tt <- tps[tps <= tmax]
  scan <- sapply(levels, function(q) apply(rr[, tps <= tmax, drop = FALSE], 2, quantile, q))
  act <- sapply(levels, function(q) quantile(apply(rr[, tps <= tmax, drop = FALSE], 1, max), q, names = FALSE))
  sgn <- if(m == "add") c("D", "B") else c("BD > midpoint", "BD < midpoint")
  cat(sprintf("\n%s, %s: thresholds over the stable range T20 to T%d\n", d, c(add = "additive", dom = "dominance")[m], tmax))
  for(l in names(levels)) cat(sprintf("  %-11s one T-age: median %.2f (range %.2f to %.2f over T-ages); any T-age: %.2f\n", l, median(scan[, l]), min(scan[, l]), max(scan[, l]), act[l]))
  x <- read.table(paste0("output/qtl2_", d, "_meansF1_", m, "_Progressive.txt"), sep = "\t")
  L <- as.matrix(x[rownames(map), paste0("X", tps)]); L[is.na(L)] <- 0
  cat("  Loci reaching at least the suggestive level for their T-age (strongest marker on the chromosome):\n")
  for(ch in c(1:19, "X")){
    k <- which(map[, "Chr"] == ch)
    a <- abs(L[k, tps <= tmax, drop = FALSE]); best <- apply(a, 2, max)
    lev <- ifelse(best >= scan[, "significant"], "significant", ifelse(best >= scan[, "P 0.10"], "P 0.10", ifelse(best >= scan[, "suggestive"], "suggestive", "")))
    if(all(lev == "")) next
    w <- which(a == max(a), arr.ind = TRUE)[1, ]
    top <- if(any(lev == "significant")) "significant" else if(any(lev == "P 0.10")) "P 0.10" else "suggestive"
    cat(sprintf("    Chr %-2s %6.1f Mb %s: peak %.2f at T%d; best level %s; significant at %s; any T-age significant: %s\n", ch, map[k[w[1]], "Mb"],
                if(L[k[w[1]], w[2]] > 0) sgn[1] else sgn[2], a[w[1], w[2]], tt[w[2]], top,
                if(any(lev == "significant")) paste0("T", tt[lev == "significant"], collapse = " ") else "none", if(max(a) >= act["significant"]) "yes" else "no"))
  }
  if(tmax < max(tps)){
    cat(sprintf("  Beyond the stable range (T%d to T770), chromosomes whose strongest marker passes the any-T-age 5%% threshold (%.2f):\n", tmax + 30, act["significant"]))
    hit <- FALSE
    for(ch in c(1:19, "X")){
      k <- which(map[, "Chr"] == ch); a <- abs(L[k, tps > tmax, drop = FALSE])
      if(max(a) < act["significant"]) next
      hit <- TRUE; w <- which(a == max(a), arr.ind = TRUE)[1, ]; T <- tps[tps > tmax][w[2]]
      cat(sprintf("    Chr %-2s %6.1f Mb %s: peak %.2f at T%d (5%% threshold at T%d alone: %.2f)\n", ch, map[k[w[1]], "Mb"],
                  if(L[k[w[1]], which(tps == T)] > 0) sgn[1] else sgn[2], a[w[1], w[2]], T, T, s5[paste0("T", T)]))
    }
    if(!hit) cat("    none\n")
  }
}
sink()
