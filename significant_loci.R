#
# significant_loci.R - Genome-wide thresholds of the strain-means maps at three levels, and the loci that pass them
#
# Rob Williams's fork, 21 November 2026. Exploratory. From the 1,000 permutations of strains of qtl2_permutations.R:
#   - significant, genome-wide P = 0.05; and P = 0.10, as Rob asked;
#   - suggestive, genome-wide P = 0.63: a result expected by chance about once per genome scan (Lander & Kruglyak 1995,
#     Nat Genet 11:241), the 37th percentile of the permutation maxima.
# Each level is given per scan (the threshold for one T-age; median over T-ages, with its range) and actuarial (any marker
# at any T-age). For each diet and chromosome, the strongest marker of the strain-means map (equal weights, LOCO
# kinship; qtl2_strainmeans.R) is then set against the per-scan threshold at each T-age, and the chromosomes that reach
# at least the suggestive level at one or more T-ages are listed.
#   Rscript significant_loci.R
# Writes output/significant_loci.txt
#

source("shared.R")
tps <- seq(20, 770, 30)
levels <- c(significant = 0.95, `P 0.10` = 0.90, suggestive = 0.37)
r <- read.table("output/qtl2_permutation_maxima.txt", sep = "\t", header = TRUE)
sink("output/significant_loci.txt", split = TRUE)
cat("Strain-means maps (equal weights, LOCO kinship): genome-wide thresholds (-log10 P) from", sum(r[, "diet"] == "CD"), "permutations of strains\n")
cat("Levels: significant (genome-wide P = 0.05), P = 0.10, and suggestive (P = 0.63; Lander & Kruglyak 1995)\n")
for(d in c("CD", "HF")){
  rr <- as.matrix(r[r[, "diet"] == d, -(1:2)])
  scan <- sapply(levels, function(q) apply(rr, 2, quantile, q))                      # T-ages x levels
  act <- sapply(levels, function(q) quantile(apply(rr, 1, max), q, names = FALSE))
  cat(sprintf("\n%s thresholds:\n", d))
  for(l in names(levels)) cat(sprintf("  %-11s one T-age: median %.2f (range %.2f to %.2f over T-ages); any T-age: %.2f\n", l, median(scan[, l]), min(scan[, l]), max(scan[, l]), act[l]))
  x <- read.table(paste0("output/qtl2_", d, "_means_equal_Progressive.txt"), sep = "\t")
  L <- as.matrix(x[rownames(map), paste0("X", tps)])
  cat(sprintf("%s loci reaching at least the suggestive level for their T-age (strongest marker on the chromosome; allele that lives longer):\n", d))
  for(ch in c(1:19, "X")){
    k <- which(map[, "Chr"] == ch)
    a <- abs(L[k, , drop = FALSE]); best <- apply(a, 2, max, na.rm = TRUE)
    lev <- ifelse(best >= scan[, "significant"], "significant", ifelse(best >= scan[, "P 0.10"], "P 0.10", ifelse(best >= scan[, "suggestive"], "suggestive", "")))
    if(all(lev == "")) next
    w <- which(a == max(a, na.rm = TRUE), arr.ind = TRUE)[1, ]
    top <- if(any(lev == "significant")) "significant" else if(any(lev == "P 0.10")) "P 0.10" else "suggestive"
    cat(sprintf("  Chr %-2s %6.1f Mb %s: peak %.2f at T%d; best level %s; significant at %s; any T-age significant: %s\n", ch, map[k[w[1]], "Mb"],
                if(L[k[w[1]], w[2]] > 0) "D" else "B", a[w[1], w[2]], tps[w[2]], top,
                if(any(lev == "significant")) paste0("T", tps[lev == "significant"], collapse = " ") else "none", if(max(a, na.rm = TRUE) >= act["significant"]) "yes" else "no"))
  }
}
sink()
