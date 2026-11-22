#
# significant_loci.R - Strain-means loci that pass the genome-wide thresholds of qtl2_permutations.R
#
# Rob Williams's fork, 21 November 2026. Exploratory. For each diet and chromosome, the strongest marker of the strain-means
# map (equal weights, LOCO kinship; qtl2_strainmeans.R) at each T-age, set against the 5% threshold for that T-age (per
# scan) and against the 5% threshold for any marker at any T-age (actuarial), both from 1,000 permutations of strains.
# Lists every chromosome that passes the per-scan threshold at one or more T-ages.
#   Rscript significant_loci.R
# Writes output/significant_loci.txt
#

source("shared.R")
tps <- seq(20, 770, 30)
r <- read.table("output/qtl2_permutation_maxima.txt", sep = "\t", header = TRUE)
sink("output/significant_loci.txt", split = TRUE)
cat("Strain-means maps (equal weights, LOCO kinship): chromosomes whose strongest marker passes the 5% genome-wide threshold\n")
cat("for its T-age (per scan), and whether it also passes the 5% threshold for any T-age (actuarial). -log10 P; allele that lives longer.\n")
for(d in c("CD", "HF")){
  rr <- as.matrix(r[r[, "diet"] == d, -(1:2)])
  scan5 <- apply(rr, 2, quantile, 0.95); act5 <- quantile(apply(rr, 1, max), 0.95)
  x <- read.table(paste0("output/qtl2_", d, "_means_equal_Progressive.txt"), sep = "\t")
  L <- as.matrix(x[rownames(map), paste0("X", tps)])
  cat(sprintf("\n%s: 5%% per scan %.2f to %.2f over T-ages; 5%% actuarial %.2f\n", d, min(scan5), max(scan5), act5))
  for(ch in c(1:19, "X")){
    k <- which(map[, "Chr"] == ch)
    a <- abs(L[k, , drop = FALSE])
    pass <- which(apply(a, 2, max, na.rm = TRUE) >= scan5)
    if(length(pass) == 0) next
    w <- which(a == max(a, na.rm = TRUE), arr.ind = TRUE)[1, ]
    cat(sprintf("  Chr %-2s %6.1f Mb %s: peak %.2f at T%d (per-scan threshold there %.2f); passes per scan at %s; actuarial: %s\n", ch, map[k[w[1]], "Mb"],
                if(L[k[w[1]], w[2]] > 0) "D" else "B", a[w[1], w[2]], tps[w[2]], scan5[w[2]], paste0("T", tps[pass], collapse = " "),
                if(max(a, na.rm = TRUE) >= act5) "yes" else "no"))
  }
}
sink()
