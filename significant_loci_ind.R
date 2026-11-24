#
# significant_loci_ind.R - The loci of the individual-mouse maps that pass genome-wide thresholds from permutations of
# strains (and F1 types) run on ATHENA through Squirm
#
# Rob Williams's fork, 23 November 2026. Exploratory. Thresholds: QTL/qtl2_individuals_perm.R (BXD mice; 200
# permutations, squirm/bxd_individual_perms, ATHENA, R 4.6.2) and, once returned, QTL/qtl2_individualsF1_perm.R (BXD and
# F1 mice) in bxd_mortality_rates_2026; their maxima are copied into output/.
#   - Levels: significant (genome-wide P = 0.05), P = 0.10, suggestive (P = 0.63). 200 permutations give a crude 5% level.
#   - Each T-age is judged against its own threshold. The thresholds fall steadily with T-age, as the survivorships
#     shrink (from about 6 at T20 to 4.6 at T770), rather than becoming unstable as the strain-mean ones did, so no
#     stable range is cut. Any T-age: the largest value over all 26 T-ages.
#   - For each map and chromosome, the strongest marker against the threshold at each T-age.
#   Rscript significant_loci_ind.R
# Writes output/significant_loci_ind.txt
#

source("shared.R")
tps <- seq(20, 770, 30)
levels <- c(significant = 0.95, `P 0.10` = 0.90, suggestive = 0.37)
sets <- list(BXD = list(maxima = "output/qtl2_individuals_perm_maxima_ATHENA.txt", maps = "output/qtl2_%s_ind_Progressive.txt", label = "BXD mice"),
             F1 = list(maxima = "output/qtl2_individualsF1_perm_maxima_ATHENA.txt", maps = "output/qtl2_%s_indF1_Progressive.txt", label = "BXD and F1 mice"))
sink("output/significant_loci_ind.txt", split = TRUE)
for(st in sets){
  if(!file.exists(st$maxima)){ cat("\n", st$label, ": no permutations yet (", st$maxima, ")\n", sep = ""); next }
  r <- read.delim(st$maxima, check.names = FALSE)
  cat(sprintf("\nIndividual-mouse maps, %s: thresholds from %d permutations\n", st$label, length(unique(r$perm))))
  for(m in c("CD", "HF", "main", "int")){
    rr <- as.matrix(r[r$map == m, paste0("T", tps)])
    scan <- sapply(levels, function(q) apply(rr, 2, quantile, q)); act <- sapply(levels, function(q) quantile(apply(rr, 1, max), q, names = FALSE))
    cat(sprintf("  %s: 5%% one T-age %.2f (T20) to %.2f (T770); any T-age %.2f. 10%%: %.2f to %.2f; any %.2f\n", m, scan[1, 1], scan[nrow(scan), 1], act[1], scan[1, 2], scan[nrow(scan), 2], act[2]))
    x <- read.table(sprintf(st$maps, m), sep = "\t"); L <- as.matrix(x[rownames(map), paste0("X", tps)]); L[is.na(L)] <- 0
    sgn <- if(m == "int") c("D - B larger on HF", "D - B smaller on HF") else c("D", "B")
    for(ch in c(1:19, "X")){
      k <- which(map[, "Chr"] == ch); a <- abs(L[k, , drop = FALSE]); best <- apply(a, 2, max)
      lev <- ifelse(best >= scan[, 1], "significant", ifelse(best >= scan[, 2], "P 0.10", ifelse(best >= scan[, 3], "suggestive", "")))
      if(all(lev == "")) next
      w <- which(a == max(a), arr.ind = TRUE)[1, ]
      top <- if(any(lev == "significant")) "significant" else if(any(lev == "P 0.10")) "P 0.10" else "suggestive"
      hit <- which(lev == top); h <- hit[which.max(best[hit])]
      cat(sprintf("    Chr %-2s %6.1f Mb %s: peak %.2f at T%d; best level %s (%s, %.2f at T%d against %.2f); any T-age significant: %s\n", ch, map[k[w[1]], "Mb"],
                  if(L[k[w[1]], w[2]] > 0) sgn[1] else sgn[2], a[w[1], w[2]], tps[w[2]], top, paste0("T", tps[hit], collapse = " "), best[h], tps[h],
                  scan[h, c(significant = 1, `P 0.10` = 2, suggestive = 3)[top]], if(max(a) >= act[1]) "yes" else "no"))
    }
  }
}
sink()
