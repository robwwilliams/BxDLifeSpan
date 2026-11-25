#
# ctl_grid_summary.R - CTL scans of every weighing age (W) against every survivorship (T), with the F1s, judged together
#
# Rob Williams's fork, 25 November 2026. Exploratory. Written by Claude (Opus 5.5). Rob: "the full actuarial scan across all ages at which
# animals were weighed (overkill)". Reads output/ctl_permutations_fixed<W>_F1_maxima.txt (W from 80 to 740 in 30-day steps). Permutation
# i uses the same assignment of genotypes at every W (the seed depends on i only), so the largest value over W and T can be compared
# with its own permutation distribution.
# Writes output/ctl_grid_summary.txt
#   Rscript ctl_grid_summary.R
#

out <- c(); best <- list()
for(W in seq(80, 740, 30)){
  f <- sprintf("output/ctl_permutations_fixed%d_F1_maxima.txt", W); if(!file.exists(f)){ cat("missing W", W, "\n"); next }
  x <- read.delim(f, check.names = FALSE); stopifnot(x[1, "perm"] == 0)
  for(d in c("CD", "HF")){
    k <- paste0(d, ".any"); t <- grep(paste0("^", d, "\\.T"), colnames(x), value = TRUE); if(length(t) == 0) next
    pT <- sapply(t, function(m) (1 + sum(x[-1, m] >= x[1, m])) / nrow(x))
    out <- rbind(out, data.frame(W = W, diet = d, scans = length(t), real = round(x[1, k], 2), P05 = round(unname(quantile(x[-1, k], .95)), 2), P.any.T = round((1 + sum(x[-1, k] >= x[1, k])) / nrow(x), 3), scans.P05 = sum(pT <= 0.05)))
    best[[d]] <- cbind(best[[d]], x[, k])
  }
}
sink("output/ctl_grid_summary.txt", split = TRUE)
cat("CTL, weight at W, survivorships T from W on, BXD and F1 mice, 1,000 hierarchical permutations\n"); options(width = 120); print(out, row.names = FALSE)
for(d in names(best)){ r <- apply(best[[d]], 1, max); cat(sprintf("%s: largest over all W and T: %.2f; permutation 5%% threshold %.2f; P = %.3f\n", d, r[1], quantile(r[-1], .95), (1 + sum(r[-1] >= r[1])) / length(r))) }
cat(sprintf("All scans: %d; with P <= 0.05 on their own: %d (about %.0f expected by chance)\n", sum(out$scans), sum(out$scans.P05), 0.05 * sum(out$scans)))
sink()
