#
# compare_F1_update.R - The maps with F1s before and after the F1 table of 24 November 2026, and with ear tags 3694 and 4520
#
# Rob Williams's fork, 24 November 2026. Written by Claude (Opus 5.5). For each map: correlation over markers and T-ages with the
# maps of the December 2024 table and with the sensitivity run, the largest difference, and the strongest marker (new, old, sensitivity).
#   Rscript compare_F1_update.R <folder with the old maps> <folder with the sensitivity maps>
# Writes output/compare_F1_update.txt
#

a <- commandArgs(trailingOnly = TRUE)
g <- read.table("data/BXD_current_rev050423.geno", sep = "\t", skip = 23, header = TRUE, check.names = FALSE); rownames(g) <- g[, "Locus"]
rd <- function(d, f) as.matrix(read.table(file.path(d, f), sep = "\t"))
pk <- function(x, m){ k <- which(abs(x) == max(abs(x[, ]), na.rm = TRUE), arr.ind = TRUE)[1, ]; sprintf("Chr %s %.1f Mb %s %+.2f", g[m[k[1]], "Chr"], g[m[k[1]], "Mb"], sub("X", "T", colnames(x)[k[2]]), x[k[1], k[2]]) }
cat1 <- c()
for(f in c("qtl2_CD_meansF1_add", "qtl2_HF_meansF1_add", "qtl2_CD_indF1", "qtl2_HF_indF1", "qtl2_main_indF1", "qtl2_int_indF1")){
  F <- paste0(f, "_Progressive.txt"); n <- rd("output", F); tc <- grep("^X[0-9]+$", colnames(n)); n <- n[, tc]; mode(n) <- "numeric"
  for(w in seq_along(a)){
    o <- rd(a[w], F)[rownames(n), colnames(n)]; mode(o) <- "numeric"
    cat1 <- rbind(cat1, data.frame(map = f, against = c("old table", "with 3694 and 4520")[w], r = round(cor(as.vector(n), as.vector(o), use = "complete.obs"), 4),
                                   largest.difference = round(max(abs(n - o), na.rm = TRUE), 2), peak.new = pk(n, rownames(n)), peak.other = pk(o, rownames(n))))
  }
}
options(width = 220); print(cat1, row.names = FALSE)
write.table(cat1, "output/compare_F1_update.txt", sep = "\t", quote = FALSE, row.names = FALSE)
