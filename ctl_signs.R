#
# ctl_signs.R - The sign of the CTL correlations, by diet and T-age, for the figure legend
#
# Rob Williams's fork, 21 November 2026. Reads CTL.R's per-marker correlations (output/CTL/<diet>_Wcor<T>.txt: the rank
# correlation of body weight at T with adjusted lifespan, in B carriers and in D carriers) and reports their medians over
# markers. The CTL maps are coloured by the sign of r(D) - r(B), so warm means the correlation is more positive in D
# carriers, whether both are negative (heavier mice die sooner) or both positive.
#   Rscript ctl_signs.R
# Writes output/ctl_signs.txt
#

res <- NULL
for(diet in c("CD", "HF")){
  for(f in list.files("output/CTL", pattern = paste0("^", diet, "_Wcor"), full.names = TRUE)){
    x <- read.table(f)
    res <- rbind(res, data.frame(diet = diet, T = as.numeric(sub(".*Wcor([0-9]+).txt", "\\1", f)), markers = nrow(x),
                                 median.rB = round(median(x[, 1]), 3), median.rD = round(median(x[, 2]), 3),
                                 both.negative = round(mean(x[, 1] < 0 & x[, 2] < 0), 3), both.positive = round(mean(x[, 1] > 0 & x[, 2] > 0), 3)))
  }
}
res <- res[order(res[, "diet"], res[, "T"]), ]
write.table(res, "output/ctl_signs.txt", sep = "\t", quote = FALSE, row.names = FALSE)
print(res, row.names = FALSE)
