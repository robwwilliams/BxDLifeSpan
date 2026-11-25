#
# weight_lifespan_sign.R - When does the correlation of body weight with lifespan change sign, and between or within strains?
#
# Rob Williams's fork, 24 November 2026. Exploratory. Written by Claude (Opus 5.5). Rob, of the sign change on high fat that the CTL
# legend put at about 290 days: "That is about 150 days earlier than I would have guessed based on Arends et al. and UM-HET3 data."
# The CTL correlations (CTL.R) pool mice of all strains, so they mix differences between strains with differences among the mice of a
# strain. Here, at each age T from 80 to 740 days (mice alive at T and on their diet by T; weights interpolated as in bodyweight.R):
#   pooled   rank correlation, each mouse weighted by 1 / mice of its strain, as in CTL.R (no adjustment for vivarium and set-up)
#   within   rank correlation of weight and lifespan, each centered on its strain's mean
#   gap      the within-strain correlation for mice that lived at least 100 days past T, so that wasting mice are left out
# 95% intervals: 1,000 bootstrap resamples of strains (the strain is the unit). Fixed seed.
#   Rscript weight_lifespan_sign.R
# Writes output/weight_lifespan_sign.txt
#

source("shared.R")
set.seed(20261124)
B <- read.table("output/bodyweightInterpolatedData.txt", sep = "\t", quote = "", check.names = FALSE)
ldata[which(ldata[, "AgeAtSetUp.in.colony..days."] < 0), "AgeAtSetUp.in.colony..days."] <- NA
l <- ldata[grepl("BXD", ldata[, "StrainName"]) & rownames(ldata) %in% rownames(B), ]
pooled <- function(w, y, s){ wt <- as.numeric(1 / table(s)[s]); cov.wt(cbind(rank(w), rank(y)), wt = wt, cor = TRUE)$cor[1, 2] }
within <- function(w, y, s){ k <- ave(w, s, FUN = length) > 1; cor((w - ave(w, s))[k], (y - ave(y, s))[k], method = "spearman") }
boot <- function(f, w, y, s){
  idx <- split(seq_along(s), s)
  b <- replicate(1000, { k <- sample(names(idx), replace = TRUE); i <- unlist(idx[k], use.names = FALSE); f(w[i], y[i], as.character(rep(seq_along(k), lengths(idx[k])))) })
  c(f(w, y, s), quantile(b, c(0.025, 0.975), names = FALSE))
}
out <- c()
for(d in c("CD", "HF")) for(T in seq(80, 740, 30)){
  x <- l[which(l[, "Diet"] == d & l[, "AgeAtSetUp.in.colony..days."] <= T & l[, "AgeAtDeath..days."] >= T), ]
  w <- B[rownames(x), as.character(T)]; ok <- !is.na(w); w <- w[ok]; y <- x[ok, "AgeAtDeath..days."]; s <- as.character(x[ok, "StrainName"])
  g <- y >= T + 100
  p <- boot(pooled, w, y, s); wi <- boot(within, w, y, s); ga <- boot(within, w[g], y[g], s[g])
  out <- rbind(out, data.frame(diet = d, T = T, mice = length(w), strains = length(unique(s)), pooled = p[1], pooled.lo = p[2], pooled.hi = p[3],
                               within = wi[1], within.lo = wi[2], within.hi = wi[3], mice.gap = sum(g), gap = ga[1], gap.lo = ga[2], gap.hi = ga[3]))
}
for(k in setdiff(colnames(out), c("diet", "T", "mice", "strains", "mice.gap"))) out[, k] <- round(out[, k], 2)
options(width = 200); print(out, row.names = FALSE)
write.table(out, "output/weight_lifespan_sign.txt", sep = "\t", quote = FALSE, row.names = FALSE)
