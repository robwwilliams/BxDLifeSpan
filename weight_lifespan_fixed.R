#
# weight_lifespan_fixed.R - The correlation of an early body weight with lifespan, as the survivorship rises
#
# Rob Williams's fork, 25 November 2026. Exploratory. Written by Claude (Opus 5.5). Rob: "did you take body weight at 3 months and then
# compute lifespan:weight correlations of the survivorship up toward 800 days?" Weight is held at W (90, 180, 365 or 545 days); the mice are
# those alive at T, from W to 770 days in 30-day steps. As in weight_lifespan_sign.R: pooled (each mouse weighted by 1 / mice of its strain)
# and within strain (both centered on the strain mean); 95% intervals from 1,000 bootstrap resamples of strains. BXD mice, no F1s.
# Writes output/weight_lifespan_fixed.txt
#   Rscript weight_lifespan_fixed.R
#

source("shared.R"); set.seed(20261125)
B <- read.table("output/bodyweightInterpolatedData.txt", sep = "\t", quote = "", check.names = FALSE)
ldata[which(ldata[, "AgeAtSetUp.in.colony..days."] < 0), "AgeAtSetUp.in.colony..days."] <- NA
l <- ldata[grepl("BXD", ldata[, "StrainName"]) & rownames(ldata) %in% rownames(B), ]
pooled <- function(w, y, s){ wt <- as.numeric(1 / table(s)[s]); cov.wt(cbind(rank(w), rank(y)), wt = wt, cor = TRUE)$cor[1, 2] }
within <- function(w, y, s){ k <- ave(w, s, FUN = length) > 1; cor((w - ave(w, s))[k], (y - ave(y, s))[k], method = "spearman") }
boot <- function(f, w, y, s){ idx <- split(seq_along(s), s)
  b <- replicate(1000, { k <- sample(names(idx), replace = TRUE); i <- unlist(idx[k], use.names = FALSE); f(w[i], y[i], as.character(rep(seq_along(k), lengths(idx[k])))) })
  c(f(w, y, s), quantile(b, c(0.025, 0.975), names = FALSE, na.rm = TRUE)) }
out <- c()
for(d in c("CD", "HF")) for(W in c(90, 180, 365, 545)) for(T in unique(c(W, seq(ceiling((W - 20) / 30) * 30 + 20, 770, 30)))){
  x <- l[which(l[, "Diet"] == d & l[, "AgeAtSetUp.in.colony..days."] <= W & l[, "AgeAtDeath..days."] >= T), ]
  w <- B[rownames(x), as.character(W)]; ok <- !is.na(w); w <- w[ok]; y <- x[ok, "AgeAtDeath..days."]; s <- as.character(x[ok, "StrainName"])
  if(length(w) < 60) next
  p <- boot(pooled, w, y, s); wi <- boot(within, w, y, s)
  out <- rbind(out, data.frame(diet = d, W = W, T = T, mice = length(w), strains = length(unique(s)), pooled = p[1], pooled.lo = p[2], pooled.hi = p[3], within = wi[1], within.lo = wi[2], within.hi = wi[3]))
}
for(k in c("pooled", "pooled.lo", "pooled.hi", "within", "within.lo", "within.hi")) out[, k] <- round(out[, k], 2)
write.table(out, "output/weight_lifespan_fixed.txt", sep = "\t", quote = FALSE, row.names = FALSE)
options(width = 200); print(out[out$T %in% c(90, 180, 365, 545, 200, 290, 380, 500, 590, 680, 770), ], row.names = FALSE)
