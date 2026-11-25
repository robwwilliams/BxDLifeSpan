#
# weight_lifespan_diet.R - Does the correlation of an early body weight with lifespan differ between chow and high fat?
#
# Rob Williams's fork, 25 November 2026. Exploratory. Written by Claude (Opus 5.5). Rob: "the DIFFERENCE at the phenotype level between CD and
# HF cohorts may well be worth reporting". High fat minus chow, of the correlations of weight_lifespan_fixed.R (weight at W; all mice alive at
# W, and survivors of T), pooled and within strain. 95% intervals: 2,000 bootstrap resamples of strains, the same strains on both diets.
# Writes output/weight_lifespan_diet.txt
#   Rscript weight_lifespan_diet.R
#

source("shared.R"); set.seed(20261125)
B <- read.table("output/bodyweightInterpolatedData.txt", sep = "\t", quote = "", check.names = FALSE)
ldata[which(ldata[, "AgeAtSetUp.in.colony..days."] < 0), "AgeAtSetUp.in.colony..days."] <- NA
l <- ldata[grepl("BXD", ldata[, "StrainName"]) & rownames(ldata) %in% rownames(B) & ldata[, "Diet"] %in% c("CD", "HF"), ]
pooled <- function(w, y, s){ wt <- as.numeric(1 / table(s)[s]); cov.wt(cbind(rank(w), rank(y)), wt = wt, cor = TRUE)$cor[1, 2] }
within <- function(w, y, s){ k <- ave(w, s, FUN = length) > 1; cor((w - ave(w, s))[k], (y - ave(y, s))[k], method = "spearman") }
out <- c()
for(W in c(90, 180, 365, 545)) for(T in unique(c(W, 380, 500, 590)[c(W, 380, 500, 590) >= W])){
  x <- l[which(l[, "AgeAtSetUp.in.colony..days."] <= W & l[, "AgeAtDeath..days."] >= T), ]; w <- B[rownames(x), as.character(W)]; k <- !is.na(w)
  x <- x[k, ]; w <- w[k]; y <- x[, "AgeAtDeath..days."]; s <- as.character(x[, "StrainName"]); d <- x[, "Diet"]; idx <- split(seq_along(s), s)
  f <- function(fun, i, lab){ r <- sapply(c("CD", "HF"), function(g){ j <- d[i] == g; fun(w[i][j], y[i][j], lab[j]) }); c(r, r[2] - r[1]) }
  for(nm in c("pooled", "within")){ fun <- get(nm); e <- f(fun, seq_along(s), s)
    b <- replicate(2000, { kk <- sample(names(idx), replace = TRUE); i <- unlist(idx[kk], use.names = FALSE); f(fun, i, as.character(rep(seq_along(kk), lengths(idx[kk]))))[3] })
    out <- rbind(out, data.frame(W = W, T = T, kind = nm, mice.CD = sum(d == "CD"), mice.HF = sum(d == "HF"), CD = round(e[1], 2), HF = round(e[2], 2), diff = round(e[3], 2), lo = round(quantile(b, .025), 2), hi = round(quantile(b, .975), 2), row.names = NULL)) }
}
options(width = 140); print(out, row.names = FALSE)
write.table(out, "output/weight_lifespan_diet.txt", sep = "\t", quote = FALSE, row.names = FALSE)
