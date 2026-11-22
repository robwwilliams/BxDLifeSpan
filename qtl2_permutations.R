#
# qtl2_permutations.R - Genome-wide thresholds for the strain-means actuarial maps, by permuting strains
#
# Rob Williams's fork, 21 November 2026. Exploratory. Rob: "We should add some guidance on genome-wide significance in
# either the figure or the legend." The strain is the unit of inference, so the permutation is over strains.
#   - The maps are those of qtl2_strainmeans.R with equal weights: at each T-age (20 to 770 days, 30-day steps), the mean
#     age at death of each BXD strain with at least 3 mice alive at T, R/qtl2 scan1 with LOCO kinship (REML).
#   - Each permutation gives every strain the whole set of strain means (all T-ages) of another strain, so the nesting of
#     the survivorships is kept, and only the link between genotypes and lifespans is broken. Strains keep their
#     genotypes and kinship.
#   - Each permutation records the largest -log10 P over markers at each T-age, and over markers and T-ages together.
#     Two thresholds follow: per scan (one T-age) and actuarial (any marker at any T-age).
#   - The permutations are drawn in advance from one seed.
# Writes output/qtl2_permutation_thresholds.txt and output/qtl2_permutation_maxima.txt
#   Rscript qtl2_permutations.R [permutations, default 1000] [cores, default 4]
#

source("shared.R")
source("qtl2_cross.R")                                                              # strains, cross, pr, K
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
nperm <- if(length(args) >= 1) as.numeric(args[1]) else 1000
cores <- if(length(args) >= 2) as.numeric(args[2]) else 4
tps <- seq(20, 770, 30)
lp <- function(lod) -pchisq(2 * log(10) * pmax(lod, 0), 1, lower.tail = FALSE, log.p = TRUE) / log(10)

set.seed(20261122)
res <- NULL
for(diet in c("CD", "HF")){
  x <- ld[which(ld[, "Diet"] == diet),]
  Y <- sapply(tps, function(T){ a <- x[x[, "AgeAtDeath..days."] >= T, ]; n <- table(a[, "StrainName"]); m <- tapply(a[, "AgeAtDeath..days."], a[, "StrainName"], mean)
                                 v <- rep(NA, length(strains)); names(v) <- strains; ok <- names(n)[n >= 3]; v[ok] <- m[ok]; v })
  dimnames(Y) <- list(strains, paste0("T", tps))
  Y <- Y[rowSums(!is.na(Y)) > 0, , drop = FALSE]
  s <- rownames(Y)
  prs <- pr[s, ]; Ks <- lapply(K, function(k) k[s, s])
  perms <- lapply(1:nperm, function(p) sample.int(length(s)))
  m <- do.call(rbind, mclapply(1:nperm, function(p){
    Yp <- Y[perms[[p]], , drop = FALSE]; rownames(Yp) <- s
    out <- scan1(prs, Yp, kinship = Ks)
    c(perm = p, apply(lp(out), 2, max, na.rm = TRUE))
  }, mc.cores = cores))
  res <- rbind(res, data.frame(diet = diet, m))
  cat(diet, ":", nperm, "permutations done\n")
}
write.table(res, "output/qtl2_permutation_maxima.txt", sep = "\t", quote = FALSE, row.names = FALSE)
sink("output/qtl2_permutation_thresholds.txt", split = TRUE)
cat("Genome-wide thresholds (-log10 P) for the strain-means maps (equal weights, LOCO kinship), from", nperm, "permutations of strains\n")
for(diet in c("CD", "HF")){
  r <- as.matrix(res[res[, "diet"] == diet, -(1:2)])
  ps <- apply(r, 2, quantile, c(0.90, 0.95, 0.99))
  act <- quantile(apply(r, 1, max), c(0.90, 0.95, 0.99))
  cat(sprintf("\n%s per scan (one T-age): 10%% %.2f, 5%% %.2f, 1%% %.2f (medians over T-ages); 5%% ranges %.2f to %.2f over T-ages\n", diet,
              median(ps[1, ]), median(ps[2, ]), median(ps[3, ]), min(ps[2, ]), max(ps[2, ])))
  cat(sprintf("%s actuarial (any marker at any T-age): 10%% %.2f, 5%% %.2f, 1%% %.2f\n", diet, act[1], act[2], act[3]))
  cat("  5% per scan by T-age:", paste(sprintf("%s %.2f", colnames(ps), ps[2, ]), collapse = ", "), "\n")
}
sink()
