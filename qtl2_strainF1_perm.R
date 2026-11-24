#
# qtl2_strainF1_perm.R - Genome-wide thresholds for the maps of strain and F1 means together (qtl2_strainF1.R), additive
# and dominance, by permuting strains and F1 types
#
# Rob Williams's fork, 23 November 2026. Exploratory. Asked for by Rob ("Yes, run both", 23 November 2026, to: "genome-wide
# thresholds for the strain + F1 maps, by permuting strains and F1 types").
#   - The maps are those of qtl2_strainF1.R: at each T-age (20 to 770 days, 30-day steps), the mean age at death of each
#     BXD strain and each F1 type with at least 3 mice alive at T, equal weights, F1 class as a covariate, LOCO kinship
#     over strains and F1s. Additive: allele probabilities (1 d.f.). Dominance: the three genotypes against the additive
#     model (1 d.f.).
#   - Each permutation gives every strain the whole set of means (all T-ages) of another strain, as in
#     qtl2_permutations.R, and every F1 type those of another F1 type of its class (short x short or long x long). So the
#     nesting of the survivorships and the class covariate are kept, and only the link between genotypes and lifespans
#     is broken. Units keep their genotypes and kinship.
#   - Each permutation records the largest -log10 P over markers at each T-age, for each map.
#   - Stable range of T-ages, by the rule of the strain-means maps: up to the last T-age before the 5% single-T-age
#     threshold first moves more than 0.3 from its median over T20 to T500.
#   - Permutation 0 is the real data, as a check against qtl2_strainF1.R's results.
# Writes output/qtl2_strainF1_perm_maxima.txt and output/qtl2_strainF1_perm_thresholds.txt
#   Rscript qtl2_strainF1_perm.R [permutations, default 1000] [cores, default 4]
#

source("shared.R")
source("qtl2_cross.R")                                                              # strains, pr
source("strainF1_setup.R")                                                         # units, AP, G3, Kloco, unitMeans, classCovar
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
nperm <- if(length(args) >= 1) as.numeric(args[1]) else 1000
cores <- if(length(args) >= 2) as.numeric(args[2]) else 4
tps <- seq(20, 770, 30)
lp <- function(lod) -pchisq(2 * log(10) * pmax(lod, 0), 1, lower.tail = FALSE, log.p = TRUE) / log(10)

set.seed(20261123)
res <- NULL
for(diet in c("CD", "HF")){
  Y <- unitMeans(diet, tps); u <- rownames(Y)
  X <- classCovar(u); ap <- AP[u, ]; g3 <- G3[u, ]; Ku <- lapply(Kloco, function(k) k[u, u])
  groups <- list(which(u %in% strains), which(u %in% types[f1class == "short"]), which(u %in% types[f1class == "long"]))
  perms <- lapply(1:nperm, function(p){ o <- seq_along(u); for(g in groups) o[g] <- g[sample.int(length(g))]; o })
  m <- do.call(rbind, mclapply(0:nperm, function(p){
    Yp <- if(p == 0) Y else Y[perms[[p]], , drop = FALSE]; rownames(Yp) <- u
    add <- scan1(ap, Yp, kinship = Ku, addcovar = X)
    full <- scan1(g3, Yp, kinship = Ku, addcovar = X)
    rbind(c(perm = p, map = 1, apply(lp(add), 2, max, na.rm = TRUE)),
          c(perm = p, map = 2, apply(lp(full - add), 2, max, na.rm = TRUE)))
  }, mc.cores = cores, mc.preschedule = FALSE))
  res <- rbind(res, data.frame(diet = diet, perm = m[, "perm"], map = c("add", "dom")[m[, "map"]], m[, -(1:2)], check.names = FALSE))
  cat(diet, ":", nperm, "permutations done\n")
}
write.table(res, "output/qtl2_strainF1_perm_maxima.txt", sep = "\t", quote = FALSE, row.names = FALSE)

### Thresholds
levels <- c(`5%` = 0.95, `10%` = 0.90, `suggestive (P = 0.63)` = 0.37)
sink("output/qtl2_strainF1_perm_thresholds.txt", split = TRUE)
cat("Genome-wide thresholds (-log10 P) for the maps of strain and F1 means together (qtl2_strainF1.R), from", nperm, "permutations\n")
cat("(strains among strains, F1 types within class). Per scan: one T-age, the median over the stable range. Any T-age: any\n")
cat("marker at any T-age in the stable range.\n")
for(diet in c("CD", "HF")) for(mp in c("add", "dom")){
  r <- as.matrix(res[res$diet == diet & res$map == mp & res$perm > 0, paste0("T", tps)])
  obs <- unlist(res[res$diet == diet & res$map == mp & res$perm == 0, paste0("T", tps)])
  q5 <- apply(r, 2, quantile, 0.95)
  base <- median(q5[tps <= 500]); off <- which(abs(q5 - base) > 0.3 & tps > 500)
  last <- if(length(off) == 0) max(tps) else tps[min(off) - 1]
  keep <- paste0("T", tps[tps <= last])
  cat(sprintf("\n%s, %s: stable range T20 to T%d\n", diet, c(add = "additive", dom = "dominance")[mp], last))
  for(l in names(levels)){
    one <- apply(r[, keep, drop = FALSE], 2, quantile, levels[l])
    cat(sprintf("  %-22s one T-age %.2f (range %.2f to %.2f over T-ages); any T-age %.2f\n", l, median(one), min(one), max(one),
                quantile(apply(r[, keep, drop = FALSE], 1, max), levels[l], names = FALSE)))
  }
  cat("  5% per scan by T-age:", paste(sprintf("%s %.2f", names(q5), q5), collapse = ", "), "\n")
  cat("  Real data, largest -log10 P by T-age:", paste(sprintf("%s %.2f", names(obs), obs), collapse = ", "), "\n")
}
sink()
