#
# qtl2_strainF1_perm.R - Genome-wide thresholds for the maps of strain and F1 means together (qtl2_strainF1.R), additive and dominance,
# by permuting strains and F1 types
#
# Rob Williams's fork. First written 23 November 2026 (Rob: "Yes, run both"); rewritten 26 November 2026 at Danny Arends's request (one
# path, no command-line input, descriptive names). Exploratory. Written by Claude (Opus 5.5). Run from the fork's root:
# Rscript qtl2_strainF1_perm.R
#
# The model is that of qtl2_strainF1.R (stated in strainF1_functions.R):  unit mean lifespan ~ F1 class + genotype + LOCO kinship.
# The permutation: each strain takes the whole row of means (every T-age) of another strain, and each F1 type that of another F1 type of
# its class (short x short, long x long). So the nesting of the survivorships and the class covariate are kept; units keep their
# genotypes and kinship; only the link between genotype and lifespan is broken. Each permuted map gives its largest -log10 P over
# markers at each T-age, for the additive and the dominance map.
# Stable range of T-ages, by the rule of the strain-means maps: up to the last T-age before the 5% single-T-age threshold first moves
# more than 0.3 from its median over T20 to T500. Permutation 0 is the real data.
# Writes output/qtl2_strainF1_perm_maxima.txt, output/qtl2_strainF1_perm_thresholds.txt and output/qtl2_strainF1_perm_run_info.txt
#

# ---- 1. Settings ---------------------------------------------------------------------------------------------------------------------
F1.FILE       <- "../bxd_mortality_rates_2026/data/corrected/AgingBXD_F1_24Nov2026_corrected.txt"
T.AGES        <- seq(20, 770, 30)
MIN.MICE      <- 3
SHORT.PARENTS <- c("BXD9", "BXD170", "BXD172", "BXD197")
LONG.PARENTS  <- c("BXD65", "BXD73b", "BXD77", "BXD102")
N.PERM        <- 1000
CORES         <- 4
SEED          <- 20261123                  # set once; the high-fat permutations continue the random stream of the chow ones

source("shared.R"); source("qtl2_cross.R"); source("strainF1_functions.R")
library(parallel)
t0 <- Sys.time()

# ---- 2. Units, genotypes and kinship, as in qtl2_strainF1.R --------------------------------------------------------------------------
f1.mice <- readF1Mice(F1.FILE)
f1.types <- f1TypeTable(f1.mice, strains, SHORT.PARENTS, LONG.PARENTS)
genotypes <- unitGenotypes(pr, strains, f1.types)
kinship <- locoKinship(genotypes$allele.dosage)
all.units <- c(strains, rownames(f1.types))

# ---- 3. Real and permuted maps, diet by diet ---------------------------------------------------------------------------------------
set.seed(SEED)
maxima <- NULL                                                              # data frame: diet, perm, map, then one column per T-age
for(diet in c("CD", "HF")){
  unit.means.all <- unitMeansAllT(ld, f1.mice, diet, T.AGES, all.units, MIN.MICE)          # units x T-ages
  permutations <- lapply(1:N.PERM, function(p) drawPermutation(rownames(unit.means.all), f1.types))   # list of row orders
  per.map <- mclapply(0:N.PERM, function(p){
    means <- unit.means.all
    if(p > 0){ means <- unit.means.all[permutations[[p]], , drop = FALSE]; rownames(means) <- rownames(unit.means.all) }  # rows move, unit names stay
    permutationMaxima(means, genotypes, kinship, f1.types)
  }, mc.cores = CORES, mc.preschedule = FALSE)
  for(p in 0:N.PERM) maxima <- rbind(maxima, data.frame(diet = diet, perm = p, map = c("add", "dom"), per.map[[p + 1]], check.names = FALSE, row.names = NULL))
  cat(diet, ":", N.PERM, "permutations done\n")
}
write.table(maxima, "output/qtl2_strainF1_perm_maxima.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# ---- 4. Thresholds ------------------------------------------------------------------------------------------------------------------
levels <- c(`5%` = 0.95, `10%` = 0.90, `suggestive (P = 0.63)` = 0.37)
T.columns <- paste0("T", T.AGES)
sink("output/qtl2_strainF1_perm_thresholds.txt", split = TRUE)
cat("Genome-wide thresholds (-log10 P) for the maps of strain and F1 means together (qtl2_strainF1.R), from", N.PERM, "permutations\n")
cat("(strains among strains, F1 types within class). Per scan: one T-age, the median over the stable range. Any T-age: any\n")
cat("marker at any T-age in the stable range.\n")
for(diet in c("CD", "HF")) for(m in c("add", "dom")){
  permuted <- as.matrix(maxima[maxima$diet == diet & maxima$map == m & maxima$perm > 0, T.columns])    # permutations x T-ages
  real <- unlist(maxima[maxima$diet == diet & maxima$map == m & maxima$perm == 0, T.columns])
  threshold.5pct <- apply(permuted, 2, quantile, 0.95)                                             # per T-age
  baseline <- median(threshold.5pct[T.AGES <= 500]); drifted <- which(abs(threshold.5pct - baseline) > 0.3 & T.AGES > 500)
  last.stable <- if(length(drifted) == 0) max(T.AGES) else T.AGES[min(drifted) - 1]
  stable <- paste0("T", T.AGES[T.AGES <= last.stable])
  cat(sprintf("\n%s, %s: stable range T20 to T%d\n", diet, c(add = "additive", dom = "dominance")[m], last.stable))
  for(l in names(levels)){
    per.T <- apply(permuted[, stable, drop = FALSE], 2, quantile, levels[l])
    cat(sprintf("  %-22s one T-age %.2f (range %.2f to %.2f over T-ages); any T-age %.2f\n", l, median(per.T), min(per.T), max(per.T),
                quantile(apply(permuted[, stable, drop = FALSE], 1, max), levels[l], names = FALSE)))
  }
  cat("  5% per scan by T-age:", paste(sprintf("%s %.2f", names(threshold.5pct), threshold.5pct), collapse = ", "), "\n")
  cat("  Real data, largest -log10 P by T-age:", paste(sprintf("%s %.2f", names(real), real), collapse = ", "), "\n")
}
sink()

writeRunInfo("output/qtl2_strainF1_perm_run_info.txt", list(F1.FILE = F1.FILE, T.AGES = range(T.AGES), MIN.MICE = MIN.MICE, N.PERM = N.PERM, CORES = CORES, SEED = SEED),
             inputs = c(F1.FILE, "data/BxD_Lifespan.txt", "data/BXD_current_rev050423.geno"), code = c("qtl2_strainF1_perm.R", "strainF1_functions.R", "qtl2_cross.R", "shared.R"),
             minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
