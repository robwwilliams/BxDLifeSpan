#
# strainF1_tests.R - checks of strainF1_functions.R: the genotypes and kinship of the units, and a sample of scans and permutations
# against the committed outputs of the first versions of qtl2_strainF1.R and qtl2_strainF1_perm.R
#
# Rob Williams's fork, 26 November 2026. Written by Claude (Opus 5.5). Run from the fork's root: Rscript strainF1_tests.R
# A failed check stops the script. Writes output/strainF1_tests.txt
#

F1.FILE <- "../bxd_mortality_rates_2026/data/corrected/AgingBXD_F1_24Nov2026_corrected.txt"
T.AGES <- seq(20, 770, 30); MIN.MICE <- 3
SHORT.PARENTS <- c("BXD9", "BXD170", "BXD172", "BXD197"); LONG.PARENTS <- c("BXD65", "BXD73b", "BXD77", "BXD102")
source("shared.R"); source("qtl2_cross.R"); source("strainF1_functions.R")
sink("output/strainF1_tests.txt", split = TRUE)
check <- function(what, ok){ cat(sprintf("%-4s %s\n", if(isTRUE(ok)) "ok" else "FAIL", what)); if(!isTRUE(ok)) stop("check failed: ", what) }

f1.mice <- readF1Mice(F1.FILE); f1.types <- f1TypeTable(f1.mice, strains, SHORT.PARENTS, LONG.PARENTS)
genotypes <- unitGenotypes(pr, strains, f1.types); kinship <- locoKinship(genotypes$allele.dosage)
check(sprintf("%d F1 types (%d short x short, %d long x long), %d F1 mice", nrow(f1.types), sum(f1.types$class == "short"), sum(f1.types$class == "long"), nrow(f1.mice)), nrow(f1.types) == 13)

# ---- 1. Genotypes of the units --------------------------------------------------------------------------------------------------------
three <- genotypes$three.genotypes
check("the three genotype probabilities of every unit sum to 1 at every marker", max(sapply(three, function(g) max(abs(apply(g, c(1, 3), sum) - 1)))) < 1e-9)
check("a strain is never BD", max(sapply(three, function(g) max(g[strains, "BD", ]))) == 0)
chr1 <- three[["1"]]; calls <- t(as.matrix(geno[dimnames(chr1)[[3]], strains]))                        # strains x Chr 1 markers
both.called <- function(t) calls[f1.types[t, "parent1"], ] %in% c("B", "D") & calls[f1.types[t, "parent2"], ] %in% c("B", "D")
check("on Chr 1, an F1 is BD (P > 0.99) where its parents' calls differ, and not BD (P < 0.01) where they agree",
      all(sapply(rownames(f1.types), function(t){ k <- both.called(t); differ <- calls[f1.types[t, "parent1"], k] != calls[f1.types[t, "parent2"], k]
        all(chr1[t, "BD", k][differ] > 0.99) && all(chr1[t, "BD", k][!differ] < 0.01) })))
check("allele dosage of D = P(DD) + P(BD) / 2", max(abs(genotypes$allele.dosage[["1"]][, "D", ] - (chr1[, "DD", ] + chr1[, "BD", ] / 2))) < 1e-12)
check("LOCO kinship is symmetric, and leaving out Chr 1 changes it", isSymmetric(kinship[["1"]]) && max(abs(kinship[["1"]] - kinship[["2"]])) > 0)

# ---- 2. Units per scan, against the committed counts ---------------------------------------------------------------------------------
committed.N <- read.delim("output/qtl2_meansF1_N.txt")
new.N <- t(sapply(seq_len(nrow(committed.N)), function(i){ u <- names(unitMeansAt(ld, f1.mice, committed.N$diet[i], committed.N$T[i], MIN.MICE)); c(sum(!u %in% rownames(f1.types)), sum(u %in% rownames(f1.types))) }))
check(sprintf("strains and F1 types per scan equal the committed qtl2_meansF1_N.txt in all %d scans", nrow(committed.N)), all(new.N == as.matrix(committed.N[, c("strains", "F1")])))

# ---- 3. A sample of scans, against the committed maps ---------------------------------------------------------------------------------
for(diet in c("CD", "HF")) for(Tage in c(20, 380, 770)){
  result <- scanStrainsAndF1s(unitMeansAt(ld, f1.mice, diet, Tage, MIN.MICE), genotypes, kinship, f1.types)
  diffs <- sapply(c(add = "additive", dom = "dominance"), function(m){
    committed <- read.table(paste0("output/qtl2_", diet, "_meansF1_", names(which(c(add = "additive", dom = "dominance") == m)), "_Progressive.txt"), sep = "\t", check.names = FALSE)
    max(abs(committed[, as.character(Tage)] - result[[m]]), na.rm = TRUE) })
  check(sprintf("%s T%d: additive and dominance maps equal the committed ones (largest differences %.1e, %.1e)", diet, Tage, diffs[1], diffs[2]), max(diffs) < 1e-6)
}

# ---- 4. A sample of permutations, against the committed maxima ------------------------------------------------------------------------
committed.max <- read.delim("output/qtl2_strainF1_perm_maxima.txt", check.names = FALSE)
set.seed(20261123); all.units <- c(strains, rownames(f1.types))
for(diet in c("CD", "HF")){
  unit.means.all <- unitMeansAllT(ld, f1.mice, diet, T.AGES, all.units, MIN.MICE)
  permutations <- lapply(1:1000, function(p) drawPermutation(rownames(unit.means.all), f1.types))   # the whole draw, so HF continues CD's stream
  for(p in 0:2){
    means <- unit.means.all; if(p > 0){ means <- unit.means.all[permutations[[p]], , drop = FALSE]; rownames(means) <- rownames(unit.means.all) }
    new <- permutationMaxima(means, genotypes, kinship, f1.types)
    old <- as.matrix(committed.max[committed.max$diet == diet & committed.max$perm == p, paste0("T", T.AGES)])
    check(sprintf("%s permutation %d: the largest -log10 P at every T-age equals the committed run (largest difference %.1e)", diet, p, max(abs(new - old))), max(abs(new - old)) < 1e-6)
  }
}
sink()
