#
# qtl2_strainF1.R - Actuarial maps of BXD strain means and BXD x BXD F1 means together: additive and dominance effects
#
# Rob Williams's fork. First written 21 November 2026; rewritten 26 November 2026 at Danny Arends's request (one path, no command-line
# input, the model stated, descriptive names). Exploratory. Written by Claude (Opus 5.5). Run from the fork's root: Rscript qtl2_strainF1.R
#
# THE MODEL (details in strainF1_functions.R), at every marker, for each diet and T-age (20 to 770 days in 30-day steps):
#
#     unit mean lifespan  ~  F1 class  +  genotype at the marker  +  polygenic background (LOCO kinship)
#
#   - A unit is a BXD strain or an F1 type with at least 3 mice alive at T; its trait is their mean age at death; units count equally.
#   - F1 class is two 0/1 columns, short x short and long x long, with inbred strains as the baseline. Each class gets its own shift, so
#     the strong heterosis of short x short F1s and the weak heterosis of long x long F1s are each absorbed (f1_heterosis.R: within a
#     class the F1 types differ no more than sampling error allows).
#   - Additive scan: D allele dosage (1 d.f.). Full scan: BB, BD, DD (2 d.f.). Dominance = full - additive (1 d.f.).
#
# Steps: 1 settings; 2 strains and F1 types; 3 genotypes and kinship of the units; 4 the scans; 5 the outputs; 6 the proximal Chr 1 report.
# Writes output/qtl2_<CD|HF>_meansF1_<add|dom>_Progressive.txt (markers x T-ages), output/qtl2_meansF1_N.txt (units per scan),
# output/strainF1_chr1.txt and output/qtl2_strainF1_run_info.txt
#

# ---- 1. Settings ---------------------------------------------------------------------------------------------------------------------
F1.FILE       <- "../bxd_mortality_rates_2026/data/corrected/AgingBXD_F1_24Nov2026_corrected.txt"   # Suheeta Roy's F1 workbook, corrected
T.AGES        <- seq(20, 770, 30)                                          # survivorship ages, days
MIN.MICE      <- 3                                                         # mice alive at T that a unit needs
SHORT.PARENTS <- c("BXD9", "BXD170", "BXD172", "BXD197")                   # short-lived parents of the F1s
LONG.PARENTS  <- c("BXD65", "BXD73b", "BXD77", "BXD102")                   # long-lived parents
CORES         <- 4

source("shared.R")                                                          # Danny's: ldata (BXD mice), geno, map
source("qtl2_cross.R")                                                      # strains, ld (their mice on CD and HF), pr, signedP
source("strainF1_functions.R")
library(parallel)
t0 <- Sys.time()

# ---- 2. Strains and F1 types ---------------------------------------------------------------------------------------------------------
f1.mice <- readF1Mice(F1.FILE)                                              # data frame: one row per F1 mouse
f1.types <- f1TypeTable(f1.mice, strains, SHORT.PARENTS, LONG.PARENTS)      # data frame: one row per F1 type (parent1, parent2, class)

# ---- 3. Genotypes and kinship of the units -------------------------------------------------------------------------------------------
genotypes <- unitGenotypes(pr, strains, f1.types)                          # list: three.genotypes, allele.dosage
kinship <- locoKinship(genotypes$allele.dosage)                            # list: one units x units matrix per chromosome

# ---- 4. The scans, one per diet and T-age --------------------------------------------------------------------------------------------
scan.list <- expand.grid(Tage = T.AGES, diet = c("CD", "HF"), stringsAsFactors = FALSE)   # 52 scans
results <- mclapply(seq_len(nrow(scan.list)), function(i){
  unit.means <- unitMeansAt(ld, f1.mice, scan.list$diet[i], scan.list$Tage[i], MIN.MICE)
  scanStrainsAndF1s(unit.means, genotypes, kinship, f1.types)
}, mc.cores = CORES)

# ---- 5. The outputs: markers x T-ages for each diet and map (Danny's layout) ---------------------------------------------------------
map.files <- c(additive = "add", dominance = "dom")                         # map -> file name part
for(diet in c("CD", "HF")) for(m in names(map.files)){
  signed.log10P <- sapply(which(scan.list$diet == diet), function(i) results[[i]][[m]])         # markers x T-ages
  out <- cbind(map, signed.log10P); colnames(out) <- c("Chr", "Locus", "cM", "Mb", T.AGES)
  write.table(out, paste0("output/qtl2_", diet, "_meansF1_", map.files[m], "_Progressive.txt"), sep = "\t", quote = FALSE)
}
units.per.scan <- do.call(rbind, lapply(seq_len(nrow(scan.list)), function(i) c(diet = scan.list$diet[i], T = scan.list$Tage[i], results[[i]]$n)))
write.table(units.per.scan, "output/qtl2_meansF1_N.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# ---- 6. Proximal Chr 1: where each F1 type is heterozygous, and the peaks with and without F1s ---------------------------------------
sink("output/strainF1_chr1.txt", split = TRUE)
cat("Proximal Chr 1 (0 to 60 Mb). Segments where each F1 type's parents differ (B against D calls of BXD_current_rev050423.geno):\n")
chr1.markers <- which(map[, "Chr"] == "1" & map[, "Mb"] <= 60)             # row numbers of proximal Chr 1 markers
for(type in rownames(f1.types)){
  call1 <- geno[chr1.markers, f1.types[type, "parent1"]]; call2 <- geno[chr1.markers, f1.types[type, "parent2"]]
  heterozygous <- !is.na(call1) & !is.na(call2) & call1 %in% c("B", "D") & call2 %in% c("B", "D") & call1 != call2
  runs <- rle(heterozygous); run.end <- cumsum(runs$lengths); run.start <- run.end - runs$lengths + 1
  segments <- paste(sprintf("%.1f-%.1f", map[chr1.markers[run.start[runs$values]], "Mb"], map[chr1.markers[run.end[runs$values]], "Mb"]), collapse = ", ")
  cat(sprintf("  %-12s (%-5s x %-5s, %s x %s): %s\n", type, f1.types[type, "parent1"], f1.types[type, "parent2"], f1.types[type, "class"], f1.types[type, "class"],
              if(any(heterozygous)) segments else "none"))
}
lodFromP <- function(signed.log10P) qchisq(pmin(10^-abs(signed.log10P), 1), 1, lower.tail = FALSE) / (2 * log(10))
cat("\nStrongest marker on proximal Chr 1, with its 1.5-LOD support interval at that T-age: strains only (equal weights,\n")
cat("qtl2_strainmeans.R) against strains and F1s (additive). Positive: the D allele lives longer.\n")
for(diet in c("CD", "HF")){
  strains.only <- read.table(paste0("output/qtl2_", diet, "_means_equal_Progressive.txt"), sep = "\t")
  with.f1 <- read.table(paste0("output/qtl2_", diet, "_meansF1_add_Progressive.txt"), sep = "\t")
  for(which.map in c("strains only", "strains + F1s")){
    P <- as.matrix((if(which.map == "strains only") strains.only else with.f1)[rownames(map)[chr1.markers], paste0("X", T.AGES)])   # markers x T-ages
    peak <- which(abs(P) == max(abs(P), na.rm = TRUE), arr.ind = TRUE)[1, ]                   # (marker row, T-age column)
    lod <- lodFromP(P[, peak[2]]); inside <- which(lod >= max(lod) - 1.5)
    cat(sprintf("  %s, %-14s peak %s at %.1f Mb, T%d: -log10 P %+.2f (LOD %.2f); 1.5-LOD interval %.1f-%.1f Mb (%.1f Mb)\n", diet, which.map, rownames(P)[peak[1]],
                map[chr1.markers[peak[1]], "Mb"], T.AGES[peak[2]], P[peak[1], peak[2]], lod[peak[1]], map[chr1.markers[min(inside)], "Mb"],
                map[chr1.markers[max(inside)], "Mb"], map[chr1.markers[max(inside)], "Mb"] - map[chr1.markers[min(inside)], "Mb"]))
  }
}
sink()

writeRunInfo("output/qtl2_strainF1_run_info.txt", list(F1.FILE = F1.FILE, T.AGES = range(T.AGES), MIN.MICE = MIN.MICE, CORES = CORES),
             inputs = c(F1.FILE, "data/BxD_Lifespan.txt", "data/BXD_current_rev050423.geno"), code = c("qtl2_strainF1.R", "strainF1_functions.R", "qtl2_cross.R", "shared.R"),
             minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
