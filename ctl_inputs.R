#
# ctl_inputs.R - the two input tables of ctl_scan.R, made from the BXD data with Danny Arends's own code
#
# Rob Williams's fork, 26 November 2026. Written by Claude (Opus 5.5). Run from the fork's root, with no arguments: Rscript ctl_inputs.R
#
# Reads, through Danny's shared.R and adjust.R (both unchanged): data/BxD_Lifespan.txt, data/BxD_Bodyweights.txt, data/BDL_10013.csv,
#   data/BxD_NAM_Lifespan.txt, data/BxD_NAM_Bodyweights.txt, data/BXD_current_rev050423.geno and output/bodyweightInterpolatedData.txt
# The linear adjustments are Danny's adjust.R, fitted within each diet; set-up ages below 0 are first set to NA, as in CTL.R:
#   lifespan                AgeAtDeath ~ vivarium + set-up age class (0-180, 180-360, ... days); residual + mean, rounded to a day
#   weight at each T-age    weight ~ vivarium + set-up age class, the same way (vivarium alone where one class only has weights)
# Writes
#   output/ctl_scan/mice.txt       one row per BXD mouse on chow (CD) or high fat (HF) with weights: strain, diet, set-up age, age at
#                                  death, adjusted lifespan, and adjusted weight w<age> at each T-age of CTL.R (50 to 770 days)
#   output/ctl_scan/genotypes.txt  B, D or other at each marker, for every BXD strain with lifespan data and genotypes
#

source("shared.R")
ldata[which(ldata[, "AgeAtSetUp.in.colony..days."] < 0), "AgeAtSetUp.in.colony..days."] <- NA
BWs <- read.table("output/bodyweightInterpolatedData.txt", sep = "\t", quote = "", check.names = FALSE)
invisible(capture.output(source("adjust.R")))                               # adds ldata$Yadj; makes BWsAdj and TPStodo

bxd <- unique(as.character(ldata[grepl("BXD", ldata[, "StrainName"]), "StrainName"]))
strains <- sort(bxd[bxd %in% colnames(geno)])

ids <- rownames(ldata)[ldata[, "StrainName"] %in% strains & ldata[, "Diet"] %in% c("CD", "HF")]
ids <- ids[ids %in% rownames(BWsAdj)]
mice <- data.frame(mouse = ids, strain = ldata[ids, "StrainName"], diet = ldata[ids, "Diet"], setup = ldata[ids, "AgeAtSetUp.in.colony..days."],
                   age = ldata[ids, "AgeAtDeath..days."], lifespan = ldata[ids, "Yadj"], BWsAdj[ids, TPStodo], check.names = FALSE)
colnames(mice)[-(1:6)] <- paste0("w", TPStodo)

dir.create("output/ctl_scan", showWarnings = FALSE)
write.table(mice, "output/ctl_scan/mice.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(map[, c("Locus", "Chr", "Mb")], geno[, strains], check.names = FALSE), "output/ctl_scan/genotypes.txt", sep = "\t", quote = FALSE, row.names = FALSE)
cat(nrow(mice), "mice of", length(unique(mice$strain)), "strains;", length(strains), "genotyped strains;", nrow(geno), "markers\n")
