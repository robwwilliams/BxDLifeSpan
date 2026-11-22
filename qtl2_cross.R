#
# qtl2_cross.R - The BXD strains of shared.R as an R/qtl2 cross, with genotype probabilities and LOCO kinship
#
# Rob Williams's fork, 21 November 2026. Used by qtl2_individuals.R and qtl2_strainmeans.R. The strains are those of
# ldata on CD and HF that have genotypes, as a sib-mated RI panel (risib); genotypes from data/BXD_current_rev050423.geno,
# H and U as missing; the X treated like an autosome, as in Danny's maps. Needs shared.R sourced first.
# Makes: strains, cross, pr (genotype probabilities, one row per strain) and K (LOCO kinship, allele sharing)
#

library(qtl2)

qtl2.cores <- max(1, parallel::detectCores() / 2 - 1)              # as in Danny's scripts
ld <- ldata[which(ldata[, "Diet"] %in% c("CD", "HF") & grepl("BXD", ldata[, "StrainName"]) & ldata[, "StrainName"] %in% colnames(geno)),]
strains <- sort(unique(ld[, "StrainName"]))

dir <- file.path(tempdir(), "bxd_qtl2")
dir.create(dir, showWarnings = FALSE)
G <- t(as.matrix(geno[, strains]))
G[G == "H"] <- NA
write.csv(data.frame(id = strains, G, check.names = FALSE), file.path(dir, "geno.csv"), row.names = FALSE, quote = FALSE, na = "-")
write.csv(data.frame(marker = rownames(map), chr = map[, "Chr"], pos = map[, "cM"]), file.path(dir, "gmap.csv"), row.names = FALSE, quote = FALSE)
write.csv(data.frame(id = strains, dummy = 1), file.path(dir, "pheno.csv"), row.names = FALSE, quote = FALSE)
write_control_file(file.path(dir, "control.yaml"), crosstype = "risib", geno_file = "geno.csv", gmap_file = "gmap.csv", pheno_file = "pheno.csv",
                   alleles = c("B", "D"), geno_codes = c(B = 1, D = 2), na.strings = c("-", "NA"), overwrite = TRUE)
cross <- suppressWarnings(read_cross2(file.path(dir, "control.yaml")))             # warns that no chromosome is the X
pr <- calc_genoprob(cross, error_prob = 0.002, cores = qtl2.cores)
K <- calc_kinship(pr, "loco", cores = qtl2.cores)

# Signed -log10 P from a LOD on 1 d.f. and the sign of an effect
signedP <- function(lod, effect){ sign(effect) * -pchisq(2 * log(10) * pmax(lod, 0), 1, lower.tail = FALSE, log.p = TRUE) / log(10) }
