#
# squirm_inputs.R - The inputs of the individual-mouse permutations (qtl2_individuals_perm.R), written into the main
# repository so that they can run on ATHENA through Squirm
#
# Rob Williams's fork, 21 November 2026. Squirm runs a job in its own clone of robwwilliams/bxd_mortality_rates_2026 and
# reads nothing outside it, so this writes there the derived inputs the job needs, and nothing of Danny's code:
#   - bxd_risib/: the R/qtl2 cross of qtl2_cross.R (the BXD genotypes of GeneNetwork's BXD_current_rev050423.geno for the
#     strains mapped, H and U as missing; the marker map in cM), with its control file;
#   - mice.txt: one row per mouse mapped (BXD strains on CD and HF with genotypes): ear tag, strain, diet, age at death,
#     vivarium, age at set-up, and the permutation stratum of its strain (Taylor's strains or the newer ones, from the
#     epochs of GeneNetwork's BDL_10002 export, as in qtl2_individuals_perm.R).
#   Rscript squirm_inputs.R [output folder, default ../bxd_mortality_rates_2026/QTL/inputs] [BDL_10002.xlsx path]
#

source("shared.R")
source("qtl2_cross.R")                                                              # strains, ld, and the cross files in dir
library(readxl)

args <- commandArgs(trailingOnly = TRUE)
out <- if(length(args) >= 1) args[1] else "../bxd_mortality_rates_2026/QTL/inputs"
bdl <- if(length(args) >= 2) args[2] else "../bxd_mortality_rates_2026/data/original/BDL_10002.xlsx"
dir.create(file.path(out, "bxd_risib"), recursive = TRUE, showWarnings = FALSE)
for(f in c("control.yaml", "geno.csv", "gmap.csv", "pheno.csv")) file.copy(file.path(dir, f), file.path(out, "bxd_risib", f), overwrite = TRUE)

x <- as.data.frame(suppressMessages(suppressWarnings(read_excel(bdl, skip = 10))))
x <- x[!is.na(x[, "Epoch"]), ]
epoch <- tapply(x[, "Epoch"], x[, "Strain"], function(v) names(sort(table(v), decreasing = TRUE))[1])[strains]
num <- as.numeric(sub("^BXD([0-9]+).*", "\\1", strains))
stratum <- ifelse(is.na(epoch), ifelse(num <= 42, "Taylor", "newer"), ifelse(substr(epoch, 1, 1) %in% c("1", "2"), "Taylor", "newer"))
names(stratum) <- strains

mice <- data.frame(id = rownames(ld), strain = ld[, "StrainName"], diet = ld[, "Diet"], age = ld[, "AgeAtDeath..days."],
                   vivarium = ld[, "vivarium"], setup.age = as.numeric(ld[, "AgeAtSetUp.in.colony..days."]), stratum = stratum[ld[, "StrainName"]])
write.table(mice, file.path(out, "mice.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(nrow(mice), "mice of", length(strains), "strains;", paste(names(table(stratum)), table(stratum), "strains", collapse = ", "), "\n")
