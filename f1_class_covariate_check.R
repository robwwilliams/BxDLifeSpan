#
# f1_class_covariate_check.R - What does the F1 class covariate absorb: heterosis, or the genetic merit of the parents too?
#
# Rob Williams's fork, 26 November 2026. Exploratory, asked for by Danny Arends's question (relayed by Rob): short and long encode the
# parents, and "there's a genetic component to your parents being a long-lived or a short-lived mouse". Written by Claude (Opus 5.5).
# Run from the fork's root: Rscript f1_class_covariate_check.R. Writes output/f1_class_covariate_check.txt
#
# The null model of qtl2_strainF1.R (no marker):  unit mean lifespan ~ F1 class (short, long; inbred baseline) + polygenic (kinship)
# Its two class coefficients are compared with
#   - the same model without kinship (plain lm): the class mean against the strains' mean, which mixes heterosis with parental merit;
#   - heterosis proper: F1 mean minus the mean of its two parents' means (midparent), averaged over the F1 types of the class.
# If the kinship term took up all of the parents' genetic merit, the coefficients would equal heterosis proper.
#

F1.FILE <- "../bxd_mortality_rates_2026/data/corrected/AgingBXD_F1_24Nov2026_corrected.txt"
SHORT.PARENTS <- c("BXD9", "BXD170", "BXD172", "BXD197"); LONG.PARENTS <- c("BXD65", "BXD73b", "BXD77", "BXD102")
T.AGES <- c(20, 200, 380, 560); MIN.MICE <- 3
source("shared.R"); source("qtl2_cross.R"); source("strainF1_functions.R")

f1.mice <- readF1Mice(F1.FILE); f1.types <- f1TypeTable(f1.mice, strains, SHORT.PARENTS, LONG.PARENTS)
dosage <- unitGenotypes(pr, strains, f1.types)$allele.dosage
shared.alleles <- Reduce(`+`, lapply(dosage, function(a){ m <- matrix(a, nrow = dim(a)[1]); tcrossprod(m) }))   # units x units, all markers
kinship <- shared.alleles / sum(sapply(dosage, function(a) dim(a)[3])); dimnames(kinship) <- rep(list(dimnames(dosage[[1]])[[1]]), 2)

classShift <- function(unit.means){                                         # returns heritability and the class coefficients, with and without kinship
  units <- names(unit.means); covariate <- classCovariate(units, f1.types)
  h2 <- est_herit(matrix(unit.means, dimnames = list(units, "lifespan")), kinship[units, units], addcovar = covariate)
  V <- h2 * 2 * kinship[units, units] + (1 - h2) * diag(length(units))      # R/qtl2's covariance: h2 2K + (1 - h2) I
  design <- cbind(intercept = 1, covariate); Vinv <- solve(V)
  gls <- solve(t(design) %*% Vinv %*% design, t(design) %*% Vinv %*% unit.means)[, 1]
  c(h2 = unname(h2), gls[c("short", "long")], lm = coef(lm(unit.means ~ covariate))[c("covariateshort", "covariatelong")])
}
heterosisProper <- function(unit.means, cls){                              # mean over the class's F1 types of F1 mean minus midparent
  types <- rownames(f1.types)[f1.types$class == cls & rownames(f1.types) %in% names(unit.means)]
  mean(sapply(types, function(t) unit.means[t] - mean(unit.means[c(f1.types[t, "parent1"], f1.types[t, "parent2"])])))
}

sink("output/f1_class_covariate_check.txt", split = TRUE)
cat("F1 class coefficients (days) in the null model of qtl2_strainF1.R, with kinship (as mapped) and without (plain lm), against heterosis\n")
cat("proper (F1 mean minus midparent). Units: strains and F1 types with at least 3 mice alive at T.\n\n")
for(diet in c("CD", "HF")) for(Tage in T.AGES){
  unit.means <- unitMeansAt(ld, f1.mice, diet, Tage, MIN.MICE); s <- classShift(unit.means)
  cat(sprintf("%s T%-3d h2 %.2f | with kinship: short %+4.0f, long %+4.0f | plain lm: short %+4.0f, long %+4.0f | heterosis proper: short %+4.0f, long %+4.0f\n",
              diet, Tage, s["h2"], s["short"], s["long"], s["lm.covariateshort"], s["lm.covariatelong"], heterosisProper(unit.means, "short"), heterosisProper(unit.means, "long")))
}
sink()
