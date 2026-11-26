#
# effective_tests_simpleM.R - The effective number of independent tests in the BXD genotypes used by the CTL scans (simpleM)
#
# Rob Williams's fork, 26 November 2026. Exploratory, asked for by Danny Arends and Rob ("run a simple M scan across the genotypes of the
# strains used in the analysis"). Written by Claude (Opus 5.5). Run from the fork's root after ctl_inputs.R: Rscript effective_tests_simpleM.R
#
# simpleM (Gao, Starmer & Martin 2008, Genet Epidemiol 32:361): the eigenvalues of the correlation matrix of the markers; the effective
# number of tests is the number of principal components needed to explain 99.5% of the variance. Computed
#   - by chromosome, and summed (the usual way: unlinked chromosomes are independent), and
#   - for the whole genome at once.
# Also Li & Ji (2005, Heredity 95:221), from the same eigenvalues: sum of I(lambda >= 1) + (lambda - floor(lambda)).
# The strains: those of the CTL scans with mice on chow or high fat. Genotypes: B = 0, D = 1; H and unknown calls take the marker's mean.
# A limit to keep in mind: with n strains a correlation matrix has at most n - 1 nonzero eigenvalues, so no block can count more than
# n - 1 tests, however many markers it has.
# Writes output/effective_tests_simpleM.txt
#

INPUT.DIR <- "output/ctl_scan"
VARIANCE.EXPLAINED <- 0.995                                                  # Gao et al. 2008

genotypes <- read.delim(file.path(INPUT.DIR, "genotypes.txt"), check.names = FALSE, colClasses = "character")
mice <- read.delim(file.path(INPUT.DIR, "mice.txt"), check.names = FALSE)
strains.used <- intersect(colnames(genotypes)[-(1:3)], unique(mice$strain))

dosage <- matrix(match(as.matrix(genotypes[, strains.used]), c("B", "D")) - 1, nrow(genotypes))   # markers x strains: 0 B, 1 D, NA other
dosage <- t(apply(dosage, 1, function(m){ m[is.na(m)] <- mean(m, na.rm = TRUE); m }))                # unknown calls take the marker mean
polymorphic <- apply(dosage, 1, function(m) !any(is.na(m)) && sd(m) > 0)

# Eigenvalues of the markers' correlation matrix, through the strains x strains matrix (same nonzero eigenvalues, much smaller)
markerEigenvalues <- function(markers.by.strains){
  standardized <- t(scale(t(markers.by.strains)))                           # each marker: mean 0, sd 1 over strains
  gram <- crossprod(standardized) / (ncol(markers.by.strains) - 1)          # strains x strains
  pmax(eigen(gram, symmetric = TRUE, only.values = TRUE)$values, 0)         # sum = number of markers
}
simpleM <- function(lambda) which(cumsum(lambda) / sum(lambda) >= VARIANCE.EXPLAINED)[1]
liJi <- function(lambda) sum((lambda >= 1) + (lambda - floor(lambda)))

per.chr <- do.call(rbind, lapply(unique(genotypes$Chr), function(chr){
  k <- which(genotypes$Chr == chr & polymorphic); lambda <- markerEigenvalues(dosage[k, , drop = FALSE])
  data.frame(Chr = chr, markers = length(k), simpleM = simpleM(lambda), LiJi = round(liJi(lambda), 1))
}))
genome.lambda <- markerEigenvalues(dosage[polymorphic, ])

sink("output/effective_tests_simpleM.txt", split = TRUE)
cat("Effective number of independent tests in the genotypes of the", length(strains.used), "BXD strains of the CTL scans (simpleM, 99.5%)\n\n")
print(per.chr, row.names = FALSE)
cat(sprintf("\nSum over chromosomes: simpleM %d, Li & Ji %.0f (of %d polymorphic markers)\n", sum(per.chr$simpleM), sum(per.chr$LiJi), sum(polymorphic)))
cat(sprintf("Whole genome at once: simpleM %d, Li & Ji %.0f (at most %d: strains - 1)\n", simpleM(genome.lambda), liJi(genome.lambda), length(strains.used) - 1))
for(m in c(sum(per.chr$simpleM), simpleM(genome.lambda)))
  cat(sprintf("Bonferroni 5%% for %d tests: -log10 P %.2f (one T-age)\n", m, -log10(0.05 / m)))
sink()
