#
# ctl_fisher_simulation.R - How well does Danny Arends's CTL test of two Spearman correlations work with strain numbers like the BXD's?
#
# Rob Williams's fork, 26 November 2026. Exploratory, asked for by Rob ("run simulations ... -logP is about 5 with a scenario using 86
# strains ... and with 0.9"). Written by Claude (Opus 5.5). Run from the fork's root: Rscript ctl_fisher_simulation.R
#
# Each simulated group is a set of strains, each with a mean weight and a mean lifespan drawn from a bivariate normal with a chosen
# Spearman correlation. The two groups' Spearman correlations are compared by the CTL test with two weights for Fisher's z:
#   as coded    df <- n - (length(cor) - 1), which with 2 groups is n - 1 (CTL.R; in UM-HET3, with 4 genotypes, it is n - 3)
#   standard    n - 3, with Spearman's extra variance (1.06 / (n - 3), Fieller et al. 1957)
# 1. Rob's scenario: 86 strains, a Spearman gap of 0.9 (+0.45 in one group, -0.45 in the other), split 43/43 and 36/50.
# 2. No CTL (both groups r = 0): is the test calibrated? Inflation (median chi-square / expected) and how often P < 1e-4 (expected 1e-4).
# Writes output/ctl_fisher_simulation.txt
#

N.SIM <- 20000
SEED <- 20261126

fisherZ <- function(r) 0.5 * log((1 + r) / (1 - r))
ctlLog10P <- function(rB, nB, rD, nD, weight){                              # -log10 P of the two-group heterogeneity chi-square, 1 d.f.
  wB <- weight(nB); wD <- weight(nD)
  chisq <- wB * wD / (wB + wD) * (fisherZ(rB) - fisherZ(rD))^2              # = sum(w z^2) - sum(w z)^2 / sum(w), for two groups
  -pchisq(chisq, 1, lower.tail = FALSE, log.p = TRUE) / log(10)
}
weightAsCoded <- function(n) n - 1
weightStandard <- function(n) (n - 3) / 1.06

# Strain means with a given Spearman correlation: Pearson rho = 2 sin(pi rho.s / 6) for a bivariate normal
simulateSpearman <- function(n.strains, spearman){
  rho <- 2 * sin(pi * spearman / 6)
  weight <- rnorm(n.strains); lifespan <- rho * weight + sqrt(1 - rho^2) * rnorm(n.strains)
  cor(weight, lifespan, method = "spearman")
}

set.seed(SEED)
sink("output/ctl_fisher_simulation.txt", split = TRUE)
cat("CTL test of two Spearman correlations of strain means, simulated (", N.SIM, "replicates per scenario)\n\n")

cat("1. Rob's scenario: 86 strains, true Spearman +0.45 against -0.45 (a gap of 0.9)\n")
for(split in list(c(43, 43), c(36, 50))){
  rB <- replicate(N.SIM, simulateSpearman(split[1], 0.45)); rD <- replicate(N.SIM, simulateSpearman(split[2], -0.45))
  coded <- ctlLog10P(rB, split[1], rD, split[2], weightAsCoded); standard <- ctlLog10P(rB, split[1], rD, split[2], weightStandard)
  cat(sprintf("   %d / %d strains: -log10 P at exactly r = +0.45 / -0.45: as coded %.2f, standard %.2f\n", split[1], split[2],
              ctlLog10P(0.45, split[1], -0.45, split[2], weightAsCoded), ctlLog10P(0.45, split[1], -0.45, split[2], weightStandard)))
  cat(sprintf("     over replicates, as coded: median %.2f (middle 50%%: %.2f to %.2f); standard: median %.2f (%.2f to %.2f)\n",
              median(coded), quantile(coded, .25), quantile(coded, .75), median(standard), quantile(standard, .25), quantile(standard, .75)))
}

cat("\n2. No CTL (true r = 0 in both groups): is the test calibrated? Expected: inflation 1.00, P < 1e-4 in 0.010% of replicates\n")
for(split in list(c(43, 43), c(36, 50), c(20, 20), c(19, 8), c(10, 10))){
  rB <- replicate(N.SIM, simulateSpearman(split[1], 0)); rD <- replicate(N.SIM, simulateSpearman(split[2], 0))
  for(w in list(c("as coded", "weightAsCoded"), c("standard", "weightStandard"))){
    p <- 10^-ctlLog10P(rB, split[1], rD, split[2], get(w[2]))
    cat(sprintf("   %2d / %2d strains, %-8s inflation %.2f; P < 1e-3 in %.3f%% (expected 0.100%%); P < 1e-4 in %.3f%%\n", split[1], split[2], w[1],
                median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1), 100 * mean(p < 1e-3), 100 * mean(p < 1e-4)))
  }
}
sink()
