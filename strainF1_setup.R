#
# strainF1_setup.R - BXD strains and BXD x BXD F1 types as one set of mapping units, as in qtl2_strainF1.R: genotype and
# allele probabilities, LOCO kinship, and the unit means at every T-age
#
# Rob Williams's fork, 23 November 2026. Used by qtl2_strainF1_perm.R and chr1_twolocus.R (and chr1_f1_dropone.R). Needs
# shared.R and qtl2_cross.R sourced first. Same rules as qtl2_strainF1.R:
#   - F1 types from Suheeta Roy's F1 workbook (corrected copy in bxd_mortality_rates_2026), natural deaths on CD and HF;
#   - each F1 takes its parents' genotype probabilities: P(BB) = pB1 pB2, P(BD) = pB1 pD2 + pD1 pB2, P(DD) = pD1 pD2;
#   - kinship is allele sharing, LOCO, over strains and F1 types together;
#   - a unit enters a T-age with at least 3 mice alive at T, and its trait is their mean age at death;
#   - F1 class (inbred, short x short, long x long) is a covariate.
# Makes: f1, types, parents, f1class, units, G3 (BB, BD, DD), AP (allele probabilities), Kloco, unitMeans(diet, tps)
# and classCovar(units)
#

f1file <- if(exists("f1file")) f1file else "../bxd_mortality_rates_2026/data/corrected/AgingBXD_F1_24Nov2026_corrected.txt"
short <- c("BXD9", "BXD170", "BXD172", "BXD197"); long <- c("BXD65", "BXD73b", "BXD77", "BXD102")

f1 <- read.delim(f1file, check.names = FALSE)
f1 <- f1[f1[, "CauseOfDeath"] == "Natural" & f1[, "DietCode"] %in% c("CD", "HF"), ]
types <- sort(unique(f1[, "StrainName"]))
parents <- t(sapply(types, function(x) c(sub("x.*", "", x), paste0("BXD", sub(".*x", "", x)))))
stopifnot(all(parents %in% strains))
f1class <- ifelse(parents[, 1] %in% short & parents[, 2] %in% short, "short", ifelse(parents[, 1] %in% long & parents[, 2] %in% long, "long", NA))
names(f1class) <- types
stopifnot(!any(is.na(f1class)))

# A "unit" is one row of the maps: a BXD strain or an F1 type. Everything below is indexed by unit name, never by position.
units <- c(strains, types)
G3 <- pr; AP <- pr
for(chr in names(pr)){
  p <- pr[[chr]]; pD <- p[, "DD", ]; pB <- p[, "BB", ]; p1 <- parents[, 1]; p2 <- parents[, 2]
  g <- array(0, c(length(units), 3, dim(p)[3]), list(units, c("BB", "BD", "DD"), dimnames(p)[[3]]))
  # strains: as they are (a strain is BB or DD)
  g[strains, "BB", ] <- pB[strains, ]; g[strains, "DD", ] <- pD[strains, ]
  # an F1 takes one allele from each parent: BB if both parents give B, DD if both give D, BD otherwise. The parents are inbred, so
  # this is a product of two strains' probabilities, marker by marker (it would be an exact 0 or 1 where both calls are certain)
  g[types, "BB", ] <- pB[p1, ] * pB[p2, ]
  g[types, "BD", ] <- pB[p1, ] * pD[p2, ] + pD[p1, ] * pB[p2, ]
  g[types, "DD", ] <- pD[p1, ] * pD[p2, ]
  G3[[chr]] <- g
  # additive coding: the share of alleles that are D (strain 0 or 1; an F1 that is BD, 0.5). The coefficient is DD minus BB in days.
  a <- array(0, c(length(units), 2, dim(p)[3]), list(units, c("B", "D"), dimnames(p)[[3]]))
  a[, "B", ] <- g[, "BB", ] + g[, "BD", ] / 2; a[, "D", ] <- g[, "DD", ] + g[, "BD", ] / 2
  AP[[chr]] <- a
}
attr(G3, "crosstype") <- "f2"; attr(G3, "alleleprobs") <- FALSE
attr(AP, "alleleprobs") <- TRUE
# LOCO kinship of every pair of units: alleles shared, summed over all markers except the chromosome scanned, divided by their number
S <- lapply(AP, function(a){ m <- matrix(a, nrow = dim(a)[1]); tcrossprod(m) })      # allele sharing summed over a chromosome's markers
nm <- sapply(AP, function(a) dim(a)[3])
Kloco <- lapply(names(AP), function(chr){ k <- (Reduce(`+`, S) - S[[chr]]) / (sum(nm) - nm[chr]); dimnames(k) <- list(units, units); k })
names(Kloco) <- names(AP)

# The trait: mean age at death of a unit's mice alive at T; NA with fewer than 3. F1 class (below) is a covariate so that the average
# hybrid vigor of short x short or long x long F1s is not read as an effect of every marker.
# Unit means (units x T-ages, NA where a unit has fewer than 3 mice alive at T); strains = FALSE leaves the F1s out
unitMeans <- function(diet, tps, strains.only = FALSE){
  x <- ld[which(ld[, "Diet"] == diet), ]; z <- f1[f1[, "DietCode"] == diet, ]
  Y <- sapply(tps, function(T){
    a <- x[x[, "AgeAtDeath..days."] >= T, ]; b <- z[z[, "AgeAtDeath (days)"] >= T, ]
    n <- c(table(a[, "StrainName"]), table(b[, "StrainName"]))
    m <- c(tapply(a[, "AgeAtDeath..days."], a[, "StrainName"], mean), tapply(b[, "AgeAtDeath (days)"], b[, "StrainName"], mean))
    v <- rep(NA, length(units)); names(v) <- units; ok <- names(n)[n >= 3]; v[ok] <- m[ok]; v })
  dimnames(Y) <- list(units, paste0("T", tps))
  if(strains.only) Y <- Y[strains, , drop = FALSE]
  Y[rowSums(!is.na(Y)) > 0, , drop = FALSE]
}
classCovar <- function(u){
  cls <- ifelse(u %in% types, f1class[u], "inbred")
  X <- cbind(short = as.numeric(cls == "short"), long = as.numeric(cls == "long")); rownames(X) <- u
  X[, colSums(X) > 0, drop = FALSE]
}
