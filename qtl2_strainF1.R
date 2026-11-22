#
# qtl2_strainF1.R - Actuarial maps of BXD strain means and BXD x BXD F1 means together, with R/qtl2 and LOCO kinship:
# additive and dominance effects
#
# Rob Williams's fork, 21 November 2026. Exploratory. Rob: "I personally would love to see the strain+F1 map using the
# 'simple method' (no weighting of strain means but n at or above 3)." The 13 F1 types were chosen as short x short and
# long x long crosses with useful recombinations over the proximal 50 Mb of Chr 1, to see dominance and to sharpen the
# positions of the Chr 1 loci. Earlier, of F1 type: "the two types of F1 were by design, so just add the 'type' as a
# cofactor", and "They are genetically just like an F2 cross."
#   - Units: at each diet and T-age (20 to 770 days, 30-day steps), the BXD strains of shared.R and the F1 types of Suheeta
#     Roy's F1 workbook (data/corrected/AgingBXD_F1_BodyWeight_NaturalDeath_corrected.txt in bxd_mortality_rates_2026,
#     euthanized F1s counted as natural deaths), each with at least 3 mice alive at T. The trait is the unit's mean age
#     at death of those mice. Every unit is weighted equally. B6D2F1, D2B6F1 and the parental strains are left out, as in
#     the other maps.
#   - Genotypes: the strains' probabilities from qtl2_cross.R. Each F1 takes its parents': P(BB) = pB1 pB2,
#     P(BD) = pB1 pD2 + pD1 pB2, P(DD) = pD1 pD2. Kinship is allele sharing, LOCO, over strains and F1s together, so an
#     F1 shares half its alleles with itself at markers where its parents differ.
#   - Cofactor: F1 type (inbred strain, short x short, long x long), so that the average hybrid vigour of each type is
#     not read as dominance at every marker.
#   - Scans (scan1, REML, LOCO kinship): additive, on allele probabilities (1 d.f.); and full, on the three genotypes
#     (2 d.f.). Dominance is the full LOD minus the additive LOD (1 d.f.); only the F1s inform it.
#   - Signed -log10 P: additive signed by D - B (positive when the D allele raises lifespan); dominance signed by the
#     dominance deviation, BD - (BB + DD) / 2 (positive when heterozygotes live longer than the mean of the homozygotes).
# Writes output/qtl2_<CD|HF>_meansF1_<add|dom>_Progressive.txt (markers x T-ages, Danny's layout),
# output/qtl2_meansF1_N.txt (strains and F1 types per diet and T-age) and output/strainF1_chr1.txt (proximal Chr 1:
# where each F1 type is heterozygous, and the peaks with and without F1s)
#   Rscript qtl2_strainF1.R [path to the F1 table]
#

source("shared.R")
source("qtl2_cross.R")                                                              # strains, cross, pr, K and signedP
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
f1file <- if(length(args) >= 1) args[1] else "../bxd_mortality_rates_2026/data/corrected/AgingBXD_F1_BodyWeight_NaturalDeath_corrected.txt"
tps <- seq(20, 770, 30)
short <- c("BXD9", "BXD170", "BXD172", "BXD197"); long <- c("BXD65", "BXD73b", "BXD77", "BXD102")

f1 <- read.delim(f1file, check.names = FALSE)
f1 <- f1[f1[, "CauseOfDeath"] == "Natural" & f1[, "DietCode"] %in% c("CD", "HF"), ]
types <- sort(unique(f1[, "StrainName"]))
parents <- t(sapply(types, function(x) c(sub("x.*", "", x), paste0("BXD", sub(".*x", "", x)))))
stopifnot(all(parents %in% strains))
f1class <- ifelse(parents[, 1] %in% short & parents[, 2] %in% short, "short", ifelse(parents[, 1] %in% long & parents[, 2] %in% long, "long", NA))
names(f1class) <- types
stopifnot(!any(is.na(f1class)))

### Genotype and allele probabilities for strains and F1 types, and LOCO kinship over both
units <- c(strains, types)
G3 <- pr; AP <- pr
for(chr in names(pr)){
  p <- pr[[chr]]; pD <- p[, "DD", ]; pB <- p[, "BB", ]
  g <- array(0, c(length(units), 3, dim(p)[3]), list(units, c("BB", "BD", "DD"), dimnames(p)[[3]]))
  g[strains, "BB", ] <- pB[strains, ]; g[strains, "DD", ] <- pD[strains, ]
  p1 <- parents[, 1]; p2 <- parents[, 2]
  g[types, "BB", ] <- pB[p1, ] * pB[p2, ]
  g[types, "BD", ] <- pB[p1, ] * pD[p2, ] + pD[p1, ] * pB[p2, ]
  g[types, "DD", ] <- pD[p1, ] * pD[p2, ]
  G3[[chr]] <- g
  a <- array(0, c(length(units), 2, dim(p)[3]), list(units, c("B", "D"), dimnames(p)[[3]]))
  a[, "B", ] <- g[, "BB", ] + g[, "BD", ] / 2; a[, "D", ] <- g[, "DD", ] + g[, "BD", ] / 2
  AP[[chr]] <- a
}
attr(G3, "crosstype") <- "f2"; attr(G3, "alleleprobs") <- FALSE
attr(AP, "alleleprobs") <- TRUE
S <- lapply(AP, function(a){ m <- matrix(a, nrow = dim(a)[1]); tcrossprod(m) })      # allele sharing summed over a chromosome's markers
nm <- sapply(AP, function(a) dim(a)[3])
Kloco <- lapply(names(AP), function(chr){ k <- (Reduce(`+`, S) - S[[chr]]) / (sum(nm) - nm[chr]); dimnames(k) <- list(units, units); k })
names(Kloco) <- names(AP)

### One diet and T-age
one <- function(diet, T){
  x <- ld[which(ld[, "Diet"] == diet & ld[, "AgeAtDeath..days."] >= T),]
  n <- table(x[, "StrainName"]); s <- names(n)[n >= 3]
  y <- tapply(x[, "AgeAtDeath..days."], x[, "StrainName"], mean)[s]
  z <- f1[f1[, "DietCode"] == diet & f1[, "AgeAtDeath (days)"] >= T, ]
  n1 <- table(z[, "StrainName"]); t1 <- names(n1)[n1 >= 3]
  y1 <- tapply(z[, "AgeAtDeath (days)"], z[, "StrainName"], mean)[t1]
  u <- c(s, t1)
  Y <- matrix(c(y, y1), dimnames = list(u, "lifespan"))
  cls <- c(rep("inbred", length(s)), f1class[t1])
  X <- cbind(short = as.numeric(cls == "short"), long = as.numeric(cls == "long")); rownames(X) <- u
  X <- X[, colSums(X) > 0, drop = FALSE]
  if(ncol(X) == 0) X <- NULL
  Kd <- decomp_kinship(lapply(Kloco, function(k) k[u, u]))
  ap <- AP[u, ]; g3 <- G3[u, ]
  add <- scan1(ap, Y, kinship = Kd, addcovar = X)
  full <- scan1(g3, Y, kinship = Kd, addcovar = X)
  eff <- lapply(names(ap), function(chr){
    ca <- scan1coef(ap[, chr], Y, kinship = Kd[[chr]], addcovar = X)
    cg <- scan1coef(g3[, chr], Y, kinship = Kd[[chr]], addcovar = X)
    bd <- if("BD" %in% colnames(cg)) cg[, "BD"] else rep(0, nrow(cg))
    cbind(add = ca[, "D"] - ca[, "B"], dom = bd - (cg[, "BB"] + cg[, "DD"]) / 2)
  })
  eff <- do.call(rbind, eff)
  list(add = signedP(add[, 1], eff[, "add"]), dom = signedP(pmax(full[, 1] - add[, 1], 0), eff[, "dom"]), lod.add = add[, 1],
       n = c(strains = length(s), F1 = length(t1)))
}
jobs <- expand.grid(T = tps, diet = c("CD", "HF"), stringsAsFactors = FALSE)
all <- mclapply(seq_len(nrow(jobs)), function(i) one(jobs[i, "diet"], jobs[i, "T"]), mc.cores = 4)
for(d in c("CD", "HF")) for(m in c("add", "dom")){
  P <- sapply(which(jobs[, "diet"] == d), function(i) all[[i]][[m]])
  out <- cbind(map, P); colnames(out) <- c("Chr", "Locus", "cM", "Mb", tps)
  write.table(out, paste0("output/qtl2_", d, "_meansF1_", m, "_Progressive.txt"), sep = "\t", quote = FALSE)
}
N <- do.call(rbind, lapply(seq_len(nrow(jobs)), function(i) c(diet = jobs[i, "diet"], T = jobs[i, "T"], all[[i]]$n)))
write.table(N, "output/qtl2_meansF1_N.txt", sep = "\t", quote = FALSE, row.names = FALSE)

### Proximal Chr 1: where each F1 type is heterozygous, and the peaks with and without F1s
sink("output/strainF1_chr1.txt", split = TRUE)
cat("Proximal Chr 1 (0 to 60 Mb). Segments where each F1 type's parents differ (B against D calls of BXD_current_rev050423.geno):\n")
k <- which(map[, "Chr"] == "1" & map[, "Mb"] <= 60)
for(t in types){
  a <- geno[k, parents[t, 1]]; b <- geno[k, parents[t, 2]]
  het <- !is.na(a) & !is.na(b) & a %in% c("B", "D") & b %in% c("B", "D") & a != b
  r <- rle(het); e <- cumsum(r$lengths); st <- e - r$lengths + 1
  seg <- paste(sprintf("%.1f-%.1f", map[k[st[r$values]], "Mb"], map[k[e[r$values]], "Mb"]), collapse = ", ")
  cat(sprintf("  %-12s (%-5s x %-5s, %s x %s): %s\n", t, parents[t, 1], parents[t, 2], f1class[t], f1class[t], if(any(het)) seg else "none"))
}
lod <- function(p) qchisq(pmin(10^-abs(p), 1), 1, lower.tail = FALSE) / (2 * log(10))
cat("\nStrongest marker on proximal Chr 1, with its 1.5-LOD support interval at that T-age: strains only (equal weights,\n")
cat("qtl2_strainmeans.R) against strains and F1s (additive). Positive: the D allele lives longer.\n")
for(d in c("CD", "HF")){
  so <- read.table(paste0("output/qtl2_", d, "_means_equal_Progressive.txt"), sep = "\t")
  sf <- read.table(paste0("output/qtl2_", d, "_meansF1_add_Progressive.txt"), sep = "\t")
  for(nmap in c("strains only", "strains + F1s")){
    M <- as.matrix((if(nmap == "strains only") so else sf)[rownames(map)[k], paste0("X", tps)])
    w <- which(abs(M) == max(abs(M), na.rm = TRUE), arr.ind = TRUE)[1, ]
    L <- lod(M[, w[2]]); inside <- which(L >= max(L) - 1.5)
    cat(sprintf("  %s, %-14s peak %s at %.1f Mb, T%d: -log10 P %+.2f (LOD %.2f); 1.5-LOD interval %.1f-%.1f Mb (%.1f Mb)\n", d, nmap, rownames(M)[w[1]],
                map[k[w[1]], "Mb"], tps[w[2]], M[w[1], w[2]], L[w[1]], map[k[min(inside)], "Mb"], map[k[max(inside)], "Mb"], map[k[max(inside)], "Mb"] - map[k[min(inside)], "Mb"]))
  }
}
sink()
