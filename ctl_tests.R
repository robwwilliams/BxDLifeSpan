#
# ctl_tests.R - checks of ctl_functions.R, each against Danny Arends's own code or against a committed result
#
# Rob Williams's fork, 26 November 2026. Written by Claude (Opus 5.5). Run from the fork's root, after ctl_inputs.R, with no arguments:
# Rscript ctl_tests.R. A failed check stops the script. Takes about a minute on one core. Writes output/ctl_scan/tests.txt
#

source("ctl_functions.R")
sink("output/ctl_scan/tests.txt", split = TRUE)
check <- function(what, ok){ cat(sprintf("%-4s %s\n", if(isTRUE(ok)) "ok" else "FAIL", what)); if(!isTRUE(ok)) stop("check failed: ", what) }

# ---- Danny's code, copied from CTL.R and Cauchy.R ------------------------------------------------------------------------------------
dannyTest <- function(cor, n){                                              # CTL.R, the loop body; cor = c(rB, rD), n = c(nB, nD)
  z <- .5*log((1.0 + cor)/(1.0 - cor)); df <- n-(length(cor)-1)
  sumOfSq <- sum(df * z^2); sqOfSum <- sum(df * z); denom <- sum(df); Cv <- sumOfSq - (sqOfSum^2/ denom)
  sign(cor[2] - cor[1]) * -log10(pchisq(Cv, 1, 0, FALSE))
}
cauchycombination <- function(p){ qc <- qcauchy(p, lower.tail=FALSE); pcauchy(mean(qc), lower.tail=FALSE) }   # Cauchy.R

# ---- 1. Fisher's z and the chi-square, including the edge cases -------------------------------------------------------------------------
set.seed(1); r <- cbind(runif(5000, -0.99, 0.99), runif(5000, -0.99, 0.99)); n <- cbind(sample(20:400, 5000, TRUE), sample(20:400, 5000, TRUE))
ours <- ctlTest(r[, 1], r[, 2], n[, 1], n[, 2]); his <- sapply(1:5000, function(i) dannyTest(r[i, ], n[i, ]))
fin <- is.finite(his)
check(sprintf("ctlTest equals CTL.R on 5,000 random pairs of correlations (largest difference %.1e)", max(abs(ours[fin] - his[fin]))), max(abs(ours[fin] - his[fin])) < 1e-12)
check(sprintf("where P is below about 1e-308 (%d of 5,000, |r| near 0.99 and n in the hundreds), both give infinity, with the same sign", sum(!fin)), all(ours[!fin] == his[!fin]))
check("0.5 log((1 + r) / (1 - r)) is atanh(r)", max(abs(0.5 * log((1 + r) / (1 - r)) - atanh(r))) < 1e-12)
check("r = 0 in both groups: z = 0 and -log10 P = 0", ctlTest(0, 0, 50, 60) == 0)
check("r = 0 in one group: the same as CTL.R", abs(ctlTest(0, 0.3, 50, 60) - dannyTest(c(0, 0.3), c(50, 60))) < 1e-12)
check("equal r in both groups: -log10 P = 0", abs(ctlTest(0.4, 0.4, 50, 60)) < 1e-12)
check("r = 1 or -1: NA (CTL.R gives NaN)", all(is.na(ctlTest(c(1, -1, 0.2), c(0.2, 0.2, 1), 50, 60))) && is.nan(dannyTest(c(1, 0.2), c(50, 60))))
check("r = NA: NA", is.na(ctlTest(NA, 0.2, 50, 60)))
check("the sign is positive where r(D) > r(B)", ctlTest(0.1, 0.5, 50, 60) > 0 && ctlTest(0.5, 0.1, 50, 60) < 0)

# ---- 2. Cauchy combination ----------------------------------------------------------------------------------------------------------
P <- matrix(10^-runif(25 * 200, 0, 20), 200); P[sample(length(P), 300)] <- NA
his <- apply(P, 1, function(p) -log10(cauchycombination(p[!is.na(p)])))
check(sprintf("cauchyCombine equals Cauchy.R on 200 markers, P from 1 to 1e-20, NA left out (largest difference %.1e)", max(abs(cauchyCombine(P) - his))), max(abs(cauchyCombine(P) - his)) < 1e-9)

# ---- 3. The correlations: C++ against Danny's line, on made-up data with ties --------------------------------------------------------
w <- round(rnorm(300, 30, 4)); y <- round(rnorm(300, 700, 150), -1); s <- sample(1:40, 300, TRUE) - 1L; G <- matrix(sample(0:2, 40 * 50, TRUE, c(.1, .45, .45)), 40)
ow <- order(w); oy <- order(y)
cc <- ctlCor(G, s, 1 / tabulate(s + 1L, 40)[s + 1L], ow - 1L, match(w[ow], unique(w[ow])), oy - 1L, match(y[oy], unique(y[oy])), tabulate(s + 1L, 40))
rr <- t(sapply(1:50, function(m) ctlCorR(w, y, s, G[s + 1L, m])))
check(sprintf("ctlCor (C++) equals cov.wt on ranks (CTL.R) at 50 made-up markers with ties (largest difference %.1e)", max(abs(cc[, 1:2] - rr), na.rm = TRUE)),
      max(abs(cc[, 1:2] - rr), na.rm = TRUE) < 1e-12 && all(!is.na(cc[, 1])))

# ---- 4. The real data: Danny's design against Danny's matrices ------------------------------------------------------------------------
inp <- readInputs("output/ctl_scan")
Tages <- as.numeric(sub("^w", "", grep("^w[0-9]+$", colnames(inp$mice), value = TRUE)))
scans <- buildScans(inp$mice, inp$strains, c(actuarial = NA), Tages, 2)
check(paste(length(scans), "scans in Danny's design: 2 diets x", length(Tages), "T-ages"), length(scans) == 2 * length(Tages))
d <- scans[["actuarial CD T200"]]; L0 <- scanAll(inp$G, list(d))[, 1]; set.seed(2); k <- sample(which(!is.na(L0)), 100)
rr <- t(sapply(k, function(p) ctlCorR(d$w, d$y, d$s, inp$G[d$s + 1L, p])))
cc <- ctlCor(inp$G[, k], d$s, d$wt, d$ow, d$tw, d$oy, d$ty, d$nall)
check(sprintf("ctlCor equals cov.wt on ranks at 100 real markers, chow T200 (largest difference %.1e)", max(abs(cc[, 1:2] - rr))), max(abs(cc[, 1:2] - rr)) < 1e-12)
real <- runMap(0, inp, scans, 20261124)
for(diet in c("CD", "HF")){
  dm <- as.matrix(read.table(paste0("output/CTL_", diet, "_LODmatrix.txt"), sep = "\t", check.names = FALSE))
  m <- real$maps$mice[inp$pat, paste("actuarial", diet, paste0("T", colnames(dm)))]; rownames(m) <- inp$map$Locus; m <- m[rownames(dm), ]
  check(sprintf("%s: the real map equals CTL.R's %s matrix, rounded as he saved it (%d values; largest difference %.3f; cells missing in one only: %d)",
                diet, diet, sum(!is.na(dm)), max(abs(dm - m), na.rm = TRUE), sum(is.na(dm) != is.na(m))), max(abs(dm - round(m, 2)), na.rm = TRUE) <= 0.0100001 && all(is.na(dm) == is.na(m)))
}

# ---- 5. Permutations ------------------------------------------------------------------------------------------------------------------
o <- permuteStrains(inp$stratum, 20261125)
check("a permutation moves strains only within their stratum", all(inp$stratum[o] == inp$stratum) && !all(o == seq_along(o)))
check("the same seed gives the same permutation", all(o == permuteStrains(inp$stratum, 20261125)))

# ---- 6. Against the committed result of ctl_permutations.R (BXD mice, Danny's design, maps 0 to 5) -------------------------------
# ctl_permutations.R took the Cauchy quantile as tan((0.5 - p) pi), which loses digits for P between about 1e-15 and 1e-12, and 1 / (p pi)
# below 1e-15. qcauchy(), in Danny's Cauchy.R and here, does not. So the Cauchy column may differ by a little where a permuted map has P
# values that small; the other columns must agree to the 4 decimals saved.
p <- c(1e-15, 1e-14, 1e-13, 1e-12); oldq <- tan((0.5 - p) * pi); exact <- 1 / tan(p * pi)
cat(sprintf("     the old Cauchy quantile at P = %s is off by %s%%\n", paste(p, collapse = ", "), paste(signif(100 * (oldq / exact - 1), 2), collapse = ", ")))
check("qcauchy() gives the exact quantile at P = 1e-14", abs(qcauchy(1e-14, lower.tail = FALSE) / exact[2] - 1) < 1e-12)
old <- read.delim("output/ctl_permutations_maxima.txt", check.names = FALSE)
cols <- c(max = "any", cauchy = "cauchy", cauchy.Bonf = "cauchy.Bonf", cauchy.BH = "cauchy.BH", Bonf = "Bonf.all.markers", BH = "BH.all.markers", lambda = "lambda")
for(i in 0:5){
  s <- if(i == 0) real$summary else runMap(i, inp, scans, 20261124)$summary
  dif <- sapply(c("CD", "HF"), function(diet) abs(s[paste("actuarial", diet, names(cols), sep = ".")] - unlist(old[old$perm == i, paste(diet, cols, sep = ".")]))); rownames(dif) <- names(cols)
  ca <- rownames(dif) == "cauchy"
  check(sprintf("map %d (%s) gives the committed max, Bonferroni, B&H and lambda (largest difference %.4f), and the Cauchy maximum within 0.02 (%.4f)",
                i, if(i == 0) "real" else "permuted", max(dif[!ca, ]), max(dif[ca, ])), max(dif[!ca, ]) < 1e-3 && max(dif[ca, ]) < 0.02)
}

# ---- 7. The strain-level map: C++ on strain units against plain R on strain means ----------------------------------------------------
d <- scans[["actuarial CD T140"]]; x <- inp$mice[match(d$id, inp$mice$mouse), ]
check("the scan's mice are those of the mice table", all(x$w140 == d$w) && all(x$lifespan == d$y))
mice.per.strain <- table(x$strain); kept <- names(mice.per.strain)[mice.per.strain >= 2]                   # strains with 2 or more mice
strain.weight <- tapply(x$w140, x$strain, mean)[kept]; strain.life <- tapply(x$lifespan, x$strain, mean)[kept]
LS <- scanAllStrains(inp$G, list(d))$coded[, 1]; set.seed(3); k <- sample(which(!is.na(LS)), 100)
plainR <- sapply(k, function(p){ g <- inp$G[match(kept, inp$strains), p]
  r <- sapply(1:2, function(j) cor(strain.weight[g == j], strain.life[g == j], method = "spearman")); ctlTest(r[1], r[2], sum(g == 1), sum(g == 2)) })
check(sprintf("strain-level map equals Spearman of strain means in plain R at 100 markers, chow T140 (largest difference %.1e)", max(abs(LS[k] - plainR))), max(abs(LS[k] - plainR)) < 1e-10)
peak <- inp$pat[inp$map$Locus == "rs13481185"]
check(sprintf("chow Chr 11 (rs13481185), T140: -log10 P %.2f with mice as n, %.2f with strains as n", abs(scanAll(inp$G, list(d))[peak, 1]), abs(LS[peak])), abs(LS[peak]) < abs(scanAll(inp$G, list(d))[peak, 1]))

# ---- 8. The standard test: formula, edge cases and calibration with no CTL -----------------------------------------------------------
standardByHand <- function(rB, rD, nB, nD){ z <- 0.5 * log((1 + c(rB, rD)) / (1 - c(rB, rD))); w <- (c(nB, nD) - 3) / 1.06
  Cv <- sum(w * z^2) - sum(w * z)^2 / sum(w); sign(rD - rB) * -log10(pchisq(Cv, 1, lower.tail = FALSE)) }
r <- cbind(runif(2000, -0.95, 0.95), runif(2000, -0.95, 0.95)); n <- cbind(sample(20:120, 2000, TRUE), sample(20:120, 2000, TRUE))
d1 <- max(abs(ctlTestStandard(r[, 1], r[, 2], n[, 1], n[, 2]) - sapply(1:2000, function(i) standardByHand(r[i, 1], r[i, 2], n[i, 1], n[i, 2]))))
check(sprintf("ctlTestStandard equals the formula written out, weights (n - 3) / 1.06, on 2,000 random cases (largest difference %.1e)", d1), d1 < 1e-12)
check("ctlTestStandard: r = 0 in both gives 0; |r| = 1 gives NA; fewer than 20 strains in a group gives NA",
      ctlTestStandard(0, 0, 30, 30) == 0 && is.na(ctlTestStandard(1, 0.2, 30, 30)) && is.na(ctlTestStandard(0.2, 0.3, 19, 40)))
set.seed(4); sims <- t(replicate(20000, c(cor(rnorm(20), rnorm(20), method = "spearman"), cor(rnorm(20), rnorm(20), method = "spearman"))))
p <- 10^-abs(ctlTestStandard(sims[, 1], sims[, 2], 20, 20)); infl <- median(qchisq(p, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
check(sprintf("with no CTL and 20 strains per group, the standard test is calibrated: inflation %.2f, P < 0.01 in %.2f%% (expected 1%%)", infl, 100 * mean(p < 0.01)),
      infl > 0.85 && infl < 1.1 && mean(p < 0.01) < 0.015)
sink()
