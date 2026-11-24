#
# chr1_twolocus.R - One locus or two on proximal Chr 1 (0 to 50 Mb) on high fat: a two-locus scan of strain means alone and
# of strain and F1 means together
#
# Rob Williams's fork, 23 November 2026. Exploratory. Asked for by Rob ("Yes, run both", 23 November 2026, to: "a two-locus
# scan of proximal Chr 1 on high fat (about 23.5 and 37.6 Mb)"). chr1_zoom.R found peaks near 23.5 Mb (early T-ages) and
# 37.6 Mb (later T-ages, with the F1s), each within 1.5 LOD of the other at most T-ages.
#   - Units, trait, covariate and kinship as in qtl2_strainF1.R (strainF1_setup.R): unit means with at least 3 mice alive
#     at T, equal weights, F1 class as a covariate (strains alone: none), kinship leaving out Chr 1. T20 to T620, the
#     stable range of the high-fat strain-means thresholds.
#   - Model: the linear mixed model of scan1, written out so that two loci can be fitted at once. The heritability is
#     estimated by REML under the model without loci, as scan1 does, and held there; each model is then generalized
#     least squares, and its LOD is n/2 log10(RSS0 / RSS). Additive allele dosage only, 1 d.f. per locus. Markers with
#     the same dosages in every unit are fitted once.
#   - For each T-age: the best single locus; the best pair of loci; LOD(2 v 1), the best pair against the best single
#     locus; and the best second locus given a locus at 23.5 Mb, and given one at 37.6 Mb.
#   - P(2 v 1): parametric bootstrap under the best single-locus model, its effect, heritability and variance as fitted;
#     each simulated trait goes through the whole search (best single locus, best pair). 1,000 per T-age and map, job i
#     drawing from seed 20261124 + i. This
#     asks whether a second locus anywhere in 0-50 Mb adds more than one locus would by chance.
# Writes output/chr1_twolocus.txt and output/maps/Chr1_twolocus.svg
#   Rscript chr1_twolocus.R [bootstraps, default 1000] [cores, default 8]
#

library(svglite)
source("shared.R")
source("qtl2_cross.R")                                                              # strains, pr
source("strainF1_setup.R")                                                         # units, AP, Kloco, unitMeans, classCovar
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
nboot <- if(length(args) >= 1) as.numeric(args[1]) else 1000
cores <- if(length(args) >= 2) as.numeric(args[2]) else 8
tps <- seq(20, 620, 30)
fixMb <- c(23.52, 37.63)                                                           # the two peaks of chr1_zoom.R
col.strains <- "grey35"; col.F1 <- "#008300"
lp <- function(lod) -pchisq(2 * log(10) * pmax(lod, 0), 1, lower.tail = FALSE, log.p = TRUE) / log(10)

k <- dimnames(AP[["1"]])[[3]]; k <- k[map[k, "Mb"] <= 50]; k <- k[order(map[k, "Mb"])]
mb <- as.numeric(map[k, "Mb"])
Dall <- AP[["1"]][, "D", k]                                                        # units x markers: D allele dosage, 0 to 1

### REML heritability with the kinship eigendecomposed: V = h d + (1 - h)
remlH <- function(yr, Xr, d){
  f <- function(h){ v <- h * d + 1 - h; w <- 1 / sqrt(v); Xw <- Xr * w; yw <- yr * w
    q <- qr(Xw); r <- qr.resid(q, yw); n <- length(yr); p <- q$rank
    -0.5 * (sum(log(v)) + 2 * sum(log(abs(diag(qr.R(q))[1:p]))) + (n - p) * log(sum(r^2))) }
  optimize(f, c(0, 0.999), maximum = TRUE)$maximum
}

### One analysis of a trait y (units u): LODs of every single locus and pair, given the rotation and weights
analyse <- function(yw, X0w, Gw){
  n <- length(yw); q0 <- qr(X0w); r0 <- qr.resid(q0, yw); rss0 <- sum(r0^2)
  Gr <- qr.resid(q0, Gw)                                                           # markers with the covariates projected out
  gg <- colSums(Gr^2); gy <- colSums(Gr * r0)
  rss1 <- rss0 - gy^2 / gg; L1 <- n / 2 * log10(rss0 / rss1)
  L2 <- matrix(NA, ncol(Gw), ncol(Gw))
  for(i in seq_len(ncol(Gw))){                                                     # add locus j to locus i
    gi <- Gr[, i]; ri <- r0 - gi * gy[i] / gg[i]
    Grj <- Gr - outer(gi, colSums(Gr * gi) / gg[i])
    L2[i, ] <- n / 2 * log10(rss0 / (rss1[i] - colSums(Grj * ri)^2 / pmax(colSums(Grj^2), 1e-12)))
    L2[i, colSums(Grj^2) < 1e-8] <- L1[i]                                          # the same locus twice adds nothing
  }
  list(L1 = L1, L2 = L2)
}

one <- function(T, with.f1, seed){
  set.seed(seed)
  Y <- unitMeans("HF", T, strains.only = !with.f1); u <- rownames(Y)[!is.na(Y[, 1])]
  y <- Y[u, 1]; X0 <- cbind(rep(1, length(u)), if(with.f1) classCovar(u) else NULL)
  e <- eigen(Kloco[["1"]][u, u], symmetric = TRUE); d <- pmax(e$values, 0); U <- e$vectors
  G <- Dall[u, , drop = FALSE]
  key <- apply(round(G, 6), 2, paste, collapse = ",")
  pat <- !duplicated(key)                                                           # one column per distinct pattern
  pid <- match(key, key[pat])                                                       # each marker's pattern
  yr <- drop(crossprod(U, y)); Xr <- crossprod(U, X0); Gr <- crossprod(U, G[, pat, drop = FALSE])
  h <- remlH(yr, Xr, d); w <- 1 / sqrt(h * d + 1 - h)
  a <- analyse(yr * w, Xr * w, Gr * w)
  i1 <- which.max(a$L1); best <- which(a$L2 == max(a$L2, na.rm = TRUE), arr.ind = TRUE)[1, ]
  fx <- sapply(fixMb, function(m) pid[which.min(abs(mb - m))])
  cond <- lapply(fx, function(f) a$L2[f, ] - a$L1[f])                              # LOD of a second locus, given f
  ### Bootstrap under the best single locus
  Xa <- cbind(Xr, Gr[, i1]) * w; beta <- qr.coef(qr(Xa), yr * w)
  h1 <- remlH(yr, cbind(Xr, Gr[, i1]), d); v1 <- h1 * d + 1 - h1
  s2 <- sum(qr.resid(qr(cbind(Xr, Gr[, i1]) / sqrt(v1)), yr / sqrt(v1))^2) / (length(y) - ncol(Xa))
  mu <- drop(cbind(Xr, Gr[, i1]) %*% beta)
  boot <- sapply(seq_len(nboot), function(b){
    ys <- mu + sqrt(s2 * v1) * rnorm(length(y))                                    # simulated in the rotated space
    hb <- remlH(ys, Xr, d); wb <- 1 / sqrt(hb * d + 1 - hb)
    ab <- analyse(ys * wb, Xr * wb, Gr * wb); max(ab$L2, na.rm = TRUE) - max(ab$L1)
  })
  pm <- function(j) mb[which(pid == j)]                                             # the span of markers of a pattern
  list(T = T, f1 = with.f1, n = length(u), nF1 = sum(u %in% types), h = h, patterns = sum(pat),
       single = c(range(pm(i1)), a$L1[i1]), pair = c(range(pm(best[1])), range(pm(best[2])), a$L2[best[1], best[2]]),
       twov1 = a$L2[best[1], best[2]] - a$L1[i1], p = (1 + sum(boot >= a$L2[best[1], best[2]] - a$L1[i1])) / (nboot + 1),
       cond = lapply(cond, function(cl){ j <- which.max(cl); c(range(pm(j)), cl[j]) }),
       L1 = a$L1[pid], condL = lapply(cond, function(cl) cl[pid]), fixL = a$L1[fx])
}

jobs <- expand.grid(T = tps, f1 = c(FALSE, TRUE))
res <- mclapply(seq_len(nrow(jobs)), function(i) one(jobs$T[i], jobs$f1[i], seed = 20261124 + i), mc.cores = cores, mc.preschedule = FALSE)   # job i draws from seed 20261124 + i

### Check: single-locus LODs against scan1 (chr1_zoom.R's maps), strains alone at T20 and strains + F1s at T500
chk <- function(file, T, f1){
  x <- read.table(file.path("output", file), sep = "\t")[k, paste0("X", T)]
  r <- res[[which(jobs$T == T & jobs$f1 == f1)]]; max(abs(abs(x) - lp(r$L1)))
}
check <- c(chk("qtl2_HF_means_equal_Progressive.txt", 20, FALSE), chk("qtl2_HF_meansF1_add_Progressive.txt", 500, TRUE))

sink("output/chr1_twolocus.txt", split = TRUE)
cat("Proximal Chr 1 (0 to 50 Mb) on high fat: one locus or two. Unit means with at least 3 mice alive at T; F1 class as a\n")
cat("covariate when the F1s are in; kinship leaving out Chr 1; heritability (h) by REML without loci. LODs on 1 d.f. per locus.\n")
cat(sprintf("Check against scan1: largest difference in -log10 P of single loci %.4f (strains, T20) and %.4f (strains + F1s, T500)\n\n", check[1], check[2]))
cat("single: best single locus (Mb span of its genotype pattern) and LOD. pair: best two loci and their joint LOD.\n")
cat("2v1: LOD of the best pair over the best single locus; P: parametric bootstrap under the best single locus (", nboot, ").\n")
cat("given 23.5 / given 37.6: the best second locus given a locus at that position, and its LOD.\n\n")
tab <- do.call(rbind, lapply(res, function(r) data.frame(
  units = if(r$f1) "strains+F1" else "strains", T = r$T, n = r$n, F1 = r$nF1, h = round(r$h, 2),
  single = sprintf("%.1f-%.1f", r$single[1], r$single[2]), LOD1 = round(r$single[3], 2),
  pair = sprintf("%.1f-%.1f + %.1f-%.1f", r$pair[1], r$pair[2], r$pair[3], r$pair[4]), LOD2 = round(r$pair[5], 2),
  `2v1` = round(r$twov1, 2), P = signif(r$p, 2),
  given23.5 = sprintf("%.1f-%.1f %.2f", r$cond[[1]][1], r$cond[[1]][2], r$cond[[1]][3]),
  given37.6 = sprintf("%.1f-%.1f %.2f", r$cond[[2]][1], r$cond[[2]][2], r$cond[[2]][3]), check.names = FALSE)))
print(tab, row.names = FALSE)
sink()

### Figure: LOD profiles at three T-ages, alone and given each locus; and P(2 v 1) by T-age
views <- c(20, 260, 500)
svglite("output/maps/Chr1_twolocus.svg", width = 7, height = 2.35 * length(views) * 2 + 2.4, bg = "white", system_fonts = list(sans = "Arial"), fix_text_size = FALSE)
layout(matrix(1:(2 * length(views) + 1)), heights = c(rep(2.35, 2 * length(views)), 2.4))
par(mai = c(0.4, 0.75, 0.25, 0.2), cex = 0.75, mgp = c(2.2, 0.6, 0), tcl = -0.3)
lt <- 0
for(T in views) for(f1 in c(FALSE, TRUE)){
  r <- res[[which(jobs$T == T & jobs$f1 == f1)]]; co <- if(f1) col.F1 else col.strains; lt <- lt + 1
  plot(NA, xlim = c(0, 50), ylim = c(0, 5.5), xaxs = "i", xaxt = "n", bty = "l", xlab = "", ylab = "LOD", las = 1)
  axis(1, at = seq(0, 50, 5), mgp = c(2, 0.4, 0)); mtext("Chr 1 (Mb)", side = 1, line = 1.3, at = 50, adj = 1, cex = 0.7)
  mtext(sprintf("High fat, T%d, %s (%d units)", T, if(f1) "strains + F1 types" else "strains alone", r$n), side = 3, line = 0.1, adj = 0, font = 2, cex = 0.8)
  mtext(letters[lt], side = 3, line = 0.1, at = grconvertX(0.05, "inches", "user"), adj = 0, font = 2, cex = 1, xpd = NA)
  lines(mb, r$L1, col = co, lwd = 1.5)
  lines(mb, pmax(r$condL[[1]], 0), col = co, lwd = 1, lty = 2)
  lines(mb, pmax(r$condL[[2]], 0), col = co, lwd = 1, lty = 3)
  abline(v = fixMb, col = "grey75", lwd = 0.5)
  legend("topleft", bty = "n", cex = 0.75, lwd = c(1.5, 1, 1), lty = 1:3, col = co,
         legend = c("one locus", "second locus, given 23.5 Mb", "second locus, given 37.6 Mb"))
  text(50, 5.3, sprintf("LOD(2 v 1) %.2f, P = %.3f", r$twov1, r$p), adj = c(1, 1), cex = 0.8)
}
par(mai = c(0.55, 0.75, 0.25, 0.2))
P <- sapply(res, function(r) r$p)
plot(NA, xlim = range(tps), ylim = c(0, 1), xlab = "T-age (days)", ylab = "P (2 v 1)", las = 1, bty = "l", xaxt = "n")
axis(1, at = seq(20, 620, 60))
mtext("Two loci against one, by T-age (parametric bootstrap)", side = 3, line = 0.1, adj = 0, font = 2, cex = 0.8)
mtext(letters[lt + 1], side = 3, line = 0.1, at = grconvertX(0.05, "inches", "user"), adj = 0, font = 2, cex = 1, xpd = NA)
abline(h = 0.05, lty = 2, col = "grey55")
for(f1 in c(FALSE, TRUE)){ i <- jobs$f1 == f1; lines(jobs$T[i], P[i], type = "b", pch = 19, cex = 0.6, col = if(f1) col.F1 else col.strains) }
legend("topright", bty = "n", cex = 0.75, lwd = 1, pch = 19, col = c(col.strains, col.F1), legend = c("strains alone", "strains + F1 types"))
invisible(dev.off())
cat("Wrote output/chr1_twolocus.txt and output/maps/Chr1_twolocus.svg\n")
