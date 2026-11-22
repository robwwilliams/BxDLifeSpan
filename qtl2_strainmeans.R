#
# qtl2_strainmeans.R - Actuarial maps of strain means with R/qtl2 and LOCO kinship, under three weightings
#
# Rob Williams's fork, 21 November 2026. Exploratory. The plan is Rob's: "Let's look at the three methods of weight:
# 1. No weighting but a minimum of 3 per strain per diet; 2. Weighting using Danny's method; Weighting use the method that
# you suggested (variance or even SEM, which ever you think a real statistician would consider acceptable). I frankly
# like the first method since it has been our default". And earlier: "The mean survivorship data MUST be mapped with
# kinship to be accurate."
#   - At each T-age (20 to 770 days, 30-day steps) and diet, the BXD strains with at least 3 mice alive at T (shared.R's
#     mice, as in Danny's maps). The trait is each strain's mean age at death of those mice; no covariates.
#   - Every scan is R/qtl2's linear mixed model with LOCO kinship (scan1, REML). Only the weights differ:
#       equal  every strain the same
#       N      the number of mice, as in Danny's lifespan.R
#       var    the inverse of the variance of each strain mean, 1 / (s2.strain + s2.within / n). s2.within is the pooled
#              variance among mice of a strain; s2.strain is the variance between strains that kinship does not explain,
#              estimated by REML for the strain means with V = s2.kin K + diag(s2.strain + s2.within / n). The weights
#              of 1 / SEM^2, from each strain's own SD, are not used: an SD from 3 mice is too noisy, and they ignore
#              s2.strain.
#   - Signed -log10 P as in qtl2_individuals.R: chi-square on 1 d.f., with the sign of the D - B effect from scan1coef
#     (same weights and kinship), positive when the D allele raises lifespan.
# Writes output/qtl2_<CD|HF>_means_<equal|N|var>_Progressive.txt (markers x T-ages, Danny's layout),
# output/qtl2_means_N.txt (strains per diet and T-age) and output/qtl2_means_variance.txt (the variance components and
# the spread of the weights)
#   Rscript qtl2_strainmeans.R        (about 10 minutes)
#

source("shared.R")
source("qtl2_cross.R")
library(parallel)

tps <- seq(20, 770, 30)
schemes <- c("equal", "N", "var")
Kall <- calc_kinship(pr, "overall")                                                 # for the variance components

# REML for strain means: V = s2.kin * K + diag(s2.strain + d), d = s2.within / n known; intercept only
remlMeans <- function(y, Ks, d){
  X <- matrix(1, length(y), 1)
  nll <- function(par){
    V <- exp(par[1]) * Ks + diag(exp(par[2]) + d, length(y))
    C <- chol(V); Vi <- chol2inv(C)
    XtViX <- crossprod(X, Vi %*% X)
    r <- y - X %*% solve(XtViX, crossprod(X, Vi %*% y))
    0.5 * (2 * sum(log(diag(C))) + log(det(XtViX)) + sum(r * (Vi %*% r)))
  }
  v <- var(y)
  o <- optim(c(log(v / 2), log(v / 2)), nll, control = list(maxit = 2000))
  exp(o$par)
}

one <- function(diet, T){
  x <- ld[which(ld[, "Diet"] == diet & ld[, "AgeAtDeath..days."] >= T),]
  n <- table(x[, "StrainName"]); s <- names(n)[n >= 3]; n <- as.numeric(n[s]); names(n) <- s
  x <- x[x[, "StrainName"] %in% s,]
  y <- tapply(x[, "AgeAtDeath..days."], x[, "StrainName"], mean)[s]
  s2.within <- sum((x[, "AgeAtDeath..days."] - y[x[, "StrainName"]])^2) / (nrow(x) - length(s))
  vc <- remlMeans(as.numeric(y), Kall[s, s], s2.within / n)
  w <- list(equal = NULL, N = n, var = 1 / (vc[2] + s2.within / n))
  Y <- matrix(y, dimnames = list(s, "lifespan"))
  Ks <- lapply(K, function(k) k[s, s])
  prs <- pr[s, ]
  res <- lapply(schemes, function(m){
    wt <- if(is.null(w[[m]])) NULL else setNames(as.numeric(w[[m]]) / mean(w[[m]]), s)
    out <- scan1(prs, Y, kinship = Ks, weights = wt)
    eff <- unlist(lapply(names(prs), function(chr){ co <- scan1coef(prs[, chr], Y, kinship = Ks[[chr]], weights = wt); co[, "DD"] - co[, "BB"] }))
    list(p = signedP(out[, 1], eff), hsq = attr(out, "hsq")[1, 1])
  })
  names(res) <- schemes
  wv <- w$var
  list(res = res, strains = length(s), info = c(diet = diet, T = T, strains = length(s), mice = nrow(x), sd.within = round(sqrt(s2.within), 1),
       sd.strain = round(sqrt(vc[2]), 1), sd.kin = round(sqrt(vc[1]), 1), ratio.strain.within = round(vc[2] / s2.within, 3),
       var.weight.max.min = round(max(wv) / min(wv), 2), N.weight.max.min = round(max(n) / min(n), 2),
       hsq.equal = round(res$equal$hsq, 3), hsq.N = round(res$N$hsq, 3), hsq.var = round(res$var$hsq, 3)))
}

jobs <- expand.grid(T = tps, diet = c("CD", "HF"), stringsAsFactors = FALSE)
all <- mclapply(seq_len(nrow(jobs)), function(i) one(jobs[i, "diet"], jobs[i, "T"]), mc.cores = 4)
info <- do.call(rbind, lapply(all, function(a) a$info))
write.table(info, "output/qtl2_means_variance.txt", sep = "\t", quote = FALSE, row.names = FALSE)
nst <- sapply(c("CD", "HF"), function(d) sapply(which(jobs[, "diet"] == d), function(i) all[[i]]$strains)); rownames(nst) <- tps
write.table(nst, "output/qtl2_means_N.txt", sep = "\t", quote = FALSE)
for(d in c("CD", "HF")) for(m in schemes){
  P <- sapply(which(jobs[, "diet"] == d), function(i) all[[i]]$res[[m]]$p)
  out <- cbind(map, P); colnames(out) <- c("Chr", "Locus", "cM", "Mb", tps)
  write.table(out, paste0("output/qtl2_", d, "_means_", m, "_Progressive.txt"), sep = "\t", quote = FALSE)
}
print(info[, c("diet", "T", "strains", "sd.within", "sd.strain", "sd.kin", "var.weight.max.min", "N.weight.max.min", "hsq.equal", "hsq.N", "hsq.var")], quote = FALSE)
