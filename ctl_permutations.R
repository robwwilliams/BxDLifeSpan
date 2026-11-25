#
# ctl_permutations.R - Danny Arends's CTL maps of the BXD mice reproduced, and tested with B&H, Bonferroni and permutations of strains
#
# Rob Williams's fork, 24 November 2026. Exploratory. Written by Claude (Opus 5.5) for Rob and Danny. Rob: "We can run it both his
# way to make sure we understand his way, and then test drive permutation tests. In Arends et al. we generally just used Bonferroni or
# B&H corrections rather than permutations."
#   - His way (CTL.R, unchanged): at each marker and T-age, weight at T and lifespan, each adjusted for vivarium and set-up class, are
#     ranked within the B mice and within the D mice, and correlated with weights 1 / mice of the strain. Fisher's z of the two
#     correlations are compared by chi-square on 1 d.f., with sample sizes the numbers of mice (n - 1 each). The arithmetic is
#     rewritten in C++ for speed and checked against output/CTL_<diet>_LODmatrix.txt (2 decimals).
#   - His corrections (Arends et al. 2026, Methods): B&H at 5% FDR and Bonferroni, over the markers of one T-age, and over all markers and
#     T-ages of a diet.
#   - Permutations: each strain takes the genotypes of another strain of its stratum (Taylor's strains, BXD1 to BXD42, or the newer
#     ones), the same one at every T-age and on both diets. Every mouse keeps its weight, lifespan and strain. So the mice of a strain
#     move together, and the correlation within and between strains is kept. Recorded per permutation: the largest -log10 P for each
#     diet and T-age and over the T-ages, and how many markers B&H and Bonferroni call.
# Writes output/ctl_permutations_maxima.txt (one row per permutation; 0 = real data) and output/ctl_permutations.txt
#   Rscript ctl_permutations.R [permutations, 1000] [cores, 12] [mode, strains]
# Modes (output names take the mode as a suffix, except strains):
#   strains  as above
#   mice     a null without strain structure: the (weight, lifespan) pairs are shuffled among the mice of each diet and T-age, and the
#            strains keep their genotypes. If P values that count mice were right, this would give no inflation beyond the weights.
#   within   weight and lifespan are first centered on their strain means (strains with at least 2 mice), so that only differences among
#            the mice of a strain are compared; permutations of strains as above. A different question, not Danny's statistic.
#   spike    as strains, but a CTL is planted at one marker (chow and high fat, T200 to T500): within D mice the lifespan is
#            replaced by its value plus 25 days per SD of weight within the strain. Does the permutation test find it?
#

library(Rcpp); library(parallel)
args <- commandArgs(trailingOnly = TRUE)
nperm <- if(length(args) >= 1) as.numeric(args[1]) else 1000
cores <- if(length(args) >= 2) as.numeric(args[2]) else 12
mode <- if(length(args) >= 3) args[3] else "strains"
sfx <- if(mode == "strains") "" else paste0("_", mode)
spike <- if(length(args) >= 4) as.numeric(args[4]) else 25                                  # days of lifespan per SD of weight, in D mice
spike.marker <- "rsm10000004000"
RNGkind("L'Ecuyer-CMRG"); set.seed(20261124)

sourceCpp(code = '
#include <Rcpp.h>
using namespace Rcpp;
// ranks (ties averaged) within a subset, walking mice in increasing order; tie = id of the value in that order
static void subrank(const IntegerVector& ord, const IntegerVector& tie, const std::vector<signed char>& g, signed char k, std::vector<double>& r){
  int n = ord.size(), i = 0, seen = 0;
  while(i < n){
    int t = tie[i], j = i, c = 0;
    while(j < n && tie[j] == t){ if(g[ord[j]] == k) c++; j++; }
    double rk = seen + (c + 1) / 2.0;
    for(int q = i; q < j; q++) if(g[ord[q]] == k) r[ord[q]] = rk;
    seen += c; i = j;
  }
}
// G: strains x patterns (1 B, 2 D, 0 other). Returns, per pattern: rB, rD, nB, nD.
// [[Rcpp::export]]
NumericMatrix ctlScan(const IntegerMatrix& G, const IntegerVector& strain, const NumericVector& wt, const IntegerVector& ow, const IntegerVector& tw,
                      const IntegerVector& oy, const IntegerVector& ty, const IntegerVector& nall){
  int np = G.ncol(), n = strain.size(); NumericMatrix out(np, 4); std::vector<signed char> g(n); std::vector<double> a(n), b(n);
  for(int p = 0; p < np; p++){
    for(int i = 0; i < n; i++) g[i] = (signed char) G(strain[i], p);
    int mB = 0, mD = 0;
    for(int i = 0; i < nall.size(); i++){ if(G(i, p) == 1) mB += nall[i]; else if(G(i, p) == 2) mD += nall[i]; }    // the >= 20 rule counts mice before weights are required
    if(mB < 20 || mD < 20){ for(int k = 0; k < 4; k++) out(p, k) = NA_REAL; continue; }
    for(int k = 1; k <= 2; k++){
      subrank(ow, tw, g, k, a); subrank(oy, ty, g, k, b);
      double s = 0, ma = 0, mb = 0; int c = 0;
      for(int i = 0; i < n; i++) if(g[i] == k){ s += wt[i]; ma += wt[i] * a[i]; mb += wt[i] * b[i]; c++; }
      ma /= s; mb /= s; double va = 0, vb = 0, cv = 0;
      for(int i = 0; i < n; i++) if(g[i] == k){ va += wt[i] * (a[i] - ma) * (a[i] - ma); vb += wt[i] * (b[i] - mb) * (b[i] - mb); cv += wt[i] * (a[i] - ma) * (b[i] - mb); }
      out(p, k - 1) = cv / std::sqrt(va * vb); out(p, k + 1) = c;
    }
  }
  return out;
}')

source("shared.R")
ldata[which(ldata[, "AgeAtSetUp.in.colony..days."] < 0), "AgeAtSetUp.in.colony..days."] <- NA
BWs <- read.table("output/bodyweightInterpolatedData.txt", sep = "\t", quote = "", check.names = FALSE)
invisible(capture.output(suppressWarnings(source("adjust.R"))))
tps <- TPStodo
bxd <- unique(as.character(ldata[grepl("BXD", ldata[, "StrainName"]), "StrainName"])); st <- sort(bxd[bxd %in% colnames(geno)])
G0 <- matrix(match(as.matrix(geno[, st]), c("B", "D"), nomatch = 0L), nrow(geno)); dimnames(G0) <- list(rownames(geno), st)
key <- apply(G0, 1, paste, collapse = ""); u <- !duplicated(key); pat <- match(key, key[u]); G <- t(G0[u, ])                     # strains x distinct patterns
stratum <- ifelse(as.numeric(gsub("[^0-9].*$", "", sub("^BXD", "", st))) <= 42, "Taylor", "newer")
cat(nrow(geno), "markers,", ncol(G), "distinct strain patterns,", length(st), "strains (", sum(stratum == "Taylor"), "Taylor's )\n")

# The mice of each diet and T-age, as in CTL.R; strains without genotypes are left out, and the >= 20 rule counts mice before weights are required
dat <- list()
for(diet in c("CD", "HF")) for(tp in tps){
  iix <- which(ldata[, "AgeAtSetUp.in.colony..days."] <= as.numeric(tp) & ldata[, "AgeAtDeath..days."] >= as.numeric(tp) & grepl("BXD", ldata[, "StrainName"]) & ldata[, "Diet"] == diet)
  x <- ldata[iix, ]; x <- x[rownames(x) %in% rownames(BWsAdj) & as.character(x[, "StrainName"]) %in% st, ]
  s <- match(as.character(x[, "StrainName"]), st); nall <- tabulate(s, length(st))
  w <- BWsAdj[rownames(x), tp]; y <- x[, "Yadj"]; ok <- !is.na(w) & !is.na(y); w <- w[ok]; y <- y[ok]; s <- s[ok]
  if(mode == "within"){ k <- ave(w, s, FUN = length) > 1; w <- w[k]; y <- y[k]; s <- s[k]; w <- w - ave(w, s); y <- y - ave(y, s) }
  if(mode == "spike"){ if(as.numeric(tp) >= 200 && as.numeric(tp) <= 500){ gd <- G0[spike.marker, st[s]] == 2; zw <- (w - ave(w, s)) / sd(w - ave(w, s)); y <- y + ifelse(gd, spike * zw, 0) } }
  n <- tabulate(s, length(st)); ow <- order(w); oy <- order(y)
  dat[[paste(diet, tp)]] <- list(s = s - 1L, wt = 1 / n[s], ow = ow - 1L, tw = as.integer(match(w[ow], unique(w[ow]))), oy = oy - 1L, ty = as.integer(match(y[oy], unique(y[oy]))), nall = nall)
}
scanAll <- function(G, dat. = dat) sapply(dat., function(d){
  r <- ctlScan(G, d$s, d$wt, d$ow, d$tw, d$oy, d$ty, d$nall); z <- atanh(r[, 1:2]); df <- r[, 3:4] - 1
  Cv <- rowSums(df * z^2) - rowSums(df * z)^2 / rowSums(df); sign(r[, 2] - r[, 1]) * -pchisq(Cv, 1, lower.tail = FALSE, log.p = TRUE) / log(10) })

# Real data, against Danny's matrices
real <- scanAll(G)[pat, ]; rownames(real) <- rownames(geno)
if(mode == "strains") for(diet in c("CD", "HF")){
  d <- as.matrix(read.table(paste0("output/CTL_", diet, "_LODmatrix.txt"), sep = "\t", check.names = FALSE)); m <- real[rownames(d), paste(diet, colnames(d))]
  cat(diet, ": against CTL.R,", sum(!is.na(d)), "values; largest difference", max(abs(d - m), na.rm = TRUE), "; cells missing in one only:", sum(is.na(d) != is.na(m)), "\n")
}
mult <- tabulate(pat, ncol(G))                                                              # markers per pattern, for counts over markers
summarize <- function(x, perm){                                                            # x: patterns x scans, signed -log10 P
  ax <- abs(x); p <- 10^-ax; p[is.na(p)] <- 1; pm <- p[pat, ]
  mx <- apply(ax, 2, max, na.rm = TRUE)
  f <- function(diet){ j <- grep(paste0("^", diet, " "), colnames(x)); q <- pm[, j]
    cs <- apply(q, 2, function(v) sum(p.adjust(v, "BH") <= 0.05)); bf <- sum(q <= 0.05 / nrow(q)); bh.all <- sum(p.adjust(as.vector(q), "BH") <= 0.05); bf.all <- sum(q <= 0.05 / length(q))
    lam <- median(qchisq(q, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
    c(any = max(mx[j]), BH.scans = sum(cs > 0), BH.one.T.markers = sum(cs), BH.all.markers = bh.all, Bonf.one.T.markers = bf, Bonf.all.markers = bf.all, lambda = lam) }
  data.frame(perm = perm, t(mx), CD = t(f("CD")), HF = t(f("HF")), check.names = FALSE)
}
r0 <- summarize(scanAll(G), 0)
t0 <- Sys.time()
ch <- split(1:nperm, rep(1:cores, length.out = nperm))
one <- function(i){ set.seed(20261124 + i)
  if(mode == "mice"){ d2 <- lapply(dat, function(d){ k <- sample.int(length(d$s)); d$s <- d$s[k]; d$wt <- d$wt[k]; d }); return(summarize(scanAll(G, d2), i)) }
  o <- seq_along(st); for(k in unique(stratum)){ j <- which(stratum == k); o[j] <- j[sample.int(length(j))] }; summarize(scanAll(G[o, , drop = FALSE]), i) }
rp <- mclapply(1:nperm, one, mc.cores = cores, mc.preschedule = TRUE)
bad <- sapply(rp, inherits, "try-error"); if(any(bad)) stop(as.character(rp[bad][[1]]))
res <- rbind(r0, do.call(rbind, rp))
colnames(res) <- gsub(" ", ".T", colnames(res))
num <- sapply(res, is.numeric); res[, num] <- round(res[, num], 4)
write.table(res, paste0("output/ctl_permutations", sfx, "_maxima.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("Ran", nperm, "permutations in", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "minutes\n")

sink(paste0("output/ctl_permutations", sfx, ".txt"), split = TRUE)
cat("CTL maps of BXD mice (Danny Arends's method), with", nperm, "permutations, mode:", mode, "\n")
nul <- res[-1, ]; rl <- res[1, ]
for(diet in c("CD", "HF")){
  cat(sprintf("\n==== %s\n", diet)); m <- grep(paste0("^", diet, "\\.T"), colnames(res), value = TRUE)
  cat(sprintf("Inflation (median chi-square / expected): real %.2f; permutations, median %.2f (5%% to 95%%: %.2f to %.2f)\n", rl[, paste0(diet, ".lambda")], median(nul[, paste0(diet, ".lambda")]), quantile(nul[, paste0(diet, ".lambda")], .05), quantile(nul[, paste0(diet, ".lambda")], .95)))
  cat(sprintf("B&H at 5%% FDR over one T-age's markers: real data, %d of %d T-ages call at least one marker (%d marker x T-age calls); in the permutations, at least one call in %.0f%%, median calls %.0f\n", rl[, paste0(diet, ".BH.scans")], length(m), rl[, paste0(diet, ".BH.one.T.markers")], 100 * mean(nul[, paste0(diet, ".BH.scans")] > 0), median(nul[, paste0(diet, ".BH.one.T.markers")])))
  cat(sprintf("B&H over all markers and T-ages: real %d calls; permutations, at least one call in %.0f%% \n", rl[, paste0(diet, ".BH.all.markers")], 100 * mean(nul[, paste0(diet, ".BH.all.markers")] > 0)))
  cat(sprintf("Bonferroni over one T-age's markers: real %d calls; permutations, at least one call in %.0f%%\n", rl[, paste0(diet, ".Bonf.one.T.markers")], 100 * mean(nul[, paste0(diet, ".Bonf.one.T.markers")] > 0)))
  cat(sprintf("Bonferroni over all markers and T-ages: real %d calls; permutations, at least one call in %.0f%%\n", rl[, paste0(diet, ".Bonf.all.markers")], 100 * mean(nul[, paste0(diet, ".Bonf.all.markers")] > 0)))
  an <- paste0(diet, ".any"); cat(sprintf("Any marker, any T-age: largest real %.2f; permutation thresholds 5%% %.2f, 10%% %.2f, 63%% %.2f; genome-wide P = %.3f\n", rl[, an], quantile(nul[, an], .95), quantile(nul[, an], .9), quantile(nul[, an], .37), (1 + sum(nul[, an] >= rl[, an])) / (nperm + 1)))
  tab <- data.frame(T = sub(".*\\.T", "", m), real = as.numeric(rl[, m]), P05 = apply(nul[, m], 2, quantile, .95), P10 = apply(nul[, m], 2, quantile, .9), P63 = apply(nul[, m], 2, quantile, .37), P.genomewide = sapply(m, function(k) (1 + sum(nul[, k] >= rl[, k])) / (nperm + 1)))
  jj <- sapply(m, function(k){ x <- real[, sub("\\.T", " ", k)]; i <- which.max(abs(x)); sprintf("Chr %s %.1f Mb %s", map[i, "Chr"], map[i, "Mb"], if(x[i] > 0) "r(D) > r(B)" else "r(B) > r(D)") })
  tab$peak <- jj; print(cbind(tab[1], round(tab[2:5], 2), P.gw = round(tab$P.genomewide, 3), peak = tab$peak), row.names = FALSE)
}
sink()
