#
# chr1_f1_dropone.R - How much each BXD x BXD F1 type moves the proximal Chr 1 peak and its 1.5-LOD support interval
#
# Rob Williams's fork, 23 November 2026. Exploratory. A check on chr1_zoom.R: there the F1s narrow the chow interval at
# T770 from 4.9-10.4 Mb to 4.9-5.5 Mb, but the edge is set by a plateau (5.5-10.4 Mb) that lies close to 1.5 LOD below the
# peak. So the strain + F1 scan of Chr 1 (as in qtl2_strainF1.R: unit means with at least 3 mice alive at T, equal
# weights, F1 type as a cofactor, additive genotype, LOCO kinship over strains and F1s) is repeated with each F1 type
# left out in turn, and then each class of F1 (short x short, long x long) left out, at the T-ages of the profiles in
# chr1_zoom.R.
# Writes output/chr1_f1_dropone.txt
#   Rscript chr1_f1_dropone.R [path to the F1 table]
#

source("shared.R")
source("qtl2_cross.R")                                                              # strains, pr and the map

args <- commandArgs(trailingOnly = TRUE)
f1file <- if(length(args) >= 1) args[1] else "../bxd_mortality_rates_2026/data/corrected/AgingBXD_F1_24Nov2026_corrected.txt"
zoomMb <- 50
views <- list(c("CD", 770), c("HF", 20), c("HF", 260), c("HF", 500))
source("strainF1_setup.R")                                                         # types, f1class, AP and Kloco, as in qtl2_strainF1.R
K1 <- Kloco[["1"]]

k <- dimnames(AP[["1"]])[[3]]; k <- k[map[k, "Mb"] <= zoomMb]
mb <- as.numeric(map[k, "Mb"])
lod2p <- function(L) -pchisq(2 * log(10) * pmax(L, 0), 1, lower.tail = FALSE, log.p = TRUE) / log(10)

### One Chr 1 scan: the peak in 0-50 Mb and its contiguous 1.5-LOD interval, as in chr1_zoom.R
scanChr1 <- function(diet, T, drop = NULL){
  x <- ld[which(ld[, "Diet"] == diet & ld[, "AgeAtDeath..days."] >= T), ]
  n <- table(x[, "StrainName"]); s <- names(n)[n >= 3]
  y <- tapply(x[, "AgeAtDeath..days."], x[, "StrainName"], mean)[s]
  z <- f1[f1[, "DietCode"] == diet & f1[, "AgeAtDeath (days)"] >= T & !f1[, "StrainName"] %in% drop, ]
  n1 <- table(z[, "StrainName"]); t1 <- names(n1)[n1 >= 3]
  y1 <- tapply(z[, "AgeAtDeath (days)"], z[, "StrainName"], mean)[t1]
  u <- c(s, t1)
  Y <- matrix(c(y, y1), dimnames = list(u, "lifespan"))
  cls <- c(rep("inbred", length(s)), f1class[t1])
  X <- cbind(short = as.numeric(cls == "short"), long = as.numeric(cls == "long")); rownames(X) <- u
  X <- X[, colSums(X) > 0, drop = FALSE]; if(ncol(X) == 0) X <- NULL
  ap <- AP[u, "1"]; ap[["1"]] <- ap[["1"]][, , k, drop = FALSE]
  L <- scan1(ap, Y, kinship = decomp_kinship(K1[u, u]), addcovar = X)[, 1]
  i <- which.max(L); cut <- L[i] - 1.5
  lo <- i; while(lo > 1 && L[lo - 1] >= cut) lo <- lo - 1
  hi <- i; while(hi < length(L) && L[hi + 1] >= cut) hi <- hi + 1
  data.frame(diet = diet, T = T, left.out = if(is.null(drop)) "none" else if(length(drop) > 1) paste("all", f1class[drop[1]], "x", f1class[drop[1]]) else drop, F1 = length(t1), peak = round(mb[i], 2), p = round(lod2p(L[i]), 2),
             lo = round(mb[lo], 2), hi = round(mb[hi], 2), width = round(mb[hi] - mb[lo], 2), other = if(any(L[-(lo:hi)] >= cut)) "yes" else "no")
}

r <- do.call(rbind, lapply(views, function(v){
  d <- v[1]; T <- as.numeric(v[2])
  rbind(scanChr1(d, T), do.call(rbind, lapply(types, function(t) scanChr1(d, T, t))),
        scanChr1(d, T, types[f1class == "short"]), scanChr1(d, T, types[f1class == "long"]))
}))
sink("output/chr1_f1_dropone.txt", split = TRUE)
cat("Proximal Chr 1 (0 to", zoomMb, "Mb), strain and F1 means together (additive): the peak and its contiguous 1.5-LOD support\n")
cat("interval with all F1 types (left.out none), then with each F1 type left out in turn. p: -log10 P at the peak. F1: F1 types\n")
cat("in the scan. other: another marker in 0-50 Mb, outside the interval, also within 1.5 LOD of the peak.\n\n")
print(r, row.names = FALSE)
cat("\nSummary: the interval with all F1 types; its range with one F1 type left out; and with each class left out\n")
for(v in views){
  x <- r[r$diet == v[1] & r$T == as.numeric(v[2]), ]; a <- x[x$left.out == "none", ]
  o <- x[x$left.out %in% types, ]; cl <- x[grepl("^all", x$left.out), ]
  wider <- o$left.out[o$width > a$width + 1]
  cat(sprintf("  %s T%s: %.1f-%.1f Mb (%.1f Mb), peak %.1f Mb. One type left out: width %.1f to %.1f Mb, peak at %s Mb; more than 1 Mb wider without %s.",
              v[1], v[2], a$lo, a$hi, a$width, a$peak, min(o$width), max(o$width), paste(sort(unique(round(o$peak, 1))), collapse = " or "),
              if(length(wider) > 0) paste(wider, collapse = ", ") else "none"))
  cat(sprintf(" %s: %.1f-%.1f Mb (%.1f Mb).", cl$left.out, cl$lo, cl$hi, cl$width), "\n")
}
sink()
