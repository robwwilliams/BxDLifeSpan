#
# cauchy_strainmeans.R - One P per marker over all T-ages (the Cauchy combination), for the strain-mean maps with and
# without the F1s, with genome-wide thresholds from permutations
#
# Rob Williams's fork, 23 November 2026. Exploratory. As planned in bxd_mortality_rates_2026/notes/Cauchy_review_2026-11-21.md
# ("Using it in the BXD paper"), and asked for by Rob on 23 November 2026 ("Go ahead and do as many of the open tasks").
#   - The Cauchy combination (Liu & Xie 2020, J Am Stat Assoc 115:393), as in Danny's Cauchy.R for Arends et al. 2026:
#     at each marker the P values of its T-ages are turned into Cauchy variables, tan((0.5 - P) pi), averaged, and turned
#     back into one P. It allows for the strong correlation between neighbouring survivorships without our having to
#     know it. Here P is exact, on 1 d.f. (the maps' -log10 P, unrounded), and kept below 1 - 1e-10 so that none is
#     infinite. Over the stable range of T-ages of each map (significant_loci.R, significant_loci_F1.R).
#   - Genome-wide thresholds: the same permutations as the maps' own (strains among strains, as in qtl2_permutations.R;
#     with the F1s, F1 types within class as well, as in qtl2_strainF1_perm.R), drawn afresh from seed 20261125. Each
#     permutation records the largest Cauchy -log10 P over markers. 1,000 permutations.
#   - Chow with the F1s is also combined over T20-T770, one T-age past its stable range, because the chow Chr 1 locus acts
#     at T770 and the strains-alone map's stable range includes it.
#   - Maps: strain means alone (qtl2_strainmeans.R, equal weights) and strain and F1 means together, additive
#     (qtl2_strainF1.R), on chow (CD) and high fat (HF). Permutation 0 is the real data, as a check: its single-T-age maps
#     match the published ones.
# Writes output/cauchy_strainmeans.txt (thresholds and the loci that pass), output/cauchy_strainmeans_markers.txt (Cauchy
# -log10 P at every marker) and output/cauchy_strainmeans_maxima.txt
#   Rscript cauchy_strainmeans.R [permutations, default 1000] [cores, default 8]
#

source("shared.R")
source("qtl2_cross.R")                                                              # strains, pr
source("strainF1_setup.R")                                                         # units, AP, Kloco, unitMeans, classCovar
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
nperm <- if(length(args) >= 1) as.numeric(args[1]) else 1000
cores <- if(length(args) >= 2) as.numeric(args[2]) else 8
tps <- seq(20, 770, 30)
lastT <- list(strains = c(CD = 770, HF = 620), F1 = c(CD = 740, HF = 590), F1.T770 = c(CD = 770))   # stable ranges; and chow with F1s to T770,
                                                                                    # one T-age past its stable range, to compare with strains alone
lp <- function(lod) -pchisq(2 * log(10) * pmax(lod, 0), 1, lower.tail = FALSE, log.p = TRUE) / log(10)
cauchy <- function(L){ p <- pmin(10^-L, 1 - 1e-10); -log10(pcauchy(rowMeans(tan((0.5 - p) * pi)), lower.tail = FALSE)) }

set.seed(20261125)
maxima <- NULL; markers <- NULL; check <- NULL
for(set in names(lastT)) for(diet in names(lastT[[set]])){
  keep <- paste0("T", tps[tps <= lastT[[set]][diet]])
  Y <- unitMeans(diet, tps, strains.only = set == "strains")[, keep, drop = FALSE]; Y <- Y[rowSums(!is.na(Y)) > 0, , drop = FALSE]
  u <- rownames(Y); X <- if(set != "strains") classCovar(u) else NULL
  ap <- AP[u, ]; Ku <- lapply(Kloco, function(k) k[u, u])
  groups <- if(set != "strains") list(which(u %in% strains), which(u %in% types[f1class == "short"]), which(u %in% types[f1class == "long"])) else list(seq_along(u))
  perms <- lapply(1:nperm, function(p){ o <- seq_along(u); for(g in groups) o[g] <- g[sample.int(length(g))]; o })
  out <- mclapply(0:nperm, function(p){
    Yp <- if(p == 0) Y else Y[perms[[p]], , drop = FALSE]; rownames(Yp) <- u
    L <- lp(scan1(ap, Yp, kinship = Ku, addcovar = X))
    if(p == 0) list(c = cauchy(L), L = L) else max(cauchy(L))
  }, mc.cores = cores, mc.preschedule = FALSE)
  maxima <- rbind(maxima, data.frame(set = set, diet = diet, perm = 1:nperm, max.cauchy = unlist(out[-1])))
  markers <- cbind(markers, out[[1]]$c)
  colnames(markers)[ncol(markers)] <- paste(set, diet, sep = ".")
  f <- if(set == "strains") paste0("output/qtl2_", diet, "_means_equal_Progressive.txt") else paste0("output/qtl2_", diet, "_meansF1_add_Progressive.txt")
  pub <- abs(as.matrix(read.table(f, sep = "\t")[rownames(out[[1]]$L), sub("T", "X", keep)]))
  check <- c(check, max(abs(pub - out[[1]]$L), na.rm = TRUE))
  cat(set, diet, ": done\n")
}
write.table(maxima, "output/cauchy_strainmeans_maxima.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cbind(map[rownames(markers), c("Chr", "Mb")], round(markers, 4)), "output/cauchy_strainmeans_markers.txt", sep = "\t", quote = FALSE)

levels <- c(`5%` = 0.95, `10%` = 0.90, `suggestive (P = 0.63)` = 0.37)
sink("output/cauchy_strainmeans.txt", split = TRUE)
cat("Cauchy combination over T-ages (one -log10 P per marker) for the strain-mean maps, with genome-wide thresholds from", nperm, "permutations\n")
cat(sprintf("Check: the real data reproduce the published single-T-age maps to %.1e (largest difference in -log10 P)\n", max(check)))
for(set in names(lastT)) for(diet in names(lastT[[set]])){
  nmc <- paste(set, diet, sep = "."); m <- maxima[maxima$set == set & maxima$diet == diet, "max.cauchy"]
  th <- quantile(m, levels); cm <- markers[, nmc]
  cat(sprintf("\n%s, %s (T20 to T%d): thresholds %s\n", c(strains = "Strain means alone", F1 = "Strain and F1 means, additive", F1.T770 = "Strain and F1 means, additive, past the stable range")[set], diet, lastT[[set]][diet],
              paste(sprintf("%s %.2f", names(levels), th), collapse = "; ")))
  cat("  Strongest marker on each chromosome reaching the suggestive level:\n")
  for(ch in c(1:19, "X")){
    k <- which(map[names(cm), "Chr"] == ch); if(length(k) == 0 || max(cm[k]) < th[3]) next
    w <- k[which.max(cm[k])]
    cat(sprintf("    Chr %-2s %6.1f Mb: Cauchy -log10 P %.2f (%s)\n", ch, map[names(cm)[w], "Mb"], cm[w],
                if(cm[w] >= th[1]) "significant" else if(cm[w] >= th[2]) "P < 0.10" else "suggestive"))
  }
}
sink()
