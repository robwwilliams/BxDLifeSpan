#
# founder_effects_joint.R - The effect of each founder allele (C, B, H, D) along proximal Chr 1 and distal Chr 2, from BXD
# and UM-HET3 together
#
# Rob Williams's fork, 24 November 2026. Exploratory, after both families' maps were seen. Rob: "Just do them both" (this,
# and NAM by Vita1 genotype).
#   - Three contrasts, one per family and parent: BXD D - B (a strain-mean difference between homozygotes, so 2 (a_D - a_B)
#     for per-allele effects a), UM-HET3 paternal D - H (a_D - a_H, per allele) and maternal B - C (a_B - a_C). Three
#     contrasts fix four founder effects up to a constant: with a_B = 0, a_D = BXD / 2, a_H = a_D - paternal, a_C = -maternal.
#     The founder ranking is then signed, which Fisher's combination of P values is not. It assumes additivity and the
#     same effect in both families' backgrounds, and nothing tests that here: three contrasts leave no degree of freedom.
#   - BXD: females on chow, strain means (at least 3 mice alive at T, equal weights), a mixed model with kinship leaving
#     out the chromosome (REML heritability without the locus, then generalized least squares); the coefficient of D dosage
#     is DD - BB, with its standard error. At the BXD marker nearest each UM-HET3 marker.
#   - UM-HET3: females, lm(age at death ~ site + cohort + treatment + dosage), per-allele coefficients with standard
#     errors, at the 891 markers of Arends et al. 2026.
#   - Three life stages, rows matched on the 40-day grid (umhet3_maps.R): BXD T110 with UM-HET3 T125, T410 with T525, and
#     T740 with T965. Each family's own T-ages are used for its own fit.
# Writes output/founder_effects_joint.txt and output/maps/Founder_effects_joint.svg
#   Rscript founder_effects_joint.R
#

library(qtl)
library(svglite)
source("shared.R")
source("qtl2_cross.R")                                                              # strains, pr
source("strainF1_setup.R")                                                         # AP, Kloco, unitMeans
stages <- data.frame(name = c("early", "middle", "late"), bxd = c(110, 410, 740), um = c(125, 525, 965))
regions <- list("1" = c(0, 45), "2" = c(100, 182))
cB <- "#7d4598"; cD <- "#ed8f22"

### UM-HET3 females, as in umhet3_maps.R
dir <- "../bxd_mortality_rates_2026/data/local/Arends_Supplementary_Tables_10Apr2026"
mc <- read.cross(format = "csvr", dir = dir, file = "S_Table_7_R:qtl_v1_cross-object_for_mapping_18Aug2025.csvr", genotypes = NULL, na.strings = c("-", "NA"))
mc <- calc.genoprob(mc, step = 0); gp <- pull.genoprob(mc)
um.markers <- read.delim("output/umhet3_paternal_female_Progressive.txt", check.names = FALSE)
um.markers <- um.markers[um.markers$Chr %in% names(regions), c("Chr", "Mb")]
um.markers <- um.markers[mapply(function(c, mb) mb >= regions[[c]][1] & mb <= regions[[c]][2], um.markers$Chr, um.markers$Mb), ]
ph <- pull.pheno(mc); fem <- which(as.character(ph[, "sex"]) == "0")
cd <- data.frame(age = as.numeric(ph[fem, "longevity"]), site = factor(ph[fem, "site"]), cohort = factor(ph[fem, "cohort"]), treatment = factor(ph[fem, "treatment"]))
dos <- function(m, cols) rowSums(gp[fem, paste0(m, ":", cols), drop = FALSE])     # R/qtl 4-way: AC = CH, BC = BH, AD = CD, BD
umFit <- function(T){
  ok <- which(!is.na(cd$age) & cd$age >= T); d <- droplevels(cd[ok, ])
  X <- model.matrix(~ site + cohort + treatment, d); qx <- qr(X); r0 <- qr.resid(qx, d$age); n <- length(ok); p <- qx$rank
  t(sapply(rownames(um.markers), function(m){
    z <- sapply(list(pat = c("AD", "BD"), mat = c("BC", "BD")), function(cl){ gr <- qr.resid(qx, dos(m, cl)[ok]); gg <- sum(gr^2); b <- sum(gr * r0) / gg
      s2 <- (sum(r0^2) - b^2 * gg) / (n - p - 1); c(b, sqrt(s2 / gg)) })
    setNames(c(z[1, "pat"], z[2, "pat"], z[1, "mat"], z[2, "mat"]), c("pat.b", "pat.se", "mat.b", "mat.se"))
  }))
}

### BXD strain means on chow, GLS with LOCO kinship
remlH <- function(yr, Xr, d){ f <- function(h){ v <- h * d + 1 - h; w <- 1 / sqrt(v); q <- qr(Xr * w); r <- qr.resid(q, yr * w)
  -0.5 * (sum(log(v)) + 2 * sum(log(abs(diag(qr.R(q))))) + (length(yr) - q$rank) * log(sum(r^2))) }; optimize(f, c(0, 0.999), maximum = TRUE)$maximum }
bxdFit <- function(T, chr, mks){
  Y <- unitMeans("CD", T, strains.only = TRUE); u <- rownames(Y)[!is.na(Y[, 1])]; y <- Y[u, 1]
  e <- eigen(Kloco[[chr]][u, u], symmetric = TRUE); d <- pmax(e$values, 0); U <- e$vectors
  yr <- drop(crossprod(U, y)); X0 <- crossprod(U, rep(1, length(u))); h <- remlH(yr, X0, d); w <- 1 / sqrt(h * d + 1 - h)
  t(sapply(mks, function(k){ x <- drop(crossprod(U, AP[[chr]][u, "D", k])); Xw <- cbind(X0, x) * w; f <- lm.fit(Xw, yr * w)
    s2 <- sum(f$residuals^2) / (length(u) - 2); V <- s2 * solve(crossprod(Xw)); c(b = unname(f$coefficients[2]), se = sqrt(V[2, 2])) }))
}

res <- NULL
for(s in seq_len(nrow(stages))){
  uf <- umFit(stages$um[s])
  for(chr in names(regions)){
    um <- um.markers[um.markers$Chr == chr, ]; bk <- dimnames(AP[[chr]])[[3]]; bmb <- map[bk, "Mb"]
    near <- sapply(um$Mb, function(mb) bk[which.min(abs(bmb - mb))])
    bf <- bxdFit(stages$bxd[s], chr, near)
    r <- data.frame(stage = stages$name[s], chr = chr, Mb = um$Mb, bxd = bf[, "b"], bxd.se = bf[, "se"],
                    pat = uf[rownames(um), "pat.b"], pat.se = uf[rownames(um), "pat.se"], mat = uf[rownames(um), "mat.b"], mat.se = uf[rownames(um), "mat.se"])
    r$a.D <- r$bxd / 2; r$a.H <- r$a.D - r$pat; r$a.C <- -r$mat                   # per-allele effects, days, relative to B
    r$se.D <- r$bxd.se / 2; r$se.H <- sqrt(r$se.D^2 + r$pat.se^2); r$se.C <- r$mat.se
    res <- rbind(res, r)
  }
}
write.table(cbind(res[, 1:3], round(res[, -(1:3)], 1)), "output/founder_effects_joint.txt", sep = "\t", quote = FALSE, row.names = FALSE)

### Summary at the Vita loci
sink("output/founder_effects_joint_summary.txt", split = TRUE)
cat("Per-allele founder effects on lifespan (days, relative to B = C57BL/6), from BXD D - B and UM-HET3 paternal D - H and maternal\n")
cat("B - C (females; BXD on chow). Mean over the UM-HET3 markers in each interval; +/- 1 s.e. of that mean's components.\n")
iv <- data.frame(name = c("Vita1a (3-8 Mb)", "Vita1b (22-26 Mb)", "Vita1y (34-38 Mb)", "Mrps5 locus (124-129 Mb)", "UM-HET3 Chr 2 peak (140-153 Mb)", "Vita2b (153-159 Mb)"),
                 chr = c("1", "1", "1", "2", "2", "2"), from = c(3, 22, 34, 124, 140, 153), to = c(8, 26, 38, 129, 153, 159))
for(s in stages$name){ cat("\n", s, " life (BXD T", stages$bxd[stages$name == s], ", UM-HET3 T", stages$um[stages$name == s], ")\n", sep = "")
  for(i in seq_len(nrow(iv))){ z <- res[res$stage == s & res$chr == iv$chr[i] & res$Mb >= iv$from[i] & res$Mb <= iv$to[i], ]
    if(!nrow(z)) next
    o <- sort(c(B = 0, D = mean(z$a.D), H = mean(z$a.H), C = mean(z$a.C)), decreasing = TRUE)
    cat(sprintf("  %-32s D %+5.0f (%.0f), H %+5.0f (%.0f), C %+5.0f (%.0f); ranking %s  [%d markers]\n", iv$name[i], mean(z$a.D), mean(z$se.D), mean(z$a.H), mean(z$se.H),
                mean(z$a.C), mean(z$se.C), paste(names(o), collapse = " > "), nrow(z))) } }
sink()

### Figure
svglite("output/maps/Founder_effects_joint.svg", width = 7.2, height = 7.2, bg = "white", system_fonts = list(sans = "Arial"), fix_text_size = FALSE)
par(mfrow = c(3, 2), mai = c(0.35, 0.6, 0.3, 0.1), oma = c(1.6, 0, 0.6, 0), cex = 0.62, mgp = c(2, 0.5, 0), tcl = -0.25, las = 1)
for(s in stages$name) for(chr in names(regions)){
  z <- res[res$stage == s & res$chr == chr, ]; z <- z[order(z$Mb), ]
  yl <- range(c(0, z$a.D - z$se.D, z$a.H - z$se.H, z$a.C - z$se.C, z$a.D + z$se.D, z$a.H + z$se.H, z$a.C + z$se.C), na.rm = TRUE)
  plot(NA, xlim = regions[[chr]], ylim = yl, xaxs = "i", bty = "l", xlab = "", ylab = "days per allele, vs B")
  abline(h = 0, col = cB, lwd = 1.2)
  for(f in c("C", "H", "D")){ co <- c(C = "grey55", H = "black", D = cD)[f]
    polygon(c(z$Mb, rev(z$Mb)), c(z[[paste0("a.", f)]] - z[[paste0("se.", f)]], rev(z[[paste0("a.", f)]] + z[[paste0("se.", f)]])), col = adjustcolor(co, 0.15), border = NA)
    lines(z$Mb, z[[paste0("a.", f)]], col = co, lwd = 1.1, type = "o", pch = 16, cex = 0.35) }
  mtext(sprintf("Chr %s, %s life (BXD T%d, UM-HET3 T%d)", chr, s, stages$bxd[stages$name == s], stages$um[stages$name == s]), side = 3, line = 0.2, adj = 0, cex = 0.7, font = 2)
  if(s == "early" && chr == "1") legend("topright", bty = "n", cex = 0.8, lwd = 1.1, col = c(cD, "black", "grey55", cB), legend = c("D (DBA/2J)", "H (C3H/HeJ)", "C (BALB/cByJ)", "B (C57BL/6), reference"))
}
mtext("Chr 1 position (Mb)", side = 1, outer = TRUE, line = 0.4, at = 0.27, cex = 0.7); mtext("Chr 2 position (Mb)", side = 1, outer = TRUE, line = 0.4, at = 0.77, cex = 0.7)
invisible(dev.off())
