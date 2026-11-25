#
# harmony_F1_models.R - Do the mixed models with F1s agree with plainer models?
#
# Rob Williams's fork, 24 November 2026. Exploratory. Written by Claude (Opus 5.5). Rob: "As long as any new code gives results that are
# generally harmonious with less complex models, then all is good." Compared, over all markers, on the signed -log10 P of the D allele:
#   A  strain + F1 means: R/qtl2 mixed model with LOCO kinship (qtl2_strainF1.R) against plain lm(mean ~ F1 class + D dosage), no
#      kinship; and against the strains alone in the same two models
#   B  individual mice: R/qtl2 (qtl2_individualsF1.R) against lme4's lmer(age ~ vivarium + set-up class + F1 class + D dosage + (1 | strain)),
#      as in lifespan.R, on 300 markers spread over the genome (lmer is slow)
# Writes output/harmony_F1_models.txt
#   Rscript harmony_F1_models.R
#

source("shared.R"); source("qtl2_cross.R"); source("strainF1_setup.R")
library(lme4); set.seed(20261124)
tps <- c(20, 380, 500, 620)
dose <- do.call(cbind, lapply(AP, function(a) a[, "D", ]))                                      # units x markers, D dosage 0 to 1
lpl <- function(y, X, g){                                                                       # plain lm, signed -log10 P of the dosage
  X <- cbind(rep(1, length(y)), X); q <- qr(X); ry <- qr.resid(q, y); rg <- qr.resid(q, g); r <- colSums(rg * ry) / sqrt(colSums(rg^2) * sum(ry^2)); df <- length(y) - ncol(X) - 1
  sign(r) * -log10(2 * pt(-abs(r * sqrt(df / (1 - r^2))), df)) }
rd <- function(f, T){ x <- read.table(file.path("output", f), sep = "\t"); setNames(as.numeric(x[, paste0("X", T)]), rownames(x)) }
out <- c(); add <- function(what, diet, T, r, n, pk) out <<- rbind(out, data.frame(comparison = what, diet = diet, T = T, r = round(r, 3), n = n, peaks = pk))
peaks <- function(a, b){ ch <- paste0("Chr", map[names(a), "Chr"], "@"); top <- function(x) names(sort(tapply(abs(x), ch, max), decreasing = TRUE))[1:5]; length(intersect(top(a), top(b))) }
for(d in c("CD", "HF")){
  Y <- unitMeans(d, tps); Ys <- unitMeans(d, tps, strains.only = TRUE)
  for(T in tps){
    k <- paste0("T", T); u <- rownames(Y)[!is.na(Y[, k])]; s <- rownames(Ys)[!is.na(Ys[, k])]
    m <- rd(paste0("qtl2_", d, "_meansF1_add_Progressive.txt"), T); ms <- rd(paste0("qtl2_", d, "_means_equal_Progressive.txt"), T)
    p <- lpl(Y[u, k], classCovar(u), dose[u, names(m)]); ps <- lpl(Ys[s, k], NULL, dose[s, names(m)])
    ok <- is.finite(m) & is.finite(p) & is.finite(ms[names(m)]) & is.finite(ps)
    add("strain + F1 means: mixed model vs plain lm with F1 class", d, T, cor(m[ok], p[ok]), length(u), peaks(m[ok], p[ok]))
    add("strains alone: mixed model vs plain lm", d, T, cor(ms[names(m)][ok], ps[ok]), length(s), peaks(ms[names(m)][ok], ps[ok]))
    add("mixed model: strains + F1s vs strains alone", d, T, cor(m[ok], ms[names(m)][ok]), length(u), peaks(m[ok], ms[names(m)][ok]))
  }
}
# Individual mice, 300 markers
f <- f1; x <- ldata[which(ldata[, "Diet"] %in% c("CD", "HF") & grepl("BXD", ldata[, "StrainName"])), ]
mi <- rbind(data.frame(strain = as.character(x[, "StrainName"]), diet = x[, "Diet"], age = x[, "AgeAtDeath..days."], viv = as.character(x[, "vivarium"]), setup = as.numeric(x[, "AgeAtSetUp.in.colony..days."]), cls = "inbred"),
            data.frame(strain = f[, "StrainName"], diet = f[, "DietCode"], age = f[, "AgeAtDeath (days)"], viv = "3", setup = f[, "AgeAtSet Up (days)"], cls = unname(f1class[f[, "StrainName"]])))
mi <- mi[mi$strain %in% units & complete.cases(mi) & mi$setup >= 0, ]; mi$su <- cut(mi$setup, seq(0, 1000, 180), include.lowest = TRUE); mi$cls <- factor(mi$cls, levels = c("inbred", "short", "long"))
mk <- sample(colnames(dose), 300)
for(d in c("CD", "HF")) for(T in c(20, 500)){
  z <- droplevels(mi[mi$diet == d & mi$age >= T, ]); full <- rd(paste0("qtl2_", d, "_indF1_Progressive.txt"), T)
  f0 <- lmer(age ~ viv + su + cls + (1 | strain), z, REML = FALSE)
  v <- suppressMessages(sapply(mk, function(k){ z$g <- dose[z$strain, k]; f1m <- update(f0, . ~ . + g, data = z); sign(fixef(f1m)["g"]) * -log10(pchisq(2 * (logLik(f1m) - logLik(f0)), 1, lower.tail = FALSE)) }))
  a <- full[mk]; ok <- is.finite(a) & is.finite(v)
  add("individual mice: R/qtl2 with LOCO kinship vs lmer with a random strain intercept (300 markers)", d, T, cor(a[ok], v[ok]), nrow(z), NA)
}
options(width = 200); print(out, row.names = FALSE)
write.table(out, "output/harmony_F1_models.txt", sep = "\t", quote = FALSE, row.names = FALSE)
