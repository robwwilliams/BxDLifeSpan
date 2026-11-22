#
# f1_heterosis.R - Heterosis of lifespan in each BXD x BXD F1 type: the F1 mean against the mean of its parents
#
# Rob Williams's fork, 21 November 2026. Exploratory. Rob asked whether a categorical cofactor (short x short against
# long x long) is enough to correct the strain + F1 maps (qtl2_strainF1.R) for heterosis. This tabulates, for each F1 type
# and diet: the parents' strain means and their midparent, the F1 mean, heterosis (F1 - midparent), and the fraction of
# the genome at which the parents differ (so the F1 is heterozygous); also heterosis as a percent of the midparent, and the
# F1 mean minus its better parent's. All mice of each group (T-age 20). Within a class, a
# heterogeneity test asks whether heterosis differs among F1 types by more than sampling error: Q = sum((h - mean h)^2 / SE^2)
# on (types - 1) d.f., with SE^2 = var(F1) / n(F1) + (var(parent 1) / n + var(parent 2) / n) / 4 and the mean weighted by 1 / SE^2.
#   Rscript f1_heterosis.R [path to the F1 table]
# Writes output/f1_heterosis.txt
#

source("shared.R")
args <- commandArgs(trailingOnly = TRUE)
f1file <- if(length(args) >= 1) args[1] else "../bxd_mortality_rates_2026/data/corrected/AgingBXD_F1_BodyWeight_NaturalDeath_corrected.txt"
short <- c("BXD9", "BXD170", "BXD172", "BXD197")
f1 <- read.delim(f1file, check.names = FALSE)
f1 <- f1[f1[, "CauseOfDeath"] == "Natural" & f1[, "DietCode"] %in% c("CD", "HF") & !is.na(f1[, "AgeAtDeath (days)"]), ]
types <- sort(unique(f1[, "StrainName"]))
p1 <- sub("x.*", "", types); p2 <- paste0("BXD", sub(".*x", "", types))
het <- sapply(seq_along(types), function(i){ a <- geno[, p1[i]]; b <- geno[, p2[i]]; ok <- a %in% c("B", "D") & b %in% c("B", "D"); mean(a[ok] != b[ok]) })

out <- NULL
for(d in c("CD", "HF")){
  a <- ldata[ldata[, "Diet"] == d, ]
  m <- tapply(a[, "AgeAtDeath..days."], a[, "StrainName"], mean); v <- tapply(a[, "AgeAtDeath..days."], a[, "StrainName"], var); n <- table(a[, "StrainName"])
  z <- f1[f1[, "DietCode"] == d, ]
  fm <- tapply(z[, "AgeAtDeath (days)"], z[, "StrainName"], mean)[types]; fn <- table(z[, "StrainName"])[types]
  fv <- tapply(z[, "AgeAtDeath (days)"], z[, "StrainName"], var)[types]
  se <- sqrt(fv / fn + (v[p1] / n[p1] + v[p2] / n[p2]) / 4)
  out <- rbind(out, data.frame(diet = d, F1 = types, class = ifelse(p1 %in% short, "short x short", "long x long"), mice = as.numeric(fn),
                               parent1 = round(m[p1]), parent2 = round(m[p2]), midparent = round((m[p1] + m[p2]) / 2), F1.mean = round(fm),
                               heterosis = round(fm - (m[p1] + m[p2]) / 2), se = round(as.numeric(se)),
                               percent = round(100 * (fm - (m[p1] + m[p2]) / 2) / ((m[p1] + m[p2]) / 2)), above.better.parent = round(fm - pmax(m[p1], m[p2])),
                               heterozygous = round(het, 3)))
}
sink("output/f1_heterosis.txt", split = TRUE)
cat("Heterosis of lifespan (days) in BXD x BXD F1 types: F1 mean minus the mean of its parents' strain means (all mice)\n\n")
print(out, row.names = FALSE)
cat("\nBy class and diet: mean heterosis (SD over F1 types); and the correlation of heterosis with the heterozygous fraction\n")
for(d in c("CD", "HF")) for(k in c("short x short", "long x long")){
  x <- out[out[, "diet"] == d & out[, "class"] == k, ]
  w <- 1 / x[, "se"]^2; hbar <- sum(w * x[, "heterosis"]) / sum(w); Q <- sum(w * (x[, "heterosis"] - hbar)^2)
  cat(sprintf("  %s %-13s %d types: %+.0f days (SD %.0f over types; median SE of one type %.0f), range %+d to %+d; heterogeneity Q = %.1f on %d d.f., P = %.2g; heterozygous %.2f to %.2f; r(heterosis, heterozygous) = %.2f\n",
              d, k, nrow(x), mean(x[, "heterosis"]), sd(x[, "heterosis"]), median(x[, "se"]), min(x[, "heterosis"]), max(x[, "heterosis"]), Q, nrow(x) - 1,
              pchisq(Q, nrow(x) - 1, lower.tail = FALSE), min(x[, "heterozygous"]), max(x[, "heterozygous"]), cor(x[, "heterosis"], x[, "heterozygous"])))
}
sink()
