#
# umhet3_maps.R - Actuarial maps of lifespan in UM-HET3 females (or males), drawn exactly as the BXD maps: a maternal map (B against
# C) and a paternal map (D against H), with genome-wide thresholds from permutations
#
# Rob Williams's fork, 24 November 2026. Exploratory. Rob: "The fun next thing is to make corresponding figures (exact) for
# the UM-HET3 to the extent that this is possible so that we can evaluate location and effects of loci. For the BXDs we
# have two ranges of colors for B and D effects." And, choosing between two designs: "Yeh, let's do only option 1 for
# now" (option 1: maternal and paternal maps in the BXD geometry; option 2, a projection of D - B through founder
# sequence, is left for pangenome mapping of the BXDs).
#   - Females first (Rob, 24 November 2026: "Yep, females only. Results from Arends et al. definitely make this an
#     uncontroversial first-pass decision. Of course it would be fun to see the stratified male data too."); males with
#     the third argument. The X: a female's paternal X is always H (her father's mother was C3H/HeJ) and a male has no
#     paternal X, so the paternal X is not mapped; the maternal X is.
#   - Mice: the UM-HET3 females of Arends et al. 2026 (Nature 654:437; Supplementary Table 7, the R/qtl cross), at the
#     891 markers of its maps (Supplementary Table 9). Genotype probabilities as in Danny's code (calc.genoprob, step 0;
#     adjustXprobs). R/qtl 4-way codes: CH, BH, CD, BD (maternal C or B; paternal H or D).
#   - The two contrasts that UM-HET3 has (no mouse is BB or DD, so there is no D - B contrast):
#       maternal  B - C: dosage P(BH) + P(BD)
#       paternal  D - H: dosage P(CD) + P(BD)
#     Each is its own scan, 1 d.f., signed by its effect: lm(age at death ~ site + cohort + treatment + dosage), as in
#     Danny's ProgressiveMapping.R, without kinship (UM-HET3 mice are individuals of a segregating cross, not strains);
#     LOD = (n/2) log10(RSS0/RSS1); -log10 P from chi-square on 1 d.f. Treatment is a covariate, as in the paper: drug-
#     treated mice are kept.
#   - T-ages: T20 to T1100 in 30-day steps, so that rows T50 to T770 are the BXD maps' rows.
#   - Colors, so that the BXD colors keep their meaning: in the paternal map, warm where D lives longer (as in BXD) and
#     gray where H does; in the maternal map, cool where B lives longer (as in BXD) and gray where C does. The scale has
#     seven bins, to >8, since UM-HET3 peaks reach 9 (Rob: "We will need more dynamic range on the P value colors"); its first
#     four bins are the BXD maps'.
#   - Genome-wide thresholds: 1,000 permutations of lifespans within site (Danny's stratifiedPerms.R, one sex at a time),
#     drawn from seed 20261126. Each records the largest -log10 P over markers at each T-age, for each map, and the
#     largest Cauchy combination over T-ages (T20-T1100 and T20-T770), as for BXD.
#   - Caution: the paternal calls on proximal Chr 1 are biased toward H at UM and UT (UMHET3/chr1_marker_check.R in
#     bxd_mortality_rates_2026), which touches Vita1a and Vita1b on the paternal map.
# Inputs: ../bxd_mortality_rates_2026/data/local/Arends_Supplementary_Tables_10Apr2026 (Supplementary Tables 7 and 9; not in git)
# Writes output/umhet3_<maternal|paternal>_<sex>_Progressive.txt, output/umhet3_N_<sex>.txt, output/umhet3_perm_maxima_<sex>.txt,
# output/umhet3_thresholds_<sex>.txt and output/maps/UMHET3_<sex>_*.svg
#   Rscript umhet3_maps.R [permutations, default 1000] [cores, default 8] [female (default) | male]
#

library(qtl)
library(readxl)
library(parallel)
library(svglite)
source("shared.R")                                                                 # the BXD map, for the chromosome layout
source("actuary_map.R")

args <- commandArgs(trailingOnly = TRUE)
nperm <- if(length(args) >= 1) as.numeric(args[1]) else 1000
cores <- if(length(args) >= 2) as.numeric(args[2]) else 8
sex <- if(length(args) >= 3) args[3] else "female"
Sex <- c(female = "females", male = "males")[sex]
dir <- "../bxd_mortality_rates_2026/data/local/Arends_Supplementary_Tables_10Apr2026"
tps <- seq(20, 1100, 30)

adjustXprobs <- function(cross){                                                    # Danny's, from ActuarialMapping/adjustXprobs.R
  sex <- getsex(cross)$sex
  pr <- cross$geno[["X"]]$prob
  for(i in 1:ncol(pr)){
    pr[sex == 0, i, 3:4] <- 0
    pr[sex == 1, i, 1:2] <- 0
    pr[, i, ] <- pr[, i, ] / rowSums(pr[, i, ])
  }
  cross$geno[["X"]]$prob <- pr
  invisible(cross)
}
mcross <- read.cross(format = "csvr", dir = dir, file = "S_Table_7_R:qtl_v1_cross-object_for_mapping_18Aug2025.csvr", genotypes = NULL, na.strings = c("-", "NA"))
mcross <- calc.genoprob(mcross, step = 0)
mcross <- adjustXprobs(mcross)
gtsp <- pull.genoprob(mcross)
f9 <- as.data.frame(suppressMessages(read_excel(file.path(dir, "S_Table_9_Vita_mapping_results_18Aug2025.xlsx"), sheet = "S_Table_9.3_Female", col_names = FALSE)))
i9 <- which(!is.na(suppressWarnings(as.numeric(f9[, 3]))))[1]
markers <- intersect(colnames(pull.geno(mcross)), unlist(f9[i9 - 1, -(1:2)]))       # the 891 markers of the paper's maps
cols <- lapply(markers, function(m) which(sub(":.*", "", colnames(gtsp)) == m))
matB <- sapply(cols, function(i) gtsp[, i[2]] + gtsp[, i[4]])                       # BH + BD
patD <- sapply(cols, function(i) gtsp[, i[3]] + gtsp[, i[4]])                       # CD + BD
colnames(matB) <- colnames(patD) <- markers
ph <- pull.pheno(mcross)
fem <- which(as.character(ph[, "sex"]) == c(female = "0", male = "1")[sex])          # the mice of one sex
cdata <- data.frame(age = as.numeric(ph[fem, "longevity"]), site = as.factor(ph[fem, "site"]), cohort = as.factor(ph[fem, "cohort"]), treatment = as.factor(ph[fem, "treatment"]))
matB <- matB[fem, ]; patD <- patD[fem, ]
ok <- !is.na(cdata$age); cdata <- cdata[ok, ]; matB <- matB[ok, ]; patD <- patD[ok, ]
umap <- data.frame(Chr = sub("_.*", "", markers), Mb = as.numeric(sub(".*_", "", markers)) / 1e6, row.names = markers)
cat(nrow(cdata), " UM-HET3 ", Sex, "; ", length(markers), " markers\n", sep = "")
lp <- function(lod) -pchisq(2 * log(10) * pmax(lod, 0), 1, lower.tail = FALSE, log.p = TRUE) / log(10)

# Signed -log10 P at every marker, for the mice alive at T, with lifespans y
scanT <- function(y, T){
  idx <- which(y >= T); d <- droplevels(cdata[idx, ])
  X <- model.matrix(~ site + cohort + treatment, d); qx <- qr(X)
  r0 <- qr.resid(qx, y[idx]); rss0 <- sum(r0^2); n <- length(idx)
  sapply(list(maternal = matB, paternal = patD), function(G){
    Gr <- qr.resid(qx, G[idx, , drop = FALSE]); gg <- colSums(Gr^2); gy <- colSums(Gr * r0)
    lod <- n / 2 * log10(rss0 / (rss0 - gy^2 / pmax(gg, 1e-12))); lod[gg < 1e-8] <- 0
    sign(gy) * lp(lod)                                                              # positive: B (maternal) or D (paternal) lives longer
  }, simplify = "array")
}
stackAll <- function(y){ s <- lapply(tps, function(T) scanT(y, T)); list(maternal = sapply(s, function(x) x[, "maternal"]), paternal = sapply(s, function(x) x[, "paternal"])) }
cauchy <- function(L){ p <- pmin(10^-abs(L), 1 - 1e-10); -log10(pcauchy(rowMeans(tan((0.5 - p) * pi)), lower.tail = FALSE)) }

### The real maps
real <- stackAll(cdata$age)
N <- sapply(tps, function(T) sum(cdata$age >= T))
for(m in names(real)){
  out <- cbind(umap, real[[m]]); colnames(out) <- c("Chr", "Mb", tps)
  write.table(out, paste0("output/umhet3_", m, "_", sex, "_Progressive.txt"), sep = "\t", quote = FALSE)
}
write.table(data.frame(T = tps, mice = N), paste0("output/umhet3_N_", sex, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)

### Permutations within site
set.seed(20261126)
perms <- lapply(1:nperm, function(p){ o <- seq_len(nrow(cdata)); for(s in levels(cdata$site)){ i <- which(cdata$site == s); o[i] <- i[sample.int(length(i))] }; o })
mx <- do.call(rbind, mclapply(1:nperm, function(p){
  L <- stackAll(cdata$age[perms[[p]]])
  do.call(rbind, lapply(names(L), function(m) data.frame(perm = p, map = m, t(apply(abs(L[[m]]), 2, max)),
                                                         cauchy.T1100 = max(cauchy(L[[m]])), cauchy.T770 = max(cauchy(L[[m]][, tps <= 770])))))
}, mc.cores = cores, mc.preschedule = FALSE))
colnames(mx)[3:(2 + length(tps))] <- paste0("T", tps)
write.table(mx, paste0("output/umhet3_perm_maxima_", sex, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE)

### Thresholds, and the loci that pass
levels <- c(significant = 0.95, `P 0.10` = 0.90, suggestive = 0.37)
th <- list()
sink(paste0("output/umhet3_thresholds_", sex, ".txt"), split = TRUE)
cat(paste0("UM-HET3 ", Sex, ","), "maternal (B - C) and paternal (D - H) actuarial maps: genome-wide thresholds (-log10 P) from", nperm, "permutations within site\n")
for(m in names(real)){
  r <- as.matrix(mx[mx$map == m, paste0("T", tps)])
  scan <- sapply(levels, function(q) apply(r, 2, quantile, q)); rownames(scan) <- tps
  act <- sapply(levels, function(q) quantile(apply(r[, tps <= 770, drop = FALSE], 1, max), q, names = FALSE))
  cy <- sapply(c("cauchy.T770", "cauchy.T1100"), function(k) quantile(mx[mx$map == m, k], 0.95, names = FALSE))
  th[[m]] <- list(scan = scan, act = act)
  cat(sprintf("\n%s: 5%% one T-age %.2f (T20) to %.2f (T770) and %.2f (T1100); any T-age to T770 %.2f; Cauchy 5%% %.2f (T20-T770), %.2f (T20-T1100)\n",
              m, scan[1, 1], scan[as.character(770), 1], scan[nrow(scan), 1], act[1], cy[1], cy[2]))
  L <- real[[m]]; cl <- cauchy(L[, tps <= 770])
  sgn <- if(m == "maternal") c("B", "C") else c("D", "H")
  for(ch in c(1:19, "X")){
    k <- which(umap$Chr == ch); if(length(k) == 0) next
    a <- abs(L[k, , drop = FALSE]); best <- apply(a, 2, max)
    lev <- ifelse(best >= scan[, 1], "significant", ifelse(best >= scan[, 2], "P 0.10", ifelse(best >= scan[, 3], "suggestive", "")))
    if(all(lev == "")) next
    w <- which(a == max(a), arr.ind = TRUE)[1, ]
    top <- if(any(lev == "significant")) "significant" else if(any(lev == "P 0.10")) "P 0.10" else "suggestive"
    cat(sprintf("  Chr %-2s %6.1f Mb %s: peak %.2f at T%d; best level %s at %s; Cauchy (T20-T770) at the chromosome's best marker %.2f\n", ch, umap$Mb[k[w[1]]],
                if(L[k[w[1]], w[2]] > 0) sgn[1] else sgn[2], a[w[1], w[2]], tps[w[2]], top, paste0("T", tps[lev == top], collapse = " "), max(cl[k])))
  }
}
sink()

### Figures, in the BXD maps' geometry (actuary_map.R; the chromosome layout of the BXD map)
bmap <- map; bmap[, "Mb"] <- 1e6 * bmap[, "Mb"]; lay <- mapLayout(bmap)
umapb <- data.frame(Chr = umap$Chr, Mb = umap$Mb * 1e6, row.names = markers)
gw <- function(m, T){ s <- th[[m]]$scan; sprintf("Genome-wide 5%% threshold, −log10 *P* at one T-age: %.1f at T20 to %.1f at T%d (any T-age to T770 %.1f). From %s permutations of lifespans within site",
                                                  s[1, 1], s[as.character(T), 1], T, th[[m]]$act[1], format(nperm, big.mark = ",")) }
note <- function(m) c(paste0("Model: lm(age at death ~ site + cohort + treatment + ", if(m == "maternal") "maternal B dosage, P(BH) + P(BD)" else "paternal D dosage, P(CD) + P(BD)",
                             "), likelihood-ratio test (1 d.f.), as in Danny Arends's ProgressiveMapping.R. UM-HET3 ", Sex, " alive at T; 891 markers (Arends et al. 2026, Supplementary Table 7)."),
                      paste("Mapping: umhet3_maps.R (robwwilliams/BxDLifeSpan).", if(m == "paternal") paste("Paternal calls on proximal Chr 1 are biased toward H at UM and UT.", if(sex == "female") "A female's paternal X is always H (her father's mother was C3H/HeJ)" else "A male has no paternal X", "so the paternal X is not mapped.") else ""))
panel <- function(m, tt, L = NULL, Nn = NULL, how = "", gwn = NULL){
  if(is.null(L)){ L <- real[[m]]; colnames(L) <- paste0("X", tps); L <- L[, paste0("X", tt), drop = FALSE]; Nn <- N[match(tt, tps)] }
  list(lods = L, N = Nn, N.title = "Mice", tps = tt, map = umapb, layout = lay,
       main = paste0("Lifespan of UM-HET3 ", Sex, ": ", if(m == "maternal") "maternal allele, B against C" else "paternal allele, D against H", how),
       sign.labels = if(m == "maternal") c("*B*: longer life", "*C*: longer life") else c("*D*: longer life", "*H*: longer life"),
       colz = if(m == "maternal") list(map.colz.n.wide, map.colz.g.wide) else list(map.colz.c.wide, map.colz.g.wide),   # the wider scale: UM-HET3 peaks reach 9
       thresholds = map.thresholds.wide, bins = map.bins.wide, gw.note = if(is.null(gwn)) gw(m, max(tt)) else gwn, model = note(m))
}
page <- function(file, tt, labs, extra = NULL){
  svglite(file, width = 24, height = 5.5 * 2 + 0.6, bg = "white", system_fonts = list(sans = "Arial"), fix_text_size = FALSE)
  layout(matrix(1:2), heights = c(5.5, 6.1))                                          # the bottom panel's extra 0.6 inch is all margin, for its notes
  par(omi = c(0, 0, 0, 0))                                                            # once: setting omi again would start a new page
  for(i in 1:2){
    par(mai = map.mai + c(if(i == 2) 0.6 else 0, 0, 0, 0), cex = 1)
    m <- c("maternal", "paternal")[i]
    p <- if(is.null(extra)) panel(m, tt) else panel(m, tt, extra$L[[m]], extra$N, extra$how, extra$gw[[m]]); p$letter <- letters[i]; p$chr.numbers <- i == 2; p$labels <- labs; do.call(actuaryMap, p)
  }
  invisible(dev.off())
}
page(paste0("output/maps/UMHET3_", sex, "_T50-T770.svg"), seq(50, 770, 30), labelTPS)               # the BXD maps' rows
page(paste0("output/maps/UMHET3_", sex, "_T50-T1100.svg"), seq(50, 1100, 30), seq(50, 1100, 60))    # the whole UM-HET3 lifespan

### Rows matched one to one with the BXD maps' rows (T50 to T770 in 30-day steps), so that each panel can be laid over a
### BXD panel. Rob, 24 November 2026: "perhaps some kind of percentage normalization would be best so that the plots can be
### laid over each other". Two matchings, against BXD females on chow (UM-HET3 mice eat chow):
###   grid40: T45 to T1005 in 40-day steps, one row per BXD row (T50 to T770 in 30-day steps). 40/30 = 1.33 is close to
###           the ratio of median ages at death (UM-HET3 females 900 days, BXD females on chow 680, 1.32), and the grid
###           lands on T365 and T725, as Arends et al.'s 15-day grid does (Rob: "I like your 40-day step size idea");
###   alive:  the UM-HET3 T-age at which the same fraction of mice is still alive as in BXD at that T-age.
### The left axis gives the UM-HET3 T-ages; the rows sit where the BXD maps' rows sit.
bxdT <- seq(50, 770, 30)
bxdAges <- ldata[ldata$Diet == "CD" & grepl("BXD", ldata$StrainName), "AgeAtDeath..days."]; bxdAges <- bxdAges[!is.na(bxdAges)]
matched <- list(grid40 = seq(45, 1005, 40),
                alive = round(sapply(bxdT, function(T) quantile(cdata$age, 1 - mean(bxdAges >= T), names = FALSE))))
how <- c(grid40 = " (rows: 40-day steps, one per BXD row)",
         alive = " (rows: the same fraction alive as BXD females on chow at each BXD row)")
for(k in names(matched)){
  Tm <- matched[[k]]; s <- lapply(Tm, function(T) scanT(cdata$age, T))
  L <- lapply(c(maternal = "maternal", paternal = "paternal"), function(m){ x <- sapply(s, function(z) z[, m]); colnames(x) <- paste0("X", Tm); x })
  gwk <- lapply(c(maternal = "maternal", paternal = "paternal"), function(m){ q <- approx(tps, th[[m]]$scan[, 1], xout = Tm, rule = 2)$y
    sprintf("Genome-wide 5%% threshold, −log10 *P* at one T-age: %.1f to %.1f over these rows. From %s permutations of lifespans within site", min(q), max(q), format(nperm, big.mark = ",")) })
  page(paste0("output/maps/UMHET3_", sex, "_rows_matched_", k, ".svg"), Tm, Tm[seq(1, length(Tm), 2)],
       list(L = L, N = sapply(Tm, function(T) sum(cdata$age >= T)), how = how[k], gw = gwk))
  cat("Rows matched by", k, ":", paste0("T", Tm, collapse = " "), "\n")
}
cat("Wrote the maps, thresholds and figures\n")
