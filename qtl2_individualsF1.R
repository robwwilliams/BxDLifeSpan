#
# qtl2_individualsF1.R - Actuarial mapping of individual mice, BXD and BXD x BXD F1 together, with R/qtl2 and LOCO
# kinship: chow, high fat, both diets (main effect), and genotype x diet
#
# Rob Williams's fork, 23 November 2026. Exploratory. For Rob's Fig. 2 (23 November 2026): "b. individuals with F1s
# stratified by CD and HFD; c. individuals with F1s, all combined, all cofactors, main effect only; d. individuals all
# combined, diet x locus interaction effect only (QUESTION if F1s should be included?)". The int map is made with the F1s
# here; qtl2_individuals.R has it without them.
#   - Mice: the BXD mice of qtl2_individuals.R (shared.R's ldata, CD and HF, with a set-up age), and the F1 mice of
#     strainF1_setup.R (Suheeta Roy's F1 workbook, natural deaths on CD and HF). All are female.
#   - Genotypes: each mouse takes its unit's allele probabilities (strain, or F1 type from both parents, so 1/2 D where
#     the parents differ) and its unit's LOCO kinship over strains and F1 types (strainF1_setup.R). Additive (1 d.f.).
#   - Covariates, as in qtl2_individuals.R, plus F1 class:
#       CD, HF  age at death ~ vivarium + set-up class + F1 class + genotype
#       main    age at death ~ diet + vivarium + set-up class + F1 class + F1 class x diet + genotype
#       int     main + genotype x diet; -log10 P from the difference of the two LODs (1 d.f.)
#     F1 class is two indicators, short x short and long x long, so the F1s' mean departure from the strains (hybrid
#     vigor, and anything else they share) is not read as a genotype effect. The F1 workbook has no vivarium column;
#     the F1s (born 2020-2022) were all housed at TSRB (Rob, 23 November 2026: "I believe all F1s were born in TSRB"),
#     vivarium 3. With the F1 class indicators in the model, any vivarium coding of the F1s gives the same fit. In
#     main and int, F1 class x diet does the same for the difference in hybrid vigor between diets (on chow short x short
#     F1s gain about 196 days over their parents, on high fat about 146), which would otherwise leak into genotype x diet
#     wherever the F1s are heterozygous.
#   - Set-up class: set-up age in 180-day classes, as in Danny's code. F1s were set up at 20-44 days.
#   - Signed -log10 P: chi-square on 1 d.f., with the sign of the D - B effect from scan1coef (for int: of D - B being
#     larger on high fat).
# Writes output/qtl2_<CD|HF|main|int>_indF1_Progressive.txt, output/qtl2_N_indF1_Progressive.txt (mice, and F1 mice, per
# model and T-age) and output/qtl2_hsq_indF1.txt
#   Rscript qtl2_individualsF1.R        (about 30 minutes on 8 cores)
#

source("shared.R")
source("qtl2_cross.R")                                                              # strains, pr, signedP
source("strainF1_setup.R")                                                         # f1, types, f1class, AP, Kloco
library(parallel)

tps <- seq(20, 770, 30)
cores <- qtl2.cores

### One table of mice: BXD and F1
bxd <- data.frame(id = rownames(ld), unit = ld[, "StrainName"], diet = ld[, "Diet"], age = ld[, "AgeAtDeath..days."],
                  vivarium = as.character(ld[, "vivarium"]), setup.age = as.numeric(ld[, "AgeAtSetUp.in.colony..days."]))
fm <- data.frame(id = paste0("F1_", f1[, "EarTagNumberCurrent"]), unit = f1[, "StrainName"], diet = f1[, "DietCode"], age = f1[, "AgeAtDeath (days)"],
                 vivarium = "3", setup.age = as.numeric(f1[, "AgeAtSet Up (days)"]))                     # TSRB (Rob, 23 November 2026)
mice <- rbind(bxd, fm)
mice <- mice[!is.na(mice$age) & !is.na(mice$setup.age) & !is.na(mice$vivarium) & mice$unit %in% units, ]
stopifnot(!any(duplicated(mice$id)))
rownames(mice) <- mice$id
mice$class <- ifelse(mice$unit %in% types, f1class[mice$unit], "inbred")
mice$setup <- cut(mice$setup.age, seq(0, 1000, 180))
cat(sum(mice$class == "inbred"), "BXD mice and", sum(mice$class != "inbred"), "F1 mice\n")

### Each mouse takes its unit's allele probabilities and kinship
ids <- mice$id; idx <- match(mice$unit, units)
apm <- AP
for(chr in names(AP)){ a <- AP[[chr]][idx, , , drop = FALSE]; dimnames(a)[[1]] <- ids; apm[[chr]] <- a }
Km <- lapply(Kloco, function(k){ m <- k[idx, idx]; dimnames(m) <- list(ids, ids); m })

effects <- function(Y, Kd, X, int = NULL){
  unlist(lapply(names(apm), function(chr){
    co <- scan1coef(apm[, chr], Y, kinship = Kd[[chr]], addcovar = X, intcovar = int)
    if(is.null(int)) return(co[, "D"] - co[, "B"])
    d <- if("D:dietHF" %in% colnames(co)) co[, "D:dietHF"] else 0
    b <- if("B:dietHF" %in% colnames(co)) co[, "B:dietHF"] else 0
    d - b
  }))
}
covars <- function(cm, m){
  f <- if(m == "main") ~ diet + vivarium + setup + short + long + short:diet + long:diet else ~ vivarium + setup + short + long
  cm$short <- as.numeric(cm$class == "short"); cm$long <- as.numeric(cm$class == "long")
  cm$vivarium <- factor(cm$vivarium); cm$setup <- droplevels(cm$setup)
  X <- model.matrix(f, cm)[, -1, drop = FALSE]
  X <- X[, colSums(X) > 0, drop = FALSE]
  q <- qr(X); X[, q$pivot[seq_len(q$rank)], drop = FALSE]                           # no collinear columns
}

res <- list(CD = NULL, HF = NULL, main = NULL, int = NULL)
nInd <- NULL; hsq <- NULL
for(T in tps){
  cv <- mice[mice$age >= T, ]
  ns <- c(); hs <- c()
  for(m in c("main", "CD", "HF")){
    cm <- if(m == "main") cv else cv[cv$diet == m, ]
    X <- covars(cm, m)
    Y <- matrix(cm$age, dimnames = list(cm$id, "age"))
    Kd <- decomp_kinship(lapply(Km, function(k) k[cm$id, cm$id]), cores = cores)
    out <- scan1(apm, Y, kinship = Kd, addcovar = X, cores = cores)
    res[[m]] <- cbind(res[[m]], signedP(out[, 1], effects(Y, Kd, X)))
    ns <- c(ns, nrow(cm), sum(cm$class != "inbred")); hs <- c(hs, attr(out, "hsq")[1, 1])
    if(m == "main"){
      int <- X[, "dietHF", drop = FALSE]
      out.int <- scan1(apm, Y, kinship = Kd, addcovar = X, intcovar = int, cores = cores)
      res[["int"]] <- cbind(res[["int"]], signedP(out.int[, 1] - out[, 1], effects(Y, Kd, X, int)))
    }
  }
  nInd <- rbind(nInd, ns); hsq <- rbind(hsq, hs)
  cat("T", T, ": mice (F1)", paste(c("both", "CD", "HF"), sprintf("%d (%d)", ns[c(1, 3, 5)], ns[c(2, 4, 6)]), collapse = ", "), "; hsq", paste(round(hs, 2), collapse = " "), "\n")
}
for(m in names(res)){
  out <- cbind(map, res[[m]]); colnames(out) <- c("Chr", "Locus", "cM", "Mb", tps)
  write.table(out, file = paste0("output/qtl2_", m, "_indF1_Progressive.txt"), sep = "\t", quote = FALSE)
}
dimnames(nInd) <- list(tps, c("main", "main.F1", "CD", "CD.F1", "HF", "HF.F1"))
dimnames(hsq) <- list(tps, c("main", "CD", "HF"))
write.table(nInd, file = "output/qtl2_N_indF1_Progressive.txt", sep = "\t", quote = FALSE)
write.table(round(hsq, 4), file = "output/qtl2_hsq_indF1.txt", sep = "\t", quote = FALSE)
