#
# qtl2_individuals.R - Actuarial mapping of individual mice with R/qtl2: linear mixed models with LOCO kinship
#
# Rob Williams's fork, 21 November 2026. Exploratory. Rob: "even better is to simply map the individuals with a full
# kinship correction." The mice, T-ages, genotypes and cofactors are those of Danny's interactions.R; the kinship matrix
# takes the place of his random strain intercept.
#   - The cross: the BXD strains of shared.R's ldata as a sib-mated RI panel (risib), genotypes from
#     data/BXD_current_rev050423.geno, H and U as missing. The X is treated like an autosome, as in Danny's maps.
#   - Genotype probabilities and LOCO kinship (allele sharing) come from the strains. Each mouse takes its strain's, so
#     mice of one strain have kinship 1 with each other.
#   - At each T-age (20 to 770 days, 30-day steps), BXD mice alive at T:
#       main  both diets: age at death ~ diet + vivarium + set-up class + genotype, with kinship (scan1, REML)
#       int   both diets: the same + diet x genotype; -log10 P from the difference of the two LODs (1 d.f.)
#       CD    chow alone:      age at death ~ vivarium + set-up class + genotype, with kinship
#       HF    high fat alone:  the same
#     Set-up class: the mouse's set-up age in 180-day classes, as in Danny's code. Mice without a set-up age are left out.
#   - Signed -log10 P: chi-square on 1 d.f. of 2 ln(10) LOD, with the sign of the D - B effect from scan1coef, positive
#     when the D allele raises lifespan (for int: when D - B is larger on high fat), as in Danny's maps.
# Writes output/qtl2_<CD|HF|main|int>_ind_Progressive.txt (markers x T-ages, Danny's layout), output/qtl2_N_ind_Progressive.txt
# (mice per model and T-age) and output/qtl2_hsq.txt (the kinship heritability of each null model, chromosome 1 left out)
#   Rscript qtl2_individuals.R        (about 20 minutes on 8 cores)
#

source("shared.R")
source("qtl2_cross.R")                                                              # strains, cross, pr, K and signedP

tps <- seq(20, 770, 30)
cores <- qtl2.cores

### Each mouse takes its strain's genotype probabilities and kinship
ids <- rownames(ld)
idx <- match(ld[, "StrainName"], strains)
prm <- pr
for(chr in names(pr)){ a <- pr[[chr]][idx, , , drop = FALSE]; dimnames(a)[[1]] <- ids; prm[[chr]] <- a }
Km <- lapply(K, function(k){ m <- k[idx, idx]; dimnames(m) <- list(ids, ids); m })

# Effect of D against B at every marker (and of D against B times high fat, when int is given)
effects <- function(Y, Kd, X, int = NULL){
  unlist(lapply(names(prm), function(chr){
    co <- scan1coef(prm[, chr], Y, kinship = Kd[[chr]], addcovar = X, intcovar = int)
    if(is.null(int)) return(co[, "DD"] - co[, "BB"])
    ci <- colnames(co)[grepl("dietHF", colnames(co)) & grepl("BB|DD", colnames(co))]
    d <- if("DD:dietHF" %in% ci) co[, "DD:dietHF"] else 0
    b <- if("BB:dietHF" %in% ci) co[, "BB:dietHF"] else 0
    return(d - b)
  }))
}

res <- list(CD = NULL, HF = NULL, main = NULL, int = NULL)
nInd <- NULL
hsq <- NULL
for(T in tps){
  x <- ld[which(ld[, "AgeAtDeath..days."] >= T),]
  cv <- data.frame(diet = x[, "Diet"], vivarium = as.factor(x[, "vivarium"]), row.names = rownames(x),
                   setup = as.factor(cut(as.numeric(x[, "AgeAtSetUp.in.colony..days."]), seq(0, 1000, 180))))
  cv <- cv[complete.cases(cv), ]
  ns <- c()
  hs <- c()
  for(m in c("main", "CD", "HF")){
    cm <- if(m == "main") cv else cv[cv[, "diet"] == m, ]
    X <- model.matrix(if(m == "main") ~ diet + vivarium + setup else ~ vivarium + setup, cm)[, -1, drop = FALSE]
    X <- X[, colSums(X) > 0, drop = FALSE]
    Y <- matrix(x[rownames(cm), "AgeAtDeath..days."], dimnames = list(rownames(cm), "age"))
    Kd <- decomp_kinship(lapply(Km, function(k) k[rownames(cm), rownames(cm)]), cores = cores)
    out <- scan1(prm, Y, kinship = Kd, addcovar = X, cores = cores)
    res[[m]] <- cbind(res[[m]], signedP(out[, 1], effects(Y, Kd, X)))
    ns <- c(ns, nrow(cm)); hs <- c(hs, attr(out, "hsq")[1, 1])
    if(m == "main"){
      int <- X[, "dietHF", drop = FALSE]
      out.int <- scan1(prm, Y, kinship = Kd, addcovar = X, intcovar = int, cores = cores)
      res[["int"]] <- cbind(res[["int"]], signedP(out.int[, 1] - out[, 1], effects(Y, Kd, X, int)))
    }
  }
  nInd <- rbind(nInd, ns)
  hsq <- rbind(hsq, hs)
  cat("T", T, ": mice", paste(c("both", "CD", "HF"), ns, collapse = ", "), "; hsq", paste(round(hs, 2), collapse = " "), "\n")
}

for(m in names(res)){
  out <- cbind(map, res[[m]])
  colnames(out) <- c("Chr", "Locus", "cM", "Mb", tps)
  write.table(out, file = paste0("output/qtl2_", m, "_ind_Progressive.txt"), sep = "\t", quote = FALSE)
}
dimnames(nInd) <- list(tps, c("main", "CD", "HF"))
dimnames(hsq) <- list(tps, c("main", "CD", "HF"))
write.table(nInd, file = "output/qtl2_N_ind_Progressive.txt", sep = "\t", quote = FALSE)
write.table(round(hsq, 4), file = "output/qtl2_hsq.txt", sep = "\t", quote = FALSE)
