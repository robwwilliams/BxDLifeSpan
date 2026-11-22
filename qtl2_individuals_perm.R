#
# qtl2_individuals_perm.R - Genome-wide thresholds for the individual-mouse maps (qtl2_individuals.R), by permuting
# strains within epochs
#
# Rob Williams's fork, 21 November 2026. Exploratory. Rob: "what are the exchangable units? Do you need a hierarchical
# permutation? Important not to over-randomize across the intrinsic structure of the cross or the P values will be way too
# lenient." And: "200 permutations is good enough to get a sense of the 0.1 genome wide level and a crude sense of 0.05."
#   - The exchangeable unit is the strain, not the mouse. In each permutation every strain takes the genotypes (and so the
#     kinship) of another strain, and keeps all its mice, on both diets, at every T-age, with their covariates (diet,
#     vivarium, set-up class). The survivorships and their nesting are kept; only the link between genotype and lifespan is
#     broken.
#   - Strains are permuted only within two strata: Taylor's strains (BXD1 to BXD42, epochs 1 and 2) and the newer ones
#     (BXD43 on, epochs 3 to 6), which differ in their relatedness. Epochs are from GeneNetwork's BDL_10002 export; strains
#     without one there are placed by number.
#   - Each permutation reruns the four maps at all 26 T-ages (chow, high fat, both diets, and genotype x diet), as in
#     qtl2_individuals.R but without the effect signs, and records the largest -log10 P over markers at each T-age.
#   - The interaction test under this permutation has a stricter null (no genetic effect at all), so its thresholds lean
#     conservative.
#   - Permutations are drawn in advance from one seed; results are appended to the output as they finish, so the run can
#     be stopped and resumed.
# Writes output/qtl2_individuals_perm_maxima.txt (one row per permutation and map)
#   Rscript qtl2_individuals_perm.R [permutations, default 200] [parallel workers, default 8] [BDL_10002.xlsx path]
#

source("shared.R")
source("qtl2_cross.R")                                                              # strains, cross, pr, K
library(parallel)
library(readxl)

args <- commandArgs(trailingOnly = TRUE)
nperm <- if(length(args) >= 1) as.numeric(args[1]) else 200
workers <- if(length(args) >= 2) as.numeric(args[2]) else 8
bdl <- if(length(args) >= 3) args[3] else "../bxd_mortality_rates_2026/data/original/BDL_10002.xlsx"
tps <- seq(20, 770, 30)
out.file <- "output/qtl2_individuals_perm_maxima.txt"
lp <- function(lod) -pchisq(2 * log(10) * pmax(lod, 0), 1, lower.tail = FALSE, log.p = TRUE) / log(10)

# Strata: Taylor's strains (epochs 1 and 2) and the newer ones
x <- as.data.frame(suppressMessages(suppressWarnings(read_excel(bdl, skip = 10))))
x <- x[!is.na(x[, "Epoch"]), ]
epoch <- tapply(x[, "Epoch"], x[, "Strain"], function(v) names(sort(table(v), decreasing = TRUE))[1])[strains]
num <- as.numeric(sub("^BXD([0-9]+).*", "\\1", strains))
stratum <- ifelse(is.na(epoch), ifelse(num <= 42, "Taylor", "newer"), ifelse(substr(epoch, 1, 1) %in% c("1", "2"), "Taylor", "newer"))
names(stratum) <- strains
cat("Strata:", paste(names(table(stratum)), table(stratum), collapse = ", "), "\n")

set.seed(20261123)
perms <- lapply(1:nperm, function(p){ o <- strains; for(s in unique(stratum)){ i <- which(stratum == s); o[i] <- strains[i][sample.int(length(i))] }; names(o) <- strains; o })

ids <- rownames(ld)
cv.all <- data.frame(diet = ld[, "Diet"], vivarium = as.factor(ld[, "vivarium"]), row.names = ids,
                     setup = as.factor(cut(as.numeric(ld[, "AgeAtSetUp.in.colony..days."]), seq(0, 1000, 180))))

one <- function(p){
  gs <- perms[[p]][ld[, "StrainName"]]                                               # each mouse takes the genotypes of its strain's partner
  idx <- match(gs, strains)
  prm <- pr
  for(chr in names(pr)){ a <- pr[[chr]][idx, , , drop = FALSE]; dimnames(a)[[1]] <- ids; prm[[chr]] <- a }
  Km <- lapply(K, function(k){ m <- k[idx, idx]; dimnames(m) <- list(ids, ids); m })
  res <- matrix(NA, 4, length(tps), dimnames = list(c("main", "int", "CD", "HF"), paste0("T", tps)))
  for(j in seq_along(tps)){
    keep <- ids[ld[, "AgeAtDeath..days."] >= tps[j]]
    cv <- cv.all[keep, ]; cv <- cv[complete.cases(cv), ]
    for(m in c("main", "CD", "HF")){
      cm <- if(m == "main") cv else cv[cv[, "diet"] == m, ]
      X <- model.matrix(if(m == "main") ~ diet + vivarium + setup else ~ vivarium + setup, cm)[, -1, drop = FALSE]
      X <- X[, colSums(X) > 0, drop = FALSE]
      Y <- matrix(ld[rownames(cm), "AgeAtDeath..days."], dimnames = list(rownames(cm), "age"))
      Kd <- decomp_kinship(lapply(Km, function(k) k[rownames(cm), rownames(cm)]))
      out <- scan1(prm, Y, kinship = Kd, addcovar = X)
      res[m, j] <- max(lp(out[, 1]))
      if(m == "main"){
        out.int <- scan1(prm, Y, kinship = Kd, addcovar = X, intcovar = X[, "dietHF", drop = FALSE])
        res["int", j] <- max(lp(out.int[, 1] - out[, 1]))
      }
    }
  }
  r <- data.frame(perm = p, map = rownames(res), res, check.names = FALSE)
  write.table(r, out.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = !file.exists(out.file), append = TRUE)
  invisible(r)
}

done <- if(file.exists(out.file)) unique(read.table(out.file, sep = "\t", header = TRUE)[, "perm"]) else c()
todo <- setdiff(1:nperm, done)
cat(length(done), "permutations already done;", length(todo), "to run on", workers, "workers\n")
t0 <- Sys.time()
invisible(mclapply(todo, one, mc.cores = workers, mc.preschedule = FALSE))
cat("Done after", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "minutes\n")
