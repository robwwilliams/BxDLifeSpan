#
# strainF1_functions.R - the steps of the maps of BXD strain and F1 means together (qtl2_strainF1.R) and of their permutations
# (qtl2_strainF1_perm.R), one function per step. Sourcing this file defines the functions and runs nothing else.
#
# Rob Williams's fork, 26 November 2026. Written by Claude (Opus 5.5) at Danny Arends's request. Checked by strainF1_tests.R.
#
# THE MODEL, fitted at every marker, separately for each diet and each T-age, by R/qtl2 scan1 (REML, LOCO kinship):
#
#     unit mean lifespan  ~  F1 class  +  genotype at the marker  +  polygenic background (kinship)
#
#   unit                 one BXD strain or one BXD x BXD F1 type, with at least 3 mice alive at the T-age
#   unit mean lifespan   the mean age at death (days) of the unit's mice alive at T. Each unit counts once, however many mice it has.
#                        Not adjusted for vivarium or set-up age (Rob's "simple method").
#   F1 class             two 0/1 columns: short x short F1, long x long F1. Inbred strains are the baseline. They absorb the average
#                        hybrid vigor of each class, so that it is not read as an effect of every marker.
#   genotype             additive scan: the D allele dosage, 0 (a BB strain), 0.5 (an F1 that is BD) or 1 (a DD strain); 1 d.f.
#                        full scan: BB, BD and DD as three classes; 2 d.f.
#                        dominance: full LOD minus additive LOD; 1 d.f. Only the F1s inform it, as only they can be BD.
#   polygenic background a random effect whose covariance is proportional to kinship: the alleles two units share over all markers
#                        except those of the chromosome scanned (LOCO, leave one chromosome out)
#
# Each LOD becomes a signed -log10 P on 1 d.f.: additive, positive when the D allele lives longer; dominance, positive when BD lives
# longer than the mean of BB and DD.
#

library(qtl2)

# ==== 1. The F1 mice and F1 types =====================================================================================================
# The F1 workbook: natural deaths on chow (CD) or high fat (HF). Returns a data frame, one row per F1 mouse.
readF1Mice <- function(f1.file){
  f1.mice <- read.delim(f1.file, check.names = FALSE)
  f1.mice[f1.mice[, "CauseOfDeath"] == "Natural" & f1.mice[, "DietCode"] %in% c("CD", "HF"), ]
}

# The F1 types and their parents. A type is named <parent 1>x<parent 2 without "BXD">, as "BXD102x73b" (BXD102 x BXD73b).
# Returns a data frame, one row per F1 type (row names: the type): parent1, parent2, class ("short" or "long").
f1TypeTable <- function(f1.mice, strains, short.parents, long.parents){
  type.names <- sort(unique(f1.mice[, "StrainName"]))
  parent1 <- sub("x.*", "", type.names); parent2 <- paste0("BXD", sub(".*x", "", type.names))
  stopifnot(all(c(parent1, parent2) %in% strains))
  f1.class <- ifelse(parent1 %in% short.parents & parent2 %in% short.parents, "short",
              ifelse(parent1 %in% long.parents & parent2 %in% long.parents, "long", NA))
  stopifnot(!any(is.na(f1.class)))                                          # every F1 type is short x short or long x long
  data.frame(parent1 = parent1, parent2 = parent2, class = f1.class, row.names = type.names)
}

# ==== 2. Genotype probabilities of every unit =========================================================================================
# strain.probs: the strains' genotype probabilities from qtl2_cross.R (calc_genoprob; per chromosome, strains x BB/DD x markers).
# An F1 takes one allele from each inbred parent, so at each marker:
#   P(BB) = P1(B) P2(B);   P(BD) = P1(B) P2(D) + P1(D) P2(B);   P(DD) = P1(D) P2(D)
# Returns a list of two genotype-probability objects, rows = units (strains, then F1 types):
#   three.genotypes  per chromosome, units x (BB, BD, DD) x markers, for the full scan
#   allele.dosage    per chromosome, units x (B, D) x markers: the share of B and of D alleles, for the additive scan
unitGenotypes <- function(strain.probs, strains, f1.types){
  units <- c(strains, rownames(f1.types)); p1 <- f1.types$parent1; p2 <- f1.types$parent2
  three.genotypes <- strain.probs; allele.dosage <- strain.probs
  for(chr in names(strain.probs)){
    prob.B <- strain.probs[[chr]][, "BB", ]; prob.D <- strain.probs[[chr]][, "DD", ]          # strains x markers
    marker.names <- dimnames(strain.probs[[chr]])[[3]]
    g <- array(0, c(length(units), 3, length(marker.names)), list(units, c("BB", "BD", "DD"), marker.names))
    g[strains, "BB", ] <- prob.B[strains, ]; g[strains, "DD", ] <- prob.D[strains, ]           # a strain is BB or DD
    g[rownames(f1.types), "BB", ] <- prob.B[p1, ] * prob.B[p2, ]
    g[rownames(f1.types), "BD", ] <- prob.B[p1, ] * prob.D[p2, ] + prob.D[p1, ] * prob.B[p2, ]
    g[rownames(f1.types), "DD", ] <- prob.D[p1, ] * prob.D[p2, ]
    three.genotypes[[chr]] <- g
    a <- array(0, c(length(units), 2, length(marker.names)), list(units, c("B", "D"), marker.names))
    a[, "B", ] <- g[, "BB", ] + g[, "BD", ] / 2; a[, "D", ] <- g[, "DD", ] + g[, "BD", ] / 2
    allele.dosage[[chr]] <- a
  }
  attr(three.genotypes, "crosstype") <- "f2"; attr(three.genotypes, "alleleprobs") <- FALSE
  attr(allele.dosage, "alleleprobs") <- TRUE
  list(three.genotypes = three.genotypes, allele.dosage = allele.dosage)
}

# ==== 3. Kinship of every pair of units, leaving out one chromosome at a time (LOCO) ==================================================
# Kinship of two units = alleles they share, summed over the markers of all other chromosomes, divided by the number of those markers.
# Returns a list, one units x units matrix per chromosome scanned.
locoKinship <- function(allele.dosage){
  shared.per.chr <- lapply(allele.dosage, function(a){ m <- matrix(a, nrow = dim(a)[1]); tcrossprod(m) })   # units x units, summed
  markers.per.chr <- sapply(allele.dosage, function(a) dim(a)[3])
  units <- dimnames(allele.dosage[[1]])[[1]]
  kinship <- lapply(names(allele.dosage), function(chr){
    k <- (Reduce(`+`, shared.per.chr) - shared.per.chr[[chr]]) / (sum(markers.per.chr) - markers.per.chr[chr])
    dimnames(k) <- list(units, units); k })
  setNames(kinship, names(allele.dosage))
}

# ==== 4. The trait at one diet and T-age ==============================================================================================
# strain.mice: BXD mice of qtl2_cross.R (ld); f1.mice: from readF1Mice(). Returns a named numeric vector: unit -> mean age at death of
# its mice alive at T, for the units with at least min.mice such mice (strains first, then F1 types).
unitMeansAt <- function(strain.mice, f1.mice, diet, Tage, min.mice){
  alive.strain <- strain.mice[which(strain.mice[, "Diet"] == diet & strain.mice[, "AgeAtDeath..days."] >= Tage), ]
  alive.f1 <- f1.mice[f1.mice[, "DietCode"] == diet & f1.mice[, "AgeAtDeath (days)"] >= Tage, ]
  mice.per.strain <- table(alive.strain[, "StrainName"])                       # strain -> mice alive at T
  mice.per.f1type <- table(alive.f1[, "StrainName"])                           # F1 type -> mice alive at T
  strains.kept <- names(mice.per.strain)[mice.per.strain >= min.mice]
  f1types.kept <- names(mice.per.f1type)[mice.per.f1type >= min.mice]
  c(tapply(alive.strain[, "AgeAtDeath..days."], alive.strain[, "StrainName"], mean)[strains.kept],
    tapply(alive.f1[, "AgeAtDeath (days)"], alive.f1[, "StrainName"], mean)[f1types.kept])
}

# The F1 class covariate for a set of units. Returns a matrix, units x (short, long) of 0/1; a column no unit has is dropped, and NULL
# is returned when no F1 type is among the units.
classCovariate <- function(unit.names, f1.types){
  unit.class <- ifelse(unit.names %in% rownames(f1.types), f1.types[unit.names, "class"], "inbred")
  covariate <- cbind(short = as.numeric(unit.class == "short"), long = as.numeric(unit.class == "long")); rownames(covariate) <- unit.names
  covariate <- covariate[, colSums(covariate) > 0, drop = FALSE]
  if(ncol(covariate) == 0) NULL else covariate
}

# ==== 5. The two scans at one diet and T-age ==========================================================================================
# unit.means: from unitMeansAt(). genotypes: from unitGenotypes(). kinship: from locoKinship(). f1.types: from f1TypeTable().
# Returns a list: additive and dominance (signed -log10 P per marker), additive.lod, and the number of strains and F1 types used.
scanStrainsAndF1s <- function(unit.means, genotypes, kinship, f1.types){
  units <- names(unit.means)
  lifespan <- matrix(unit.means, dimnames = list(units, "lifespan"))           # units x 1, the response
  covariate <- classCovariate(units, f1.types)
  kinship.units <- decomp_kinship(lapply(kinship, function(k) k[units, units]))
  dosage <- genotypes$allele.dosage[units, ]; three <- genotypes$three.genotypes[units, ]
  lod.additive <- scan1(dosage, lifespan, kinship = kinship.units, addcovar = covariate)    # markers x 1
  lod.full <- scan1(three, lifespan, kinship = kinship.units, addcovar = covariate)
  effects <- do.call(rbind, lapply(names(dosage), function(chr){                # the direction of each effect, marker by marker
    coef.additive <- scan1coef(dosage[, chr], lifespan, kinship = kinship.units[[chr]], addcovar = covariate)
    coef.three <- scan1coef(three[, chr], lifespan, kinship = kinship.units[[chr]], addcovar = covariate)
    coef.BD <- if("BD" %in% colnames(coef.three)) coef.three[, "BD"] else rep(0, nrow(coef.three))
    cbind(additive = coef.additive[, "D"] - coef.additive[, "B"],                # days, D allele against B allele
          dominance = coef.BD - (coef.three[, "BB"] + coef.three[, "DD"]) / 2)  # days, BD against the mean of BB and DD
  }))
  list(additive = signedP(lod.additive[, 1], effects[, "additive"]),
       dominance = signedP(pmax(lod.full[, 1] - lod.additive[, 1], 0), effects[, "dominance"]),
       additive.lod = lod.additive[, 1], n = c(strains = sum(!units %in% rownames(f1.types)), F1 = sum(units %in% rownames(f1.types))))
}

# ==== 6. Permutations ================================================================================================================
# unit.means.all: units x T-ages of unit means (NA where a unit has fewer than min.mice mice alive at T).
unitMeansAllT <- function(strain.mice, f1.mice, diet, Tages, units, min.mice){
  means <- sapply(Tages, function(Tage){ m <- unitMeansAt(strain.mice, f1.mice, diet, Tage, min.mice); v <- setNames(rep(NA_real_, length(units)), units); v[names(m)] <- m; v })
  dimnames(means) <- list(units, paste0("T", Tages))
  means[rowSums(!is.na(means)) > 0, , drop = FALSE]                            # units with a mean at one T-age at least
}

# One permutation: each strain takes the whole row of means (every T-age) of another strain, and each F1 type that of another F1 type of
# its class. Units keep their genotypes and kinship. Returns the new order of the rows of unit.means.all. Uses the current random stream.
drawPermutation <- function(unit.names, f1.types){
  groups <- list(strains = which(!unit.names %in% rownames(f1.types)),
                 short = which(unit.names %in% rownames(f1.types)[f1.types$class == "short"]),
                 long = which(unit.names %in% rownames(f1.types)[f1.types$class == "long"]))
  new.order <- seq_along(unit.names)
  for(g in groups) new.order[g] <- g[sample.int(length(g))]
  new.order
}

# The largest -log10 P over markers at each T-age, for the additive and dominance maps of one (permuted) set of unit means.
# Returns a 2 x T-ages matrix (rows additive, dominance).
permutationMaxima <- function(unit.means.all, genotypes, kinship, f1.types){
  units <- rownames(unit.means.all)
  covariate <- classCovariate(units, f1.types); kinship.units <- lapply(kinship, function(k) k[units, units])
  lod.additive <- scan1(genotypes$allele.dosage[units, ], unit.means.all, kinship = kinship.units, addcovar = covariate)
  lod.full <- scan1(genotypes$three.genotypes[units, ], unit.means.all, kinship = kinship.units, addcovar = covariate)
  log10P <- function(lod) -pchisq(2 * log(10) * pmax(lod, 0), 1, lower.tail = FALSE, log.p = TRUE) / log(10)
  rbind(additive = apply(log10P(lod.additive), 2, max, na.rm = TRUE), dominance = apply(log10P(lod.full - lod.additive), 2, max, na.rm = TRUE))
}

# A run record: settings, input and code fingerprints (md5), R and qtl2 versions, time
writeRunInfo <- function(file, settings, inputs, code, minutes){
  writeLines(c(paste("Run:", format(Sys.time(), "%Y-%m-%d %H:%M %Z"), "on", Sys.info()[["nodename"]], "-", R.version.string, "- qtl2", as.character(packageVersion("qtl2"))),
               paste("Settings:", paste(names(settings), sapply(settings, paste, collapse = ","), sep = " = ", collapse = "; ")),
               paste("Inputs (md5):", paste(inputs, tools::md5sum(inputs), collapse = "; ")),
               paste("Code (md5):", paste(code, tools::md5sum(code), collapse = "; ")),
               paste("Minutes:", round(minutes, 1))), file)
}
