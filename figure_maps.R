#
# figure_maps.R - Danny's actuarial maps (Fig. 2) redrawn with one fixed geometry (actuary_map.R), so that every panel
# superimposes with pixel accuracy
#
# Rob Williams's fork, 21 November 2026. Reads Danny's result files in output/ (from lifespan.R, interactions.R,
# bodyweight.R and CTL.R; run_CTL.R makes the CTL ones; qtl2_individuals.R the kinship maps) when they are there. Otherwise reads the two workbooks Danny sent
# on 5 November 2026, which hold the same results: give their paths.
#   Rscript figure_maps.R [QTL_Lifespan_05Nov2026.xlsx QTL_BW_5Nov.xlsx]
# The numbers on the right axis are Danny's _N files when present; otherwise they are counted here with his rules: BXD
# strains with at least 3 mice alive at T (lifespan strain means), BXD mice on CD and HF alive at T (individuals), and BXD
# strains with at least 3 mice weighed at T (body weight, once bodyweight.R or run_CTL.R has interpolated the weights).
# Writes output/maps/: one SVG per panel (24 x 5.5 inches); Danny's pages of two; Figure 2 (lifespan), Figure 3 (body
# weight and CTL) and the R/qtl2 kinship maps of Figure 2; and every panel on one page, for alignment by eye
#

library(readxl)
source("shared.R")
source("actuary_map.R")

args <- commandArgs(trailingOnly = TRUE)
xlsx.lifespan <- if(length(args) >= 1) args[1] else "QTL_Lifespan_05Nov2026.xlsx"
xlsx.bw <- if(length(args) >= 2) args[2] else "QTL_BW_5Nov.xlsx"

map[, "Mb"] <- 1e6 * map[, "Mb"]
layout <- mapLayout(map)
dir.create("output/maps", showWarnings = FALSE, recursive = TRUE)

# The last commit that changed a file, for the note under each map
commitOf <- function(f){ h <- suppressWarnings(system(paste("git log -1 --format=%h --", f), intern = TRUE)); if(length(h) == 0) "uncommitted" else h }
drawn <- paste0("Drawn by actuary_map.R (robwwilliams/BxDLifeSpan ", commitOf("actuary_map.R"), ").")
mapping <- function(f, repo = "DannyArends/BxDLifeSpan"){ paste0("Mapping: ", f, " (", repo, " ", commitOf(f), ").") }

# A result as markers x T-ages, columns X20 ... X770: Danny's text file if it exists, else the workbook sheet
readResult <- function(txt, xlsx, sheet){
  if(file.exists(txt)){ res <- read.table(txt, sep = "\t"); attr(res, "source") <- txt; return(res) }
  if(!file.exists(xlsx)) return(NULL)
  x <- as.data.frame(suppressMessages(read_excel(xlsx, sheet = sheet, col_names = FALSE)))
  hdr <- unlist(x[1, ])
  res <- suppressWarnings(sapply(x[-1, -(1:5)], as.numeric))
  colnames(res) <- paste0("X", hdr[-(1:5)])
  res <- res[, which(colSums(!is.na(res)) > 0 & !duplicated(colnames(res))), drop = FALSE]   # the repeated last header is empty
  rownames(res) <- x[-1, 1]
  attr(res, "source") <- paste0(basename(xlsx), ", sheet ", sheet)
  return(res)
}
readN <- function(txt, counted){ if(file.exists(txt)) read.table(txt, sep = "\t")[match(plotTPS, seq(20, 780, 30)), 1] else counted }

nStrains <- function(diet){ sapply(plotTPS, function(day){
  isDIET <- ldata[which(ldata[, "Diet"] == diet & ldata[, "AgeAtDeath..days."] >= day),]
  mtable <- table(isDIET[, "StrainName"])
  length(grep("BXD", names(mtable)[mtable >= 3]))
})}
nInd <- sapply(plotTPS, function(day){
  isDIET <- ldata[which(ldata[, "Diet"] %in% c("CD", "HF") & grepl("BXD", ldata[, "StrainName"]) & ldata[, "AgeAtDeath..days."] >= day),]
  sum(isDIET[, "StrainName"] %in% colnames(geno))
})
nBW <- function(diet, lods){
  if(!file.exists("output/bodyweightInterpolatedData.txt")) return(NULL)
  BWs <- read.table("output/bodyweightInterpolatedData.txt", sep = "\t", quote = "", check.names = FALSE)
  bxd <- bdata[grep("^BXD", bdata[, "StrainName"]),]
  sapply(plotTPS, function(day){
    if(all(is.na(lods[, paste0("X", day)]))) return("")
    ids <- intersect(rownames(bxd)[bxd[, "Diet"] == diet], rownames(BWs)[!is.na(BWs[, as.character(day)])])
    mtable <- table(bxd[ids, "StrainName"])
    length(mtable[mtable >= 3])
  })
}

lifespanModel <- c("Model: lm(strain mean ~ set-up class + genotype, weights = mice per strain), F test of genotype. Strains with at least 3 mice alive at T; strain mean of age at death, rounded to days; set-up class: the strain's mean set-up age, in 180-day classes; genotype B −1, H 0, D +1. No kinship.",
                   mapping("lifespan.R"))
individualModel <- function(term){ c(paste0("Model: lmer(age at death ~ diet + vivarium + set-up class + genotype", if(term == "interaction") " + diet × genotype" else "",
                                            " + (1 | strain)), maximum likelihood; likelihood-ratio test of ", if(term == "interaction") "diet × genotype" else "genotype",
                                            ". BXD mice on chow and high fat alive at T; set-up class: the mouse's set-up age, in 180-day classes; genotype B −1, H 0, D +1; random strain intercept, no kinship matrix."),
                                     mapping("interactions.R")) }
bwModel <- c("Model: lm(strain mean body weight ~ genotype, weights = mice per strain), F test. Each mouse's weight interpolated linearly between its weighings; BXD strains with at least 3 mice weighed at T; T-ages with more than 400 mice weighed; genotype B −1, H 0, D +1. No kinship.",
             mapping("bodyweight.R"))
ctlModel <- c("Model: rank correlation of body weight at T with lifespan (both adjusted for vivarium and set-up class) in B and in D mice, each mouse weighted 1 / mice of its strain; χ² test (1 d.f.) of the difference on Fisher's z. Mice set up by T and alive at T; markers with at least 20 mice of each genotype.",
              paste(mapping("CTL.R"), "Inputs from bodyweight.R and adjust.R, run by run_CTL.R."))
qtl2Model <- function(m){
  covs <- if(m %in% c("main", "int")) "diet + vivarium + set-up class" else "vivarium + set-up class"
  c(paste0("Model: linear mixed model with LOCO kinship (REML): age at death ~ ", covs, " + genotype", if(m == "int") " + diet × genotype" else "",
           " + polygenic effect; likelihood-ratio test of ", if(m == "int") "diet × genotype" else "genotype", " (1 d.f.). BXD mice ",
           c(CD = "on chow", HF = "on high fat", main = "on both diets", int = "on both diets")[m],
           " alive at T; set-up class: set-up age in 180-day classes; kinship: allele sharing of strain genotype probabilities."),
    mapping("qtl2_individuals.R", "robwwilliams/BxDLifeSpan"))
}
qtl2N <- function(m){ if(!file.exists("output/qtl2_N_ind_Progressive.txt")) return(NULL); n <- read.table("output/qtl2_N_ind_Progressive.txt", sep = "\t"); n[match(plotTPS, as.numeric(rownames(n))), c(CD = "CD", HF = "HF", main = "main", int = "main")[m]] }
withSource <- function(model, lods){ c(model[1], paste(model[2], "Results:", attr(lods, "source"), "—", drawn)) }

lifespan.signs <- c("*D*: longer life", "*B*: longer life")
panels <- list(
  lifespan_CD = list(lods = readResult("output/CD_wMeans_Y_Progressive.txt", xlsx.lifespan, "CD"), N = readN("output/CD_wMeans_Y_Progressive_N.txt", nStrains("CD")), N.title = "Strains",
                     main = "Lifespan: strain means on chow diet", model = lifespanModel, sign.labels = lifespan.signs),
  lifespan_HF = list(lods = readResult("output/HF_wMeans_Y_Progressive.txt", xlsx.lifespan, "HF"), N = readN("output/HF_wMeans_Y_Progressive_N.txt", nStrains("HF")), N.title = "Strains",
                     main = "Lifespan: strain means on high-fat diet", model = lifespanModel, sign.labels = lifespan.signs),
  lifespan_main = list(lods = readResult("output/main_ind_Progressive.txt", xlsx.lifespan, "Gene (Main)"), N = readN("output/N_ind_Progressive.txt", nInd), N.title = "Mice",
                       main = "Lifespan of individual mice on both diets: main genetic effect", model = individualModel("main"), sign.labels = lifespan.signs),
  lifespan_int = list(lods = readResult("output/int_ind_Progressive.txt", xlsx.lifespan, "Gene x Diet (Interaction)"), N = readN("output/N_ind_Progressive.txt", nInd), N.title = "Mice",
                      main = "Lifespan of individual mice on both diets: genotype × diet interaction", model = individualModel("interaction"),
                      sign.labels = c("*D* − *B* larger on HF", "*D* − *B* smaller on HF")),
  BW_CD = list(lods = readResult("output/BW_CD_wMeans_Progressive.txt", xlsx.bw, "CD"), N.title = "Strains", main = "Body weight: strain means on chow diet", model = bwModel,
               sign.labels = c("*D*: heavier", "*B*: heavier")),
  BW_HF = list(lods = readResult("output/BW_HF_wMeans_Progressive.txt", xlsx.bw, "HF"), N.title = "Strains", main = "Body weight: strain means on high-fat diet", model = bwModel,
               sign.labels = c("*D*: heavier", "*B*: heavier")),
  CTL_CD = list(lods = readResult("output/CTL_CD_LODmatrix.txt", "", ""), main = "Correlated trait loci (CTL), body weight with lifespan: chow diet", model = ctlModel,
                sign.labels = c("*r*(*D*) > *r*(*B*)", "*r*(*B*) > *r*(*D*)")),
  CTL_HF = list(lods = readResult("output/CTL_HF_LODmatrix.txt", "", ""), main = "Correlated trait loci (CTL), body weight with lifespan: high-fat diet", model = ctlModel,
                sign.labels = c("*r*(*D*) > *r*(*B*)", "*r*(*B*) > *r*(*D*)")),
  qtl2_CD = list(lods = readResult("output/qtl2_CD_ind_Progressive.txt", "", ""), N = qtl2N("CD"), N.title = "Mice", model = qtl2Model("CD"), sign.labels = lifespan.signs,
                 main = "Lifespan of individual mice on chow diet: R/qtl2 with LOCO kinship"),
  qtl2_HF = list(lods = readResult("output/qtl2_HF_ind_Progressive.txt", "", ""), N = qtl2N("HF"), N.title = "Mice", model = qtl2Model("HF"), sign.labels = lifespan.signs,
                 main = "Lifespan of individual mice on high-fat diet: R/qtl2 with LOCO kinship"),
  qtl2_main = list(lods = readResult("output/qtl2_main_ind_Progressive.txt", "", ""), N = qtl2N("main"), N.title = "Mice", model = qtl2Model("main"), sign.labels = lifespan.signs,
                   main = "Lifespan of individual mice on both diets: main genetic effect, R/qtl2 with LOCO kinship"),
  qtl2_int = list(lods = readResult("output/qtl2_int_ind_Progressive.txt", "", ""), N = qtl2N("int"), N.title = "Mice", model = qtl2Model("int"),
                  sign.labels = c("*D* − *B* larger on HF", "*D* − *B* smaller on HF"),
                  main = "Lifespan of individual mice on both diets: genotype × diet interaction, R/qtl2 with LOCO kinship")
)
missing <- names(panels)[sapply(panels, function(p) is.null(p$lods))]
if(length(missing) > 0) cat("No results for", paste(missing, collapse = ", "), "\n")
panels <- panels[!names(panels) %in% missing]
for(nm in names(panels)){
  p <- panels[[nm]]
  if(grepl("^CTL", nm)) attr(p$lods, "source") <- paste0("CTL.R run in this fork (", format(file.mtime(paste0("output/", nm, "_LODmatrix.txt")), "%d %B %Y"), ")")
  p$model <- withSource(p$model, p$lods)
  p$lods <- byTPS(p$lods)
  if(grepl("^BW", nm)) p$N <- nBW(sub("BW_", "", nm), p$lods)
  panels[[nm]] <- p
}

for(nm in names(panels)){ actuaryPage(paste0("output/maps/", nm, ".svg"), panels[nm], map = map, layout = layout) }
pages <- list(actuaryQTL_lifespan_wMeans = c("lifespan_CD", "lifespan_HF"), actuaryQTL_lifespan_main_interactions = c("lifespan_main", "lifespan_int"),
              actuaryQTL_BW = c("BW_CD", "BW_HF"), actuaryCTL_lifespan = c("CTL_CD", "CTL_HF"),
              Figure2_lifespan = c("lifespan_CD", "lifespan_HF", "lifespan_main", "lifespan_int"), Figure3_bodyweight_CTL = c("BW_CD", "BW_HF", "CTL_CD", "CTL_HF"),
              Figure2_lifespan_kinship = c("qtl2_CD", "qtl2_HF", "qtl2_main", "qtl2_int"))
for(pg in names(pages)){
  if(all(pages[[pg]] %in% names(panels))){ actuaryPage(paste0("output/maps/", pg, ".svg"), panels[pages[[pg]]], map = map, layout = layout) }
}
# Everything on one page, each of Danny's lifespan maps followed by its kinship counterpart, for alignment by eye
all.order <- c("lifespan_CD", "qtl2_CD", "lifespan_HF", "qtl2_HF", "lifespan_main", "qtl2_main", "lifespan_int", "qtl2_int", "BW_CD", "BW_HF", "CTL_CD", "CTL_HF")
actuaryPage("output/maps/all_maps.svg", panels[intersect(all.order, names(panels))], map = map, layout = layout)
cat("Wrote", length(panels), "panels to output/maps/\n")
