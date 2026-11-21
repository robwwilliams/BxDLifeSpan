#
# figure_maps.R - Danny's actuarial maps (Fig. 2) redrawn with one fixed geometry (actuary_map.R), so that every panel
# superimposes with pixel accuracy
#
# Rob Williams's fork, 21 November 2026. Reads Danny's result files in output/ (from lifespan.R, interactions.R,
# bodyweight.R and CTL.R) when they are there. Otherwise reads the two workbooks Danny sent on 5 November 2026, which
# hold the same results: give their paths.
#   Rscript figure_maps.R [QTL_Lifespan_05Nov2026.xlsx QTL_BW_5Nov.xlsx]
# The numbers on the right axis are Danny's _N files when present; otherwise they are counted here with his rules: BXD
# strains with at least 3 mice alive at T (strain means), and BXD mice on CD and HF alive at T (individuals).
# Writes output/maps/: one SVG per panel (24 x 5.5 inches), the pages of two as Danny drew them, and all panels on one page
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

# A result as markers x T-ages, columns X20 ... X770: Danny's text file if it exists, else the workbook sheet
readResult <- function(txt, xlsx, sheet){
  if(file.exists(txt)){ return(read.table(txt, sep = "\t")) }
  if(!file.exists(xlsx)) return(NULL)
  x <- as.data.frame(suppressMessages(read_excel(xlsx, sheet = sheet, col_names = FALSE)))
  hdr <- unlist(x[1, ])
  res <- suppressWarnings(sapply(x[-1, -(1:5)], as.numeric))
  colnames(res) <- paste0("X", hdr[-(1:5)])
  res <- res[, which(colSums(!is.na(res)) > 0 & !duplicated(colnames(res))), drop = FALSE]   # the repeated last header is empty
  rownames(res) <- x[-1, 1]
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

panels <- list(
  lifespan_CD = list(lods = readResult("output/CD_wMeans_Y_Progressive.txt", xlsx.lifespan, "CD"),
                     N = readN("output/CD_wMeans_Y_Progressive_N.txt", nStrains("CD")), main = "Lifespan: Strain means on chow diet"),
  lifespan_HF = list(lods = readResult("output/HF_wMeans_Y_Progressive.txt", xlsx.lifespan, "HF"),
                     N = readN("output/HF_wMeans_Y_Progressive_N.txt", nStrains("HF")), main = "Lifespan: Strain means on high fat diet"),
  lifespan_main = list(lods = readResult("output/main_ind_Progressive.txt", xlsx.lifespan, "Gene (Main)"),
                       N = readN("output/N_ind_Progressive.txt", nInd), main = "Lifespan of individuals on both diets mapping combined main genetic effects"),
  lifespan_int = list(lods = readResult("output/int_ind_Progressive.txt", xlsx.lifespan, "Gene x Diet (Interaction)"),
                      N = readN("output/N_ind_Progressive.txt", nInd), main = "Lifespan of individuals mapping locus by diet interaction effects"),
  BW_CD = list(lods = readResult("output/BW_CD_wMeans_Progressive.txt", xlsx.bw, "CD"), main = "Body weight of individuals on chow diet"),
  BW_HF = list(lods = readResult("output/BW_HF_wMeans_Progressive.txt", xlsx.bw, "HF"), main = "Body weight of individuals on high fat diet"),
  CTL_CD = list(lods = readResult("output/CTL_CD_LODmatrix.txt", "", ""), main = "CTL: Chow Diet", mb.ticks = TRUE),
  CTL_HF = list(lods = readResult("output/CTL_HF_LODmatrix.txt", "", ""), main = "CTL: High Fat Diet", mb.ticks = TRUE)
)
missing <- names(panels)[sapply(panels, function(p) is.null(p$lods))]
if(length(missing) > 0) cat("No results for", paste(missing, collapse = ", "), "\n")
panels <- panels[!names(panels) %in% missing]
panels <- lapply(panels, function(p){ p$lods <- byTPS(p$lods); p })

for(nm in names(panels)){ actuaryPage(paste0("output/maps/", nm, ".svg"), panels[nm], map = map, layout = layout) }
pages <- list(actuaryQTL_lifespan_wMeans = c("lifespan_CD", "lifespan_HF"), actuaryQTL_lifespan_main_interactions = c("lifespan_main", "lifespan_int"),
              actuaryQTL_BW = c("BW_CD", "BW_HF"), actuaryCTL_lifespan = c("CTL_CD", "CTL_HF"))
for(pg in names(pages)){
  if(all(pages[[pg]] %in% names(panels))){ actuaryPage(paste0("output/maps/", pg, ".svg"), panels[pages[[pg]]], map = map, layout = layout) }
}
actuaryPage("output/maps/all_maps.svg", panels, map = map, layout = layout)
cat("Wrote", length(panels), "panels to output/maps/\n")
cat("Right-axis numbers:\n"); print(rbind(T = plotTPS, CD = panels$lifespan_CD$N, HF = panels$lifespan_HF$N, individuals = panels$lifespan_main$N))
