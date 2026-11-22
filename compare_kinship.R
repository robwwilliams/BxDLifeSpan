#
# compare_kinship.R - Danny's maps of individual mice (random strain intercept) against R/qtl2 with LOCO kinship
#
# Rob Williams's fork, 21 November 2026. Exploratory. Reads Danny's results (output/main_ind_Progressive.txt and
# int_ind_Progressive.txt, or the sheets "Gene (Main)" and "Gene x Diet (Interaction)" of his 5 November workbook) and
# the maps of qtl2_individuals.R. Reports, at each T-age, the correlation of the two maps over markers, signed and as
# |-log10 P|, the kinship heritability, and the strongest marker of each map.
#   Rscript compare_kinship.R [QTL_Lifespan_05Nov2026.xlsx]
# Writes output/compare_kinship.txt
#

library(readxl)
source("shared.R")

args <- commandArgs(trailingOnly = TRUE)
xlsx <- if(length(args) >= 1) args[1] else "QTL_Lifespan_05Nov2026.xlsx"
tps <- seq(20, 770, 30)

readDanny <- function(txt, sheet){
  if(file.exists(txt)) return(read.table(txt, sep = "\t"))
  x <- as.data.frame(suppressMessages(read_excel(xlsx, sheet = sheet, col_names = FALSE)))
  res <- suppressWarnings(sapply(x[-1, -(1:5)], as.numeric))
  colnames(res) <- paste0("X", unlist(x[1, -(1:5)]))
  res <- res[, which(colSums(!is.na(res)) > 0 & !duplicated(colnames(res))), drop = FALSE]
  rownames(res) <- x[-1, 1]
  return(res)
}
peak <- function(v){ i <- which.max(abs(v)); sprintf("Chr %s %.1f Mb %s %.2f", map[i, "Chr"], map[i, "Mb"], if(v[i] > 0) "D" else "B", abs(v[i])) }

hsq <- read.table("output/qtl2_hsq.txt", sep = "\t")
sink("output/compare_kinship.txt", split = TRUE)
cat("Danny's maps of individual mice (lmer, random strain intercept) against R/qtl2 with LOCO kinship (qtl2_individuals.R).\n")
cat("r: correlation over markers of signed -log10 P; r.abs: of |-log10 P|. hsq: kinship heritability of the R/qtl2 null model (chromosome 1 left out).\n")
cat("Peaks: strongest marker, with the allele that raises lifespan (for the interaction: D when D - B is larger on high fat).\n")
for(m in c("main", "int")){
  danny <- readDanny(paste0("output/", m, "_ind_Progressive.txt"), c(main = "Gene (Main)", int = "Gene x Diet (Interaction)")[m])
  q <- read.table(paste0("output/qtl2_", m, "_ind_Progressive.txt"), sep = "\t")
  d <- danny[rownames(map), paste0("X", tps)]; k <- q[rownames(map), paste0("X", tps)]
  cat("\n=====", c(main = "Main genetic effect, both diets", int = "Genotype x diet interaction")[m], "\n")
  print(data.frame(T = tps, hsq = round(hsq[as.character(tps), "main"], 3),
                   r = sapply(seq_along(tps), function(i) round(cor(d[, i], k[, i], use = "complete"), 3)),
                   r.abs = sapply(seq_along(tps), function(i) round(cor(abs(d[, i]), abs(k[, i]), use = "complete"), 3)),
                   peak.Danny = sapply(seq_along(tps), function(i) peak(d[, i])),
                   peak.kinship = sapply(seq_along(tps), function(i) peak(k[, i]))), row.names = FALSE)
}
cat("\n===== Kinship heritability of each R/qtl2 null model, by T-age\n")
print(round(hsq, 3))
sink()
