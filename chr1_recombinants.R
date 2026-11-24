#
# chr1_recombinants.R - BXD strains with recombinations in the two proximal Chr 1 intervals, and whether they have
# lifespan data: the strains (and F1s) that could split the intervals
#
# Rob Williams's fork, 23 November 2026. Exploratory. Rob: "Wish we had more BXDs with strategically placed recombinations.
# Be my guest to look for any you think we should add (although that is a 3 years project ;-)".
#   - Chow interval: at T770 every marker from 4.86 to 5.32 Mb has the same genotypes in all the strains and F1 types
#     mapped (chr1_zoom.R), so the peak is a flat block; the 1.5-LOD interval is 4.86-5.52 Mb with the F1s. Genes (mm10
#     RefSeq): Tcea1, Rgs20, Atp6v1h inside the block; Mrpl15 and Lypla1 just proximal, Oprk1 just distal.
#   - High-fat interval: the T500 peak at 37.6 Mb, with markers within 1.5 LOD from about 34.4 to 38.7 Mb.
#   - For every BXD strain in the genotype file (mm10): each switch between B and D within a window around each
#     interval, placed between its two flanking markers; and the strain's number of mice with lifespans on chow and on
#     high fat (shared.R). Heterozygous (H) and unknown calls are skipped when finding switches.
#   - Strains on the extinct list of the genotype file's header (as of December 2016) are flagged; that list may be
#     out of date.
# Writes output/chr1_recombinants.txt
#   Rscript chr1_recombinants.R
#

source("shared.R")
hdr <- readLines("data/BXD_current_rev050423.geno", n = 23)
extinct <- strsplit(sub(".*including: ", "", grep("extinct as of", hdr, value = TRUE)), " ")[[1]]
extinct <- trimws(gsub("\t", "", extinct))
bxd <- grep("^BXD[0-9]+[a-z]?$", colnames(geno), value = TRUE)
nCD <- table(factor(ldata[ldata$Diet == "CD", "StrainName"], levels = bxd))
nHF <- table(factor(ldata[ldata$Diet == "HF", "StrainName"], levels = bxd))
f1parents <- c("BXD9", "BXD170", "BXD172", "BXD197", "BXD65", "BXD73b", "BXD77", "BXD102")

windows <- list(chow = c(4.0, 7.0, 4.8466, 5.5266), highfat = c(30.0, 42.0, 34.4, 38.7))   # window; interval of interest
blockChow <- c(4.8583, 5.3155)                                                     # the flat peak block at chow T770

switches <- function(s, w){
  k <- rownames(map)[map$Chr == "1" & map$Mb >= w[1] & map$Mb <= w[2]]; k <- k[order(map[k, "Mb"])]
  g <- as.character(geno[k, s]); ok <- g %in% c("B", "D"); g <- g[ok]; mb <- map[k, "Mb"][ok]
  i <- which(g[-1] != g[-length(g)])
  if(length(i) == 0) return(NULL)
  data.frame(strain = s, from = mb[i], to = mb[i + 1], change = paste0(g[i], ">", g[i + 1]))
}
sink("output/chr1_recombinants.txt", split = TRUE)
cat("BXD strains with a B/D switch near the two proximal Chr 1 intervals (mm10). from-to: the flanking markers of the switch (Mb).\n")
cat("CD, HF: mice with lifespans on chow and on high fat (all ages). F1 parent: a parent of the 13 F1 types. Extinct: on the\n")
cat("genotype file's list of extinct strains (December 2016; may be out of date).\n")
for(w in names(windows)){
  W <- windows[[w]]
  r <- do.call(rbind, lapply(bxd, switches, w = W[1:2]))
  r <- r[r$to > W[3] & r$from < W[4], ]                                             # switches that could fall inside the interval
  r$CD <- as.numeric(nCD[r$strain]); r$HF <- as.numeric(nHF[r$strain])
  r$F1.parent <- ifelse(r$strain %in% f1parents, "yes", ""); r$extinct <- ifelse(r$strain %in% extinct, "yes", "")
  if(w == "chow") r$splits.block <- ifelse(r$to > blockChow[1] & r$from < blockChow[2], "yes", "")
  r <- r[order(r$from, r$strain), ]
  cat(sprintf("\n%s interval %.2f-%.2f Mb: %d switches in %d strains\n", w, W[3], W[4], nrow(r), length(unique(r$strain))))
  r$from <- round(r$from, 3); r$to <- round(r$to, 3)
  print(r, row.names = FALSE)
}
sink()
