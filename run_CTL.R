#
# run_CTL.R - Danny's CTL inputs and CTL maps, run in order
#
# Rob Williams's fork, 21 November 2026. CTL.R needs output/bodyweightInterpolatedData.txt, which bodyweight.R writes
# before it maps (lines 1-91). This runs those lines, then CTL.R, so that figure_maps.R can draw the CTL panels.
#   Rscript run_CTL.R
#
pdf("output/run_CTL_diagnostics.pdf")
eval(parse(text = readLines("bodyweight.R")[1:91]))
source("CTL.R")
invisible(dev.off())
