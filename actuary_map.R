#
# actuary_map.R - Actuarial heat maps (markers x T-ages) with one fixed geometry, so that panels from every script
# superimpose with pixel accuracy
#
# Rob Williams's fork, 21 November 2026. Danny's drawing code from QTL_plots.R, interaction_plot.R, bodyweight.R and
# CTL_plots.R, in one function. His colours, bins, row lines and chromosome boxes are kept. What is fixed:
#   - Margins in inches (par(mai)), the same for every panel. The right margin always has room for 4-digit Ns, so a
#     panel with 2-digit Ns, or none, has exactly the same chromosome area. (QTL_plots.R used R's default right margin of
#     2.1 lines, interaction_plot.R 4.1 lines.)
#   - The x layout (chromosome lengths from the full map, 15 Mb gaps) and the rows (plotTPS, whatever T-ages a result
#     holds; CTL_plots.R drew its matrices without byTPS).
#   - Text size: cex = 1 on every page, however many panels it has.
#   - A white background over the whole page, so that the bounding box Illustrator sees is the page itself.
#   - Runs of neighbouring markers with the same colour are drawn as one rectangle: the same picture, a smaller file, and
#     no hairline seams between rectangles.
# SVG through svglite (Arial, each label one text element), or PDF when the file name ends in .pdf
#

library(svglite)

plotTPS <- seq(50, 770, 30)
map.mai <- c(1.02, 0.82, 0.82, 0.82)                       # bottom, left, top, right (inches): mar = c(5.1, 4.1, 4.1, 4.1)
map.colz.c <- c("white", "#fceea7", "#fdb863", "#e66101", "#d01c8b", "#c42e2f")
map.colz.n <- c("white", "#c6dbef", "#6baed6", "#3182bd", "#08519c", "#08306b")
map.thresholds <- c(0, 1, 2, 3, 4, 5, 100)
map.bins <- c("0 - 1", "1 - 2", "2 - 3", "3 - 4", "4 - 5", ">5")

# Chromosome layout from the full genetic map (Mb in bp, as the plot scripts scale it), shared by every panel
mapLayout <- function(map, chrs = c(1:19, "X"), gap = 15000000){
  chrs.length <- unlist(lapply(chrs, function(x){ max(map[map[, "Chr"] == x, "Mb"]) }))
  names(chrs.length) <- chrs
  chr.x <- c(0, cumsum(chrs.length + gap))[1:length(chrs)]
  names(chr.x) <- chrs
  return(list(chrs = chrs, length = chrs.length, start = chr.x, l.x = sum(chrs.length) + gap * (length(chrs)-1)))
}

# One panel: lods is markers x T-ages (columns X50 ... X770, as byTPS returns), N the numbers shown on the right axis
actuaryMap <- function(lods, map, layout, N = NULL, main = "", tps = plotTPS, threshold = 1, legend = TRUE, mb.ticks = FALSE){
  lods <- lods[, paste0("X", tps), drop = FALSE]
  nr <- length(tps)
  plot(c(1, layout$l.x), c(1, 1 + nr), t = 'n', bty = 'n', xlab = "", xaxt = 'n', ylab = "", yaxt = 'n', xaxs = "i", yaxs = "i", main = main)
  abline(h = 1:nr, col = rgb(0.5, 0.5, 0.5))
  chr.pe <- 0
  for(chr in layout$chrs){
    chr.s <- layout$start[chr]; chr.l <- layout$length[chr]
    rect(chr.s, 0, chr.s + chr.l, 1 + nr, col = rgb(1, 1, 1, 0.5), border = NA)
    if(mb.ticks){
      axis(1, at = chr.s + seq(0, chr.l, 50e6), labels = seq(0, chr.l, 50e6) / 1e6, cex.axis = 0.4, mgp = c(3, 0.2, 0))
      mtext(chr, side = 1, line = 1.3, at = chr.s + chr.l / 2)
    }else{
      axis(1, at = chr.s + chr.l / 2, chr)
    }
    rect(chr.pe, 0, chr.s, 1 + nr, col = "white", border = NA)
    abline(v = chr.s)
    abline(v = chr.s + chr.l)
    rr <- rownames(map)[which(map[, "Chr"] == chr)]
    rr <- rr[order(map[rr, "Mb"])]
    rr.pos <- map[rr, "Mb"] + chr.s
    x.s <- c(chr.s, (rr.pos[-1] + rr.pos[-length(rr)]) / 2)                  # each marker reaches halfway to its neighbours
    x.e <- c((rr.pos[-1] + rr.pos[-length(rr)]) / 2, chr.s + chr.l)
    for(tp in 1:nr){
      lod <- as.numeric(lods[rr, tp])
      cN <- findInterval(abs(lod), map.thresholds, left.open = TRUE)        # 2 for 1 - 2, ..., 6 for > 5, as Danny's max(which(thresholds < abs(lod)))
      cN[is.na(lod) | abs(lod) <= threshold] <- 0
      code <- sign(lod) * cN; code[is.na(code)] <- 0
      run <- cumsum(c(TRUE, diff(code) != 0))
      for(g in unique(run[code != 0])){
        i <- which(run == g)
        rect(x.s[min(i)], tp, x.e[max(i)], tp + 1, col = (if(code[i[1]] < 0) map.colz.n else map.colz.c)[abs(code[i[1]])], border = NA)
      }
    }
    chr.pe <- chr.s + chr.l
  }
  rect(layout$start, 1, layout$start + layout$length, 1 + nr, lwd = 2, xpd = NA)
  axis(2, at = 0.5 + 1:nr, tps, las = 2)
  if(!is.null(N)) axis(4, at = 0.5 + 1:nr, N, las = 2)
  if(legend) legend("topright", fill = c(map.colz.c, map.colz.n), c(map.bins, map.bins), ncol = 2, bg = "white", title = "-log10P")
}

# A page of panels, each a list of actuaryMap's arguments (lods, N, main, ...). Every page has the same panel size:
# 24 x 5.5 inches per panel, as in Danny's 24 x 11 inch pages of two
actuaryPage <- function(file, panels, map, layout = mapLayout(map), width = 24, panel.height = 5.5){
  height <- panel.height * length(panels)
  if(grepl("\\.pdf$", file)){ pdf(file, width = width, height = height, bg = "white") }else{
    svglite(file, width = width, height = height, bg = "white", system_fonts = list(sans = "Arial"))
  }
  op <- par(mfrow = c(length(panels), 1))
  par(mai = map.mai, omi = c(0, 0, 0, 0), cex = 1)
  for(p in panels){ do.call(actuaryMap, c(p, list(map = map, layout = layout))) }
  invisible(dev.off())
}
