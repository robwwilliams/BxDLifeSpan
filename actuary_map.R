#
# actuary_map.R - Actuarial heat maps (markers x T-ages) with one fixed geometry, so that panels from every script
# superimpose with pixel accuracy
#
# Rob Williams's fork, 21 November 2026. Danny's drawing code from QTL_plots.R, interaction_plot.R, bodyweight.R and
# CTL_plots.R, in one function. His warm colors, bins, row lines and chromosome boxes are kept. What is fixed:
#   - Margins in inches (par(mai)), the same for every panel, so that every panel has exactly the same chromosome area.
#     (QTL_plots.R used R's default right margin of 2.1 lines, interaction_plot.R 4.1 lines.)
#   - The x layout (chromosome lengths from the full map, 15 Mb gaps) and the rows (plotTPS, whatever T-ages a result
#     holds; CTL_plots.R drew its matrices without byTPS).
#   - Text size: cex = 1 on every page, however many panels it has.
#   - A white background over the whole page, so that the bounding box Illustrator sees is the page itself.
#   - Runs of neighboring markers with the same color are drawn as one rectangle: the same picture, a smaller file, and
#     no hairline seams between rectangles.
# Rob's changes of 21 November 2026:
#   - Titles left-justified and in plain words; each panel's model, mapping code and results file printed small under it
#     (10 pt, for working copies; in the paper it belongs in the legend).
#   - The legend in the right margin, clear of the maps, saying which allele each color favors, with B, D, P and r in
#     italic (written *B* in labels).
#   - Type sizes from Rob's working Illustrator file: titles 24 pt bold, chromosome numbers 30 pt. Axis numbers and the
#     legend are 17 pt and the axis titles 23 pt, so that all are at least 5 pt when a 24-inch panel is scaled to
#     Nature's 183 mm.
#   - Ticks every 25 Mb along each chromosome, longer every 50 Mb, without numbers.
#   - The B side shaded from pale aqua through green and blue to navy, so that its five bins differ in hue as well as in
#     lightness, as the warm side does (Danny's were five blues).
#   - The right-axis title (Strains, Mice) beside its numbers, reading top to bottom, so that it is not read as part of
#     the legend.
#   - Nature-style panel letters (a, b, ...) at the far left of each panel on a page of several; ticks along the top of
#     each chromosome as well as the bottom; chromosome numbers set lower, clear of the ticks; the model note close under
#     its own panel; and, where permutations have given one, the genome-wide threshold as the first line under the panel.
# SVG through svglite (Arial, each label an editable text element, no text stretching), or PDF when the file name ends
# in .pdf
#

library(svglite)

plotTPS <- seq(50, 770, 30)
labelTPS <- seq(50, 770, 60)
map.mai <- c(1.45, 1.5, 0.7, 3.5)                         # bottom, left, top, right (inches)
map.colz.c <- c("white", "#fceea7", "#fdb863", "#e66101", "#d01c8b", "#c42e2f")
map.colz.n <- c("white", "#bfe6d8", "#66c2a4", "#2b8cbe", "#253494", "#081d58")
map.colz.g <- c("white", "#e3e3e3", "#bdbdbd", "#8c8c8c", "#595959", "#262626")   # grays: the UM-HET3 founder not shared with BXD (C or H)
map.thresholds <- c(0, 1, 2, 3, 4, 5, 100)
map.bins <- c("1–2", "2–3", "3–4", "4–5", ">5")
# A wider scale, for maps whose peaks run far above 5 (UM-HET3: Rob, 24 November 2026, "We will need more dynamic range on
# the P value colors"). The first four bins and colors are the same as above, so the two scales read alike up to 5.
map.thresholds.wide <- c(0, 1, 2, 3, 4, 5, 6, 8, 100)
map.bins.wide <- c("1–2", "2–3", "3–4", "4–5", "5–6", "6–8", ">8")
map.colz.c.wide <- c(map.colz.c[1:5], "#c42e2f", "#8e1622", "#4a0a12")
map.colz.n.wide <- c(map.colz.n[1:5], "#081d58", "#2c0b5e", "#12002e")
map.colz.g.wide <- c("white", "#ebebeb", "#cccccc", "#a8a8a8", "#858585", "#636363", "#3d3d3d", "#141414")
map.cex <- 1.42                                            # axis numbers and legend, 17 pt
map.cex.title <- 1.9                                       # axis titles, 23 pt

# Chromosome layout from the full genetic map (Mb in bp, as the plot scripts scale it), shared by every panel
mapLayout <- function(map, chrs = c(1:19, "X"), gap = 15000000){
  chrs.length <- unlist(lapply(chrs, function(x){ max(map[map[, "Chr"] == x, "Mb"]) }))
  names(chrs.length) <- chrs
  chr.x <- c(0, cumsum(chrs.length + gap))[1:length(chrs)]
  names(chr.x) <- chrs
  return(list(chrs = chrs, length = chrs.length, start = chr.x, l.x = sum(chrs.length) + gap * (length(chrs)-1)))
}

# Position in user coordinates of a point given in inches from the plot region's right and top edges
inchX <- function(dx){ grconvertX(grconvertX(par("usr")[2], "user", "inches") + dx, "inches", "user") }
inchY <- function(dy){ grconvertY(grconvertY(par("usr")[4], "user", "inches") + dy, "inches", "user") }

# Text with italic pieces, left-aligned at (x, y): "*D*: longer life" sets D in italic. One text element per piece
richText <- function(x, y, label, cex = map.cex, font = 1, vadj = 1){
  pieces <- strsplit(label, "*", fixed = TRUE)[[1]]
  for(i in seq_along(pieces)){
    if(pieces[i] == "") next
    f <- if(i %% 2 == 0) (if(font == 2) 4 else 3) else font
    lead <- sub("^( *).*$", "\\1", pieces[i])                                  # SVG drops spaces at the ends of a text element,
    text(x + strwidth(lead, cex = cex, font = f), y, trimws(pieces[i]), adj = c(0, vadj), cex = cex, font = f, xpd = NA)   # so place them as offsets
    x <- x + strwidth(pieces[i], cex = cex, font = f)
  }
}

# Width of such a label, in user coordinates
richWidth <- function(label, cex = map.cex){
  pieces <- strsplit(label, "*", fixed = TRUE)[[1]]
  sum(sapply(seq_along(pieces), function(i) strwidth(pieces[i], cex = cex, font = if(i %% 2 == 0) 3 else 1)))
}

# The legend, in the right margin: -log10 P bins for each sign, with what the sign means
mapLegend <- function(sign.labels, x0 = 1.3, colz = list(map.colz.c, map.colz.n), bins = map.bins){
  bw <- 0.36; bh <- 0.2
  richText(inchX(x0), inchY(0), "−log10 *P*", font = 2)
  y <- -0.3
  for(s in 1:2){
    richText(inchX(x0), inchY(y), sign.labels[s])
    y <- y - 0.3
    for(b in seq_along(bins)){
      rect(inchX(x0), inchY(y - bh), inchX(x0 + bw), inchY(y), col = colz[[s]][b + 1], border = "grey40", lwd = 0.5, xpd = NA)
      text(inchX(x0 + bw + 0.1), inchY(y - bh / 2), bins[b], adj = c(0, 0.5), cex = map.cex, xpd = NA)
      y <- y - bh - 0.035
    }
    y <- y - 0.08
  }
  text(inchX(x0), inchY(y), "White: below 1", adj = c(0, 1), cex = map.cex, xpd = NA)
}

# One panel: lods is markers x T-ages (columns X50 ... X770, as byTPS returns), N the numbers on the right axis
actuaryMap <- function(lods, map, layout, N = NULL, N.title = NULL, main = "", model = NULL, tps = plotTPS, threshold = 1,
                       sign.labels = c("*D* allele higher", "*B* allele higher"), chr.numbers = TRUE, letter = NULL, gw.note = NULL,
                       colz = list(map.colz.c, map.colz.n), labels = labelTPS, thresholds = map.thresholds, bins = map.bins){             # colors of positive and negative values; labeled T-ages
  lods <- lods[, paste0("X", tps), drop = FALSE]
  nr <- length(tps)
  plot(c(1, layout$l.x), c(1, 1 + nr), t = 'n', bty = 'n', xlab = "", xaxt = 'n', ylab = "", yaxt = 'n', xaxs = "i", yaxs = "i")
  title(main = main, adj = 0, line = 1.3, cex.main = 2)
  if(!is.null(letter)) mtext(letter, side = 3, line = 1.3, at = grconvertX(0.15, "inches", "user"), adj = 0, font = 2, cex = 2.3)
  abline(h = 1:nr, col = rgb(0.5, 0.5, 0.5))
  chr.pe <- 0
  for(chr in layout$chrs){
    chr.s <- layout$start[chr]; chr.l <- layout$length[chr]
    rect(chr.s, 0, chr.s + chr.l, 1 + nr, col = rgb(1, 1, 1, 0.5), border = NA)
    tk <- seq(25e6, chr.l, 25e6)
    for(side in c(1, 3)){
      axis(side, at = chr.s + tk[tk %% 50e6 != 0], labels = FALSE, tcl = -0.25, lwd = 0, lwd.ticks = 1)   # 25 Mb
      axis(side, at = chr.s + tk[tk %% 50e6 == 0], labels = FALSE, tcl = -0.5, lwd = 0, lwd.ticks = 1)    # 50 Mb
    }
    if(chr.numbers) mtext(chr, side = 1, line = 1.7, at = chr.s + chr.l / 2, cex = 2.5)
    rect(chr.pe, 0, chr.s, 1 + nr, col = "white", border = NA)
    abline(v = chr.s)
    abline(v = chr.s + chr.l)
    rr <- rownames(map)[which(map[, "Chr"] == chr)]
    rr <- rr[order(map[rr, "Mb"])]
    rr.pos <- map[rr, "Mb"] + chr.s
    x.s <- c(chr.s, (rr.pos[-1] + rr.pos[-length(rr)]) / 2)                  # each marker reaches halfway to its neighbors
    x.e <- c((rr.pos[-1] + rr.pos[-length(rr)]) / 2, chr.s + chr.l)
    for(tp in 1:nr){
      lod <- as.numeric(lods[rr, tp])
      cN <- pmin(findInterval(abs(lod), thresholds, left.open = TRUE), length(bins) + 1)  # 2 for 1 - 2, ..., 6 for > 5, as Danny's bins
      cN[is.na(lod) | abs(lod) <= threshold] <- 0
      code <- sign(lod) * cN; code[is.na(code)] <- 0
      run <- cumsum(c(TRUE, diff(code) != 0))
      for(g in unique(run[code != 0])){
        i <- which(run == g)
        rect(x.s[min(i)], tp, x.e[max(i)], tp + 1, col = colz[[if(code[i[1]] < 0) 2 else 1]][abs(code[i[1]])], border = NA)
      }
    }
    chr.pe <- chr.s + chr.l
  }
  rect(layout$start, 1, layout$start + layout$length, 1 + nr, lwd = 2, xpd = NA)
  at <- 0.5 + 1:nr
  axis(2, at = at, labels = FALSE, tcl = -0.25)
  axis(2, at = at[match(labels, tps)], labels, las = 2, cex.axis = map.cex, tcl = -0.5, mgp = c(3, 0.6, 0))
  mtext("T-age (days)", side = 2, line = 3.9, cex = map.cex.title)
  if(!is.null(N)){
    labs <- N[match(labels, tps)]
    axis(4, at = at, labels = FALSE, tcl = -0.25)
    axis(4, at = at[match(labels, tps)], labs, las = 2, cex.axis = map.cex, tcl = -0.5, mgp = c(3, 0.6, 0))
    if(!is.null(N.title)){
      w <- max(strwidth(as.character(labs), units = "inches", cex = map.cex))      # beside the widest number
      text(inchX(0.22 + w + 0.2), mean(par("usr")[3:4]), N.title, srt = -90, cex = map.cex.title, xpd = NA)
    }
  }
  mapLegend(sign.labels, colz = colz, bins = bins)
  line0 <- if(chr.numbers) 4.3 else 1.1                                          # just under its own panel
  if(!is.null(gw.note)){                                                          # the genome-wide threshold first, larger than the model note
    richText(par("usr")[1], grconvertY(grconvertY(par("usr")[3], "user", "inches") - 0.2 * (line0 + 0.3), "inches", "user"), gw.note, cex = 1.1, vadj = 1)
    line0 <- line0 + 1.5
  }
  if(!is.null(model)){
    model <- unlist(lapply(model, strwrap, width = 320))                        # about 330 characters fit on a line; two lines leave room above the next title
    for(l in seq_along(model)){ mtext(model[l], side = 1, line = line0 + 1.05 * (l - 1), adj = 0, cex = 0.8, col = "grey25") }
  }
}

# A page of panels, each a list of actuaryMap's arguments (lods, N, main, ...). Every panel is 24 x 5.5 inches. On a page of
# several, chromosome numbers are under the bottom panel only, as in Rob's working file (bottom.only = FALSE puts them
# under every panel); the geometry is the same either way
actuaryPage <- function(file, panels, map, layout = mapLayout(map), width = 24, panel.height = 5.5, bottom.only = TRUE){
  height <- panel.height * length(panels)
  if(grepl("\\.pdf$", file)){ pdf(file, width = width, height = height, bg = "white") }else{
    svglite(file, width = width, height = height, bg = "white", system_fonts = list(sans = "Arial"), fix_text_size = FALSE)
  }
  op <- par(mfrow = c(length(panels), 1))
  par(mai = map.mai, omi = c(0, 0, 0, 0), cex = 1)
  for(i in seq_along(panels)){
    lt <- if(length(panels) > 1) letters[i] else NULL                                  # Nature panel letters on a page of several
    do.call(actuaryMap, c(panels[[i]], list(map = map, layout = layout, chr.numbers = !bottom.only || i == length(panels), letter = lt)))
  }
  invisible(dev.off())
}
