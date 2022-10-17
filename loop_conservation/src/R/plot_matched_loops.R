options(warn = 1)

# Read config file
config <- yaml::yaml.load_file("config/config.yml")

library(CNEr)
library(data.table)
library(GenomicRanges)
library(ggforce)
library(ggplot2)
library(karyoploteR)

source("src/R/new_aes.R")

genome <- "D_mel"
other_genome <- "D_vir"

#
#  Read long-range loops in D. melanogaster and D. virilis
#

lap <- as.data.table(read.table("data/loops/loops_D_mel_annotated.tsv",
  header = T, sep = "\t", comment.char = "#", stringsAsFactors = F))

extract_loops <- function(genome)
{
  dt <- fread(paste0("data/loops/long_range_loops_", genome, ".tsv"), header = T, sep = "\t")

  loops <- with(dt, data.table(
    order = c(seq_len(nrow(dt)), seq_len(nrow(dt))),
    loop_id = c(loop_id, loop_id),
    aggregated_id = c(aggregated_id, aggregated_id),
    anchor = c(rep("A1", nrow(dt)), rep("A2", nrow(dt))),
    chrom = c(chr1, chr2),
    start = c(x1, y1),
    end = c(x2, y2),
    midpoint = c(as.integer((x1 + x2) / 2), as.integer((y1 + y2) / 2))
  ))
  setkey(loops, order)
  loops[, order := NULL]

  return(loops)
}

lap_other <- extract_loops(other_genome)

#
#  Read loop matching
#

candidates <- fread(paste0("data/loops/loops_", genome, "_to_", other_genome, "_candidates.tsv"), header = T, sep = "\t")

lap_correspondence <- unique(candidates[!is.na(other_loop_id), list(loop_id, other_loop_id, A1_chrom, A1_other_chrom, A2_chrom, A2_other_chrom)])
# double-check that all the loops are cis-chromosomal, and that the chromosome identity is conserved
stopifnot(with(lap_correspondence, identical(A1_chrom, A2_chrom)))
stopifnot(with(lap_correspondence, identical(A1_other_chrom, A2_other_chrom)))
stopifnot(with(lap_correspondence, identical(A1_chrom, A1_other_chrom)))

lap[, match_loop_id := ifelse(loop_id %in% lap_correspondence$loop_id, loop_id, NA)]
lap_other[, match_loop_id := lap_correspondence$loop_id[match(loop_id, lap_correspondence$other_loop_id)]]

#
#  Prepare the data on CE chains for plotting
#

load(paste0("data/lastz/ceFinalAnnotated.", genome, ".", other_genome, ".Rdata"))

cneFinal_dt <- as.data.table(cneFinal)
cneFinal_dt[, first.midpoint := (first.start + first.end) / 2]
cneFinal_dt[, second.midpoint := (second.start + second.end) / 2]

# determine the chain orientation
cneFinal_dt[, chain.first.ori := unique(head(first.midpoint, -1) < tail(first.midpoint, -1)), by = "chain_id"]
cneFinal_dt[, chain.second.ori := unique(head(second.midpoint, -1) < tail(second.midpoint, -1)), by = "chain_id"]

# if collapse_CNEs = TRUE, use 'mstart' and 'mend' instead of 'start' and 'end'
cneFinal_dt[, chain.first.start := ifelse(chain.first.ori, first.end, first.start)]
cneFinal_dt[, chain.first.end := ifelse(chain.first.ori, c(tail(first.start, -1), NA), c(tail(first.end, -1), NA)), by = "chain_id"]
cneFinal_dt[, chain.first.mstart := first.midpoint]
cneFinal_dt[, chain.first.mend := c(tail(first.midpoint, -1), NA), by = "chain_id"]

cneFinal_dt[, chain.second.start := ifelse(chain.second.ori, second.end, second.start)]
cneFinal_dt[, chain.second.end := ifelse(chain.second.ori, c(tail(second.start, -1), NA), c(tail(second.end, -1), NA)), by = "chain_id"]
cneFinal_dt[, chain.second.mstart := second.midpoint]
cneFinal_dt[, chain.second.mend := c(tail(second.midpoint, -1), NA), by = "chain_id"]

# for an overview: just draw a single segment per CNE chain
chain_dt <- cneFinal_dt[, list(
  chain.first.start = min(first.start),
  chain.first.end = max(first.end),
  chain.second.start = min(second.start),
  chain.second.end = max(second.end),
  .N
), by = list(first.seqnames, first.strand, second.seqnames, second.strand, chain_id)]
stopifnot(!anyDuplicated(chain_dt$chain_id))

#
#  Main plotting function
#

syntenic_plot <- function(range, range_other, draw_exons = TRUE, collapse_CNEs = TRUE, collapse_chains = FALSE, anchor_margin = 1e3, anchors_as_lines = TRUE, anchors_other_as_lines = FALSE, label_anchors = TRUE, breaks = waiver())
{
  stopifnot(length(range) == 1)
  stopifnot(length(range_other) == 1)
  this_chrom <- as.character(seqnames(range))
  chrom_other <- as.character(seqnames(range_other))

  anchor_dt <- lap[sub("^chr", "", anchor_chr) == this_chrom, ]
  anchor_other_dt <- lap_other[chrom == chrom_other, ]

  dt <- cneFinal_dt[first.seqnames == this_chrom & second.seqnames == chrom_other, ]
  # dt <- dt[start(range) <= first.end & first.start <= end(range) & start(range_other) <= second.end & second.start <= end(range_other), ]

  p <- ggplot(dt, aes(color = second.strand)) +
    theme_bw() +
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5), axis.text.y.right = element_text(hjust = 0.5)) +
    xlab(paste0("D. mel. chromosome ", this_chrom)) +
    ylab(paste0("D. vir. chromosome ", chrom_other)) +
    coord_fixed()

  if (draw_exons)
    p <- p +
      geom_rect(data = exons_other[seqnames == chrom_other, ], aes(ymin = start, ymax = end,
        xmin = -Inf, xmax = Inf), fill = "black", inherit.aes = FALSE, alpha = 0.3) +
      geom_rect(data = exons[seqnames == this_chrom, ], aes(xmin = start, xmax = end,
        ymin = -Inf, ymax = Inf), fill = "black", inherit.aes = FALSE, alpha = 0.3)

  if (anchors_as_lines)
    p <- p +
      geom_vline(data = anchor_dt, aes(xintercept = (anchor_start + anchor_end) / 2,
        color = match_loop_id), alpha = 0.3)
  else
    p <- p +
      geom_rect(data = anchor_dt, aes(xmin = anchor_start - anchor_margin, xmax = anchor_end + anchor_margin,
        ymin = -Inf, ymax = Inf, fill = match_loop_id), inherit.aes = FALSE, alpha = 0.3)

  if (anchors_other_as_lines)
    p <- p +
      geom_hline(data = anchor_other_dt, aes(yintercept = (start + end) / 2,
        color = match_loop_id), alpha = 0.3)
  else
    p <- p +
      geom_rect(data = anchor_other_dt, aes(ymin = start - anchor_margin, ymax = end + anchor_margin,
        xmin = -Inf, xmax = Inf, fill = match_loop_id), inherit.aes = FALSE, alpha = 0.3)

  if (any(!is.na(c(anchor_dt$match_loop_id, anchor_other_dt$match_loop_id))))
    p <- p +
      scale_fill_hue(l = 70, c = 150, na.value = "grey80", guide = "none") +
      scale_color_hue(l = 70, c = 150, na.value = "grey80", guide = "none")
  else
  # in case all the matched loop IDs are equal to NA, use the na.value color throughout
    p <- p +
      scale_fill_brewer(na.value = "grey80", guide = "none") +
      scale_color_brewer(na.value = "grey80", guide = "none")

  if (label_anchors)
    p <- p +
      scale_x_continuous(label=scales::comma, breaks = breaks, limits = c(start(range), end(range)), oob = function(x, ...) x,
        sec.axis = dup_axis(breaks = anchor_dt$anchor_summit, labels = with(anchor_dt, paste0(loop_id, '\n', anchor)), name = NULL)) +
      scale_y_continuous(label=scales::comma, breaks = breaks, limits = c(start(range_other), end(range_other)), oob = function(x, ...) x,
        sec.axis = dup_axis(breaks = anchor_other_dt$midpoint, labels = with(anchor_other_dt, paste0(loop_id, '\n', anchor)), name = NULL))
  else
    p <- p +
      scale_x_continuous(label=scales::comma, breaks = breaks, limits = c(start(range), end(range)), oob = function(x, ...) x) +
      scale_y_continuous(label=scales::comma, breaks = breaks, limits = c(start(range_other), end(range_other)), oob = function(x, ...) x)

  p <- p +
    new_scale_color() +
    scale_color_manual(values = c('+' = '#0571b0', '-' = '#ca0020', '*' = '#999999'), guide = 'none')

  if (collapse_chains)
    p <- p + geom_segment(data = chain_dt[first.seqnames == this_chrom & second.seqnames == chrom_other, ],
      aes(x = chain.first.start, xend = chain.first.end, y = ifelse(second.strand == '-', chain.second.end, chain.second.start), yend = ifelse(second.strand == '-', chain.second.start, chain.second.end)), size = 1)
  else
  {
    if (collapse_CNEs)
      p <- p +
        geom_segment(aes(x = chain.first.mstart, xend = chain.first.mend, y = chain.second.mstart, yend = chain.second.mend), size = 0.5) +
        geom_point(aes(x = first.midpoint, y = second.midpoint), size = 0.5)
    else
      p <- p +
        geom_segment(aes(x = chain.first.start, xend = chain.first.end, y = chain.second.start, yend = chain.second.end), size = 0.5) +
      geom_segment(aes(x = first.start, xend = first.end, y = ifelse(second.strand == '-', second.end, second.start), yend = ifelse(second.strand == '-', second.start, second.end)), size = 2)
  }

  return(p)
}

syntenic_plot_around_anchors <- function(loop_id, anchor, loop_id_other, anchor_other, margin = 25e3, ...)
{
  this_loop_id <- loop_id
  this_anchor <- anchor
  mid <- lap[lap$loop_id == this_loop_id & lap$anchor == this_anchor, ]
  mid <- mid[, list(midpoint = mean((anchor_start + anchor_end) / 2)), by = anchor_chr]
  range <- GRanges(sub("^chr", "", mid$anchor_chr), IRanges(mid$midpoint - margin, mid$midpoint + margin))

  mid_other <- lap_other[lap_other$loop_id == loop_id_other & lap_other$anchor == anchor_other, ]
  mid_other <- mid_other[, list(midpoint = mean((start + end) / 2)), by = chrom]
  range_other <- GRanges(mid_other$chrom, IRanges(mid_other$midpoint - margin, mid_other$midpoint + margin))

  return(syntenic_plot(range, range_other, ...))
}

syntenic_plot_chromosome <- function(chrom, chrom_other, draw_exons = FALSE, collapse_chains = TRUE, ...)
{
  range <- GRanges(chrom, IRanges(1, seqlengths(first(cneFinal))[chrom]))
  range_other <- GRanges(chrom_other, IRanges(1, seqlengths(second(cneFinal))[chrom_other]))
  return(syntenic_plot(range, range_other, draw_exons = draw_exons, collapse_chains = collapse_chains, ...))
}

syntenic_plot_beziers <- function(chrom, chrom_other, label_anchors = TRUE)
{
  range <- GRanges(chrom, IRanges(1, seqlengths(first(cneFinal))[chrom]))
  range_other <- GRanges(chrom_other, IRanges(1, seqlengths(second(cneFinal))[chrom_other]))

  stopifnot(length(range) == 1)
  stopifnot(length(range_other) == 1)

  this_chrom <- chrom
  anchor_dt <- lap[sub("^chr", "", anchor_chr) == this_chrom, ]
  anchor_other_dt <- lap_other[chrom == chrom_other, ]

  beziers <- anchor_dt[, list(x = rep(anchor_summit, each = 2), y = c(0, rep(anchor_summit[2] - anchor_summit[1], 2)^0.9 + 1e6, 0), point = c('end', 'control', 'control', 'end')), by = list(loop_id, match_loop_id)]
  beziers_other <- anchor_other_dt[, list(y = rep(midpoint, each = 2), x = c(0, rep(midpoint[2] - midpoint[1], 2)^0.9 + 1e6, 0), point = c('end', 'control', 'control', 'end')), by = list(loop_id, match_loop_id)]

  p <- ggplot() +
    theme_bw() +
    theme(axis.text.y = element_text(angle = 90, hjust = 0.5), axis.text.y.right = element_text(hjust = 0.5)) +
    xlab(paste0("D. mel. chromosome ", chrom)) +
    ylab(paste0("D. vir. chromosome ", chrom_other)) +
    coord_fixed()

  if (any(!is.na(c(beziers$match_loop_id, beziers_other$match_loop_id))))
    p <- p + scale_color_hue(l = 70, c = 150, na.value = "grey80", guide = "none")
  else
    # in case all the matched loop IDs are equal to NA, use the na.value color throughout
    p <- p + scale_color_brewer(na.value = "grey80", guide = "none")

  if (label_anchors)
    p <- p +
      scale_x_continuous(label=scales::comma, limits = c(start(range), end(range)), oob = function(x, ...) x,
        sec.axis = dup_axis(breaks = anchor_dt$anchor_summit, labels = with(anchor_dt, paste0(loop_id, '\n', anchor)), name = NULL)) +
      scale_y_continuous(label=scales::comma, limits = c(start(range_other), end(range_other)), oob = function(x, ...) x,
        sec.axis = dup_axis(breaks = anchor_other_dt$midpoint, labels = with(anchor_other_dt, paste0(loop_id, '\n', anchor)), name = NULL))
  else
    p <- p +
      scale_x_continuous(label=scales::comma, limits = c(start(range), end(range)), oob = function(x, ...) x) +
      scale_y_continuous(label=scales::comma, limits = c(start(range_other), end(range_other)), oob = function(x, ...) x)

  p1 <- p +
    geom_point(data = beziers, aes(x = x, y = 10e6, color = match_loop_id), alpha = 0.6) +
    geom_point(data = beziers_other, aes(x = 0, y = 12e6, color = match_loop_id), alpha = 0.6) +
    geom_bezier(data = beziers[is.na(match_loop_id), ], aes(x = x, y = y, group = loop_id, color = match_loop_id), alpha = 0.6) +
    geom_bezier(data = beziers[!is.na(match_loop_id), ], aes(x = x, y = y, group = loop_id, color = match_loop_id), alpha = 0.6)
  p2 <- p +
    geom_point(data = beziers, aes(x = 10e6, y = 0, color = match_loop_id), alpha = 0.6) +
    geom_point(data = beziers_other, aes(x = 12e6, y = y, color = match_loop_id), alpha = 0.6) +
    geom_bezier(data = beziers_other[is.na(match_loop_id), ], aes(x = x, y = y, group = loop_id, color = match_loop_id), alpha = 0.6) +
    geom_bezier(data = beziers_other[!is.na(match_loop_id), ], aes(x = x, y = y, group = loop_id, color = match_loop_id), alpha = 0.6)

  return(list(p1, p2))
}


pdf("analysis/plots/loops_Dmel_to_Dvir.pdf", width = 6, height = 6)

dt <- unique(candidates[!is.na(other_loop_id), list(loop_id, other_loop_id, A1_other_anchor, A2_other_anchor)])
for (i in seq_len(nrow(dt)))
{
  print(syntenic_plot_around_anchors(dt$loop_id[i], "A1", dt$other_loop_id[i], dt$A1_other_anchor[i]))
  print(syntenic_plot_around_anchors(dt$loop_id[i], "A2", dt$other_loop_id[i], dt$A2_other_anchor[i]))
}

dev.off()


pdf("analysis/plots/loops_Dmel_to_Dvir_mini.pdf", width = 3, height = 3)

dt <- unique(candidates[other_loop_id == "D_vir_L27", list(loop_id, other_loop_id, A1_other_anchor, A2_other_anchor)])
for (i in seq_len(nrow(dt)))
{
  print(syntenic_plot_around_anchors(dt$loop_id[i], "A1", dt$other_loop_id[i], dt$A1_other_anchor[i], margin = 100e3, anchor_margin = 0, breaks = scales::pretty_breaks(n = 2, min.n = 2)))
  print(syntenic_plot_around_anchors(dt$loop_id[i], "A2", dt$other_loop_id[i], dt$A2_other_anchor[i], margin = 100e3, anchor_margin = 0, breaks = scales::pretty_breaks(n = 2, min.n = 2)))
}

dev.off()


pdf("analysis/plots/loops_Dmel_to_Dvir_zoom.pdf", width = 6, height = 6)

dt <- unique(candidates[!is.na(other_loop_id), list(loop_id, other_loop_id, A1_other_anchor, A2_other_anchor)])
for (i in seq_len(nrow(dt)))
{
  print(syntenic_plot_around_anchors(dt$loop_id[i], "A1", dt$other_loop_id[i], dt$A1_other_anchor[i], collapse_CNEs = F, anchors_as_lines = FALSE, anchor_margin = 0, margin = 2.5e3))
  print(syntenic_plot_around_anchors(dt$loop_id[i], "A2", dt$other_loop_id[i], dt$A2_other_anchor[i], collapse_CNEs = F, anchors_as_lines = FALSE, anchor_margin = 0, margin = 2.5e3))
}

dev.off()


pdf("analysis/plots/chrom_Dmel_to_Dvir.pdf", width = 6, height = 6)

for (chrom in c("2L", "2R", "3L", "3R", "X"))
{
  print(syntenic_plot_chromosome(chrom, chrom, anchors_other_as_lines = TRUE))
  pp <- syntenic_plot_beziers(chrom, chrom)
  print(pp)

  # xy-plot to match the karyotype plot; the numbers below are chosen so that their scales match
  lim <- c(0, 39.34094e6)
  print(pp[[1]] + lims(x = lim, y = lim))
  print(pp[[2]] + lims(x = lim, y = lim))
}

dev.off()

#
#  Plot karyotypes for D. melanogaster and D. virilis assemblies
#

# Put together Seqinfo for D. melanogaster assembly

chromosomes <- c("2L", "2R", "3L", "3R", "4", "X")

seqlengths <- read.table("data/genome/D.melanogaster/dm6/ucsc/fasta/dm6.chrom.sizes", sep="\t", stringsAsFactors=FALSE)
seqlengths <- seqlengths[order(seqlengths[[1]]), ]
seqlengths[[1]] <- sub("^chr", "", seqlengths[[1]])
seqlengths <- seqlengths[seqlengths[[1]] %in% chromosomes, ]

# add a dummy chromosome of 40 Mb size
chromInfo <- Seqinfo(seqnames = c(seqlengths[[1]], "d"),
  seqlength=c(seqlengths[[2]], 40e6),
  isCircular=rep(FALSE, nrow(seqlengths) + 1),
  genome="Dmel")

cytobands <- getCytobands("dm6")
seqlevels(cytobands) <- sub("^chr", "", seqlevels(cytobands))
cytobands <- cytobands[which(seqnames(cytobands) %in% chromosomes)]

# Put together Seqinfo for D. virilis assembly

seqlengths_other <- read.table("data/genome/Drosophila_Renschler2019/GSE120751/suppl/GSE120751_Dvir_HiC.canonical.chrom.sizes", sep="\t", stringsAsFactors=FALSE)

# add a dummy chromosome of 40 Mb size
chromInfo_other <- Seqinfo(seqnames = c(seqlengths_other[[1]], "d"),
  seqlength=c(seqlengths_other[[2]], 40e6),
  isCircular=rep(FALSE, nrow(seqlengths_other) + 1),
  genome="Dvir_HiC")

# Do the plotting using karyoploteR

pdf("analysis/plots/chrom_Dmel_to_Dvir_karyotypes.pdf", width = 6, height = 6)

plot.params <- getDefaultPlotParams(plot.type=1)
plot.params$ideogramheight <- 25

kp <- plotKaryotype(genome=chromInfo, cytobands=cytobands, plot.params=plot.params)
kpAddBaseNumbers(kp, tick.dist=5e6, tick.len=20, tick.col=NA, minor.tick.dist=1e6, minor.tick.len=10, cex=0.7)
kpPlotRegions(kp, data.frame(chr="2L", start=0, end=10e6), col="#FFAACC", xlim = c(0, 0.5))
kpPoints(kp, chr=sub("^chr", "", lap$anchor_chr), x=lap$anchor_summit, y=as.integer(as.factor(lap$loop_id)) / length(unique(lap$loop_id)))

kp <- plotKaryotype(genome=chromInfo_other, plot.params=plot.params)
kpAddBaseNumbers(kp, tick.dist=5e6, tick.len=20, tick.col=NA, minor.tick.dist=1e6, minor.tick.len=10, cex=0.7)
kpPlotRegions(kp, data.frame(chr="2L", start=0, end=10e6), col="#FFAACC")
kpPoints(kp, chr=lap_other$chrom, x=lap_other$midpoint, y=as.integer(as.factor(lap_other$loop_id)) / length(unique(lap$loop_id)))

dev.off()
