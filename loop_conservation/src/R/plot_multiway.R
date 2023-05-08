# Read config file
config <- yaml::yaml.load_file("config/config.yml")

library(data.table)
library(GenomicRanges)
library(ggplot2)
require(ggsignif)

datasets <- c("HiC_WT_L3_CNS_Rep1", "HiC_WT_L3_CNS_Rep2", "HiC_WT_L3_CNS_Rep4")
bin_size <- 2000L # for ploting 1D interaction profiles
options(mc.cores = 32)


#
#  Read the loop data
#

genome <- "D_mel"
other_genome <- "D_vir"

load(paste0("data/loops/loops_", genome, "_to_", other_genome, ".Rdata")) # match


#
#  Extract the anchors
#

anchor_dt <- unique(match$lap[, c("anchor_id", "anchor_chr", "anchor_start", "anchor_end", "anchor_summit")])
setkey(anchor_dt, anchor_chr, anchor_start)
anchor_gr <- with(anchor_dt, GRanges(anchor_chr, IRanges(anchor_start + 1L, anchor_end)))

# check that the anchors are non-overlapping (this is not required by design, but was the case here)
stopifnot(length(findOverlaps(anchor_gr, drop.self = T)) == 0)

# check that the anchors are sorted by their genomic coordinates (again, not required by design, but the case here)
stopifnot(identical(anchor_gr, sort(anchor_gr)))


#
#  Read the multi-way interactions for all the anchors, as well as Hi-C normalization factors
#

mt_dt <- list()
mt_gr <- list()

for (anchor_id in anchor_dt$anchor_id)
{
  dt <- NULL
  for (dataset in datasets)
    dt <- rbind(dt, fread(paste0("data/HiC/multiway/", config$genome$symbol, "_", dataset, "_", anchor_id, ".tsv.gz"), header = T, sep = "\t"))

  gr <- with(dt, GRanges(sub("^chr", "", chrom), IRanges(start + 1L, end)))
  gr$read_id <- dt$read_id

  mt_dt[[anchor_id]] <- dt
  mt_gr[[anchor_id]] <- gr
}


#
#  Plot some summaries
#

pdf("analysis/plots/multiway_stats_alignments_read.pdf", width = 3.5, height = 2)
dt <- mt_dt[[1]][, list(.N), by = "read_id"]
p <- ggplot(dt, aes(x = N)) +
  geom_bar() +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Number of alignments per paired\u00adend read", y = "Count") +
  NULL
print(p)
dev.off()


#
#  Functions to filter multi-way interactions, conditioning on the given anchor(s)
#

filter_interactions <- function(mt_gr, anchors, margin)
{
  this_mt_gr <- mt_gr

  for (i in seq_along(anchors))
  {
    this_mt_dt <- as.data.table(this_mt_gr)
    this_mt_dt$overlaps <- overlapsAny(this_mt_gr + margin, anchors[i])
    this_mt_dt[, overlaps := any(overlaps), by = "read_id"]
    this_mt_gr <- this_mt_gr[this_mt_dt$overlaps]
  }

  this_mt_gr$overlaps <- NULL
  return(this_mt_gr)
}

get_anchor <- function(anchor_id)
{
  return(anchor_gr[match(anchor_id, anchor_dt$anchor_id)])
}

count_interactions <- function(mt_gr)
{
  return(length(unique(mt_gr$read_id)))
}


#
#  Function to calculate interaction profiles, conditioning on the given anchor(s)
#

interaction_profile <- function(region, mt_gr, normalize = TRUE)
{
  vf <- floor(start(region) / bin_size):ceiling(end(region) / bin_size)
  int_fill <- data.table(chrom = as.character(seqnames(region)),
    start = vf * bin_size + 1L, end = vf * bin_size + bin_size, value = 0L)

  int <- as.data.table(resize(mt_gr, 1, fix = "center"))
  setnames(int, "seqnames", "chrom")
  norm_f <- count_interactions(int)
  int[, start := as.integer(floor(start / bin_size)) * bin_size + 1L]
  int[, end := as.integer(ceiling(start / bin_size)) * bin_size]
  int <- int[chrom == as.character(seqnames(region)) & start(region) <= end & start <= end(region), ]
  int <- int[, list(value = length(unique(read_id)) * ifelse(normalize, 1 / norm_f, 1L)),
    by = c("chrom", "start", "end")]

  int <- rbind(int[, c("chrom", "start", "end", "value")], int_fill)
  int <- int[, list(value = sum(value)), by = c("chrom", "start", "end")]
  setkey(int, "chrom", "start")
  int[, pos := start - 1L + as.integer(bin_size / 2)]

  return(int)
}


#
#  Find all possible and all observed three-way interactions that use the observed anchors
#

build_triplets <- function(v)
{
  n <- length(v)
  dt <- NULL
  if (n > 2)
    for (i in seq(1, n - 2))
      for (j in seq(i + 1, n - 1))
        for (k in seq(j + 1, n))
          dt <- rbind(dt, data.table(p1 = v[i], p2 = v[j], p3 = v[k]))
  stopifnot(nrow(dt) == choose(n, 3))
  return(dt)
}

# first put together all the possible triplets (within each chromosome separately)
threeway_list <- lapply(with(anchor_dt, split(anchor_id, anchor_chr)), build_triplets)
threeway_dt <- do.call(rbind, threeway_list)

# annotate the triplets with their distances
threeway_dt[, distance12 := with(anchor_dt, anchor_summit[match(p2, anchor_id)] - anchor_summit[match(p1, anchor_id)])]
threeway_dt[, distance23 := with(anchor_dt, anchor_summit[match(p3, anchor_id)] - anchor_summit[match(p2, anchor_id)])]

# now mark all the observed loops
lap1 <- match$lap[anchor == 1, ]
lap2 <- match$lap[anchor == 2, ]
stopifnot(all.equal(lap1$loop_id, lap2$loop_id))
lap_anchor_id <- paste(lap1$anchor_id, lap2$anchor_id)
threeway_dt[, loop12 := lap1$loop_id[match(paste(p1, p2), lap_anchor_id)]]
threeway_dt[, loop13 := lap1$loop_id[match(paste(p1, p3), lap_anchor_id)]]
threeway_dt[, loop23 := lap1$loop_id[match(paste(p2, p3), lap_anchor_id)]]

# limit to the triplets where at least one pair of anchors matches an observed loop
threeway_dt <- threeway_dt[!is.na(loop12) | !is.na(loop13) | !is.na(loop23), ]

# explicitly mark the observed triplets
sel <- with(threeway_dt, !is.na(loop12) & !is.na(loop13) & !is.na(loop23))
threeway_dt[, group := factor(ifelse(sel, "observed", "control"), c("observed", "control"))]
threeway_dt[, id := ifelse(sel, paste(loop12, loop13, loop23, sep = "_"), NA)]


#
#  Build another, "mirror" set of three-way interactions to use as a control
#

mirror_margin <- 100e3

for (i in which(threeway_dt$group == "observed"))
{
    p1 <- threeway_dt$p1[i]
    p2 <- threeway_dt$p2[i]
    p3 <- threeway_dt$p3[i]
    gr <- c(get_anchor(p1), get_anchor(p2), get_anchor(p3))

    ip1 <- match(p1, anchor_dt$anchor_id)
    ip2 <- match(p2, anchor_dt$anchor_id)
    ip3 <- match(p3, anchor_dt$anchor_id)

    # control shifted left
    ca <- data.table(anchor_id = paste0("C", i, "A"), anchor_chr = anchor_dt$anchor_chr[ip1],
      anchor_start = anchor_dt$anchor_start[ip1] - (anchor_dt$anchor_end[ip3] - anchor_dt$anchor_end[ip2]) - mirror_margin,
      anchor_end = anchor_dt$anchor_start[ip1] - (anchor_dt$anchor_start[ip3] - anchor_dt$anchor_end[ip2]) + mirror_margin,
      anchor_summit = anchor_dt$anchor_start[ip1] - (anchor_dt$anchor_summit[ip3] - anchor_dt$anchor_end[ip2]))
    # control shifted right
    cb <- data.table(anchor_id = paste0("C", i, "B"), anchor_chr = anchor_dt$anchor_chr[ip3],
      anchor_start = anchor_dt$anchor_end[ip3] + (anchor_dt$anchor_start[ip2] - anchor_dt$anchor_end[ip1]) - mirror_margin,
      anchor_end = anchor_dt$anchor_end[ip3] + (anchor_dt$anchor_start[ip2] - anchor_dt$anchor_start[ip1]) + mirror_margin,
      anchor_summit = anchor_dt$anchor_end[ip3] + (anchor_dt$anchor_start[ip2] - anchor_dt$anchor_summit[ip1]))
    anchor_dt <- rbind(anchor_dt, ca, cb)

    threeway_dt <- rbind(threeway_dt,
      data.table(p1 = paste0("C", i, "A"), p2 = p1, p3 = p2,
        distance12 = NA, distance23 = NA, loop12 = NA, loop13 = NA, loop23 = threeway_dt$loop12[i],
        group = "control_mirror", id = threeway_dt$id[i]),
      data.table(p1 = threeway_dt$p2[i], p2 = p3, p3 = paste0("C", i, "B"),
        distance12 = NA, distance23 = NA, loop12 = threeway_dt$loop23[i], loop13 = NA, loop23 = NA,
        group = "control_mirror", id = threeway_dt$id[i])
      )
  }


# re-annotate the triplets with their distances
threeway_dt[, distance12 := with(anchor_dt, anchor_summit[match(p2, anchor_id)] - anchor_summit[match(p1, anchor_id)])]
threeway_dt[, distance23 := with(anchor_dt, anchor_summit[match(p3, anchor_id)] - anchor_summit[match(p2, anchor_id)])]

# re-sort the anchors
setkey(anchor_dt, anchor_chr, anchor_start)
anchor_gr <- with(anchor_dt, GRanges(sub("^chr", "", anchor_chr), IRanges(anchor_start + 1L, anchor_end)))

# check that the anchors are sorted by their genomic coordinates (again, not required by design, but the case here)
stopifnot(identical(anchor_gr, sort(anchor_gr)))


#
#  Three-way interaction plotting in 1D
#

plot_1d_interaction_plot <- function(region, int, viewpoints, highlights, margin, ymax = 1)
{
  caption <- paste0(as.character(seqnames(region)), ":",
    format(as.integer(start(region)), big.mark = ","), "-",
    format(as.integer(end(region)), big.mark = ","))

  p <- ggplot(int, aes(pos, value)) +
    labs(x = caption, y = "Interaction frequency") +
    scale_x_continuous(labels = scales::comma) +
    scale_y_continuous(labels = scales::comma) +
    geom_rect(data = as.data.table(highlights + margin), aes(xmin = start - 1L, xmax = end), ymin = -Inf, ymax = Inf, inherit.aes = F, fill = "#386cb0", alpha = 0.4) +
    geom_point(data = as.data.table(viewpoints), aes(x = (start - 1L + end) / 2), y = Inf, inherit.aes = F, fill = "#386cb0", color = "#386cb0", shape = 25, size = 5) +
    coord_cartesian(xlim = c(start(region) - 1L, end(region)), ylim = c(0, ymax)) +
    geom_line() +
    NULL
  return(p)
}

print_all_1d_interaction_plots <- function(region, mt_gr, gr, margin, expand = 0.05)
{
  expand_width <- expand * (end(region) - start(region) + 1L)
  int23 <- interaction_profile(region + expand_width, filter_interactions(mt_gr, gr[2:3],
    margin = margin))
  int2 <- interaction_profile(region + expand_width, filter_interactions(mt_gr, gr[2],
    margin = margin))
  int12 <- interaction_profile(region + expand_width, filter_interactions(mt_gr, gr[1:2],
    margin = margin))
  # ymax <- max(int23$value, int2$value, int12$value) * expand
  ymax <- 0.05

  print(plot_1d_interaction_plot(region, int23, gr[2:3], gr, margin, ymax))
  print(plot_1d_interaction_plot(region, int2, gr[2], gr, margin, ymax))
  print(plot_1d_interaction_plot(region, int12, gr[1:2], gr, margin, ymax))
}

threeway_observed_dt <- threeway_dt[group == "observed", ]

# for (margin in c(1000L, 2000L, 5000L, 10000L))
for (margin in c(2000L))
  for (i in seq_len(nrow(threeway_observed_dt)))
  {
    message(threeway_observed_dt$id[i])
    pdf(paste0("analysis/plots/multiway_interactions_", threeway_observed_dt$id[i], "_margin_", margin / 1e3, "kb.pdf"), width = 8, height = 2)

    p1 <- threeway_observed_dt$p1[i]
    p2 <- threeway_observed_dt$p2[i]
    p3 <- threeway_observed_dt$p3[i]
    gr <- c(get_anchor(p1), get_anchor(p2), get_anchor(p3))

    region <- GRanges(seqnames(gr[1]), IRanges(start(gr[1]), end(gr[3]))) + 0.2 * (end(gr[3]) - start(gr[1]))
    start(region) <- max(start(region), 1)
    print_all_1d_interaction_plots(region, mt_gr[[p2]], gr, margin = margin)
    dev.off()
  }


#
#  Calculate enrichment statistics for possible and all observed three-way interactions
#

annotate_threeway <- function(i, margin)
{
  cat(".")
  dt <- threeway_dt[i, ]
  p1 <- dt$p1
  p2 <- dt$p2
  p3 <- dt$p3
  gr <- c(get_anchor(p1), get_anchor(p2), get_anchor(p3))

  n12 <- count_interactions(filter_interactions(mt_gr[[p2]], gr[1:2], margin = margin))
  n2 <- count_interactions(filter_interactions(mt_gr[[p2]], gr[2], margin = margin))
  n23 <- count_interactions(filter_interactions(mt_gr[[p2]], gr[2:3], margin = margin))
  n123 <- count_interactions(filter_interactions(mt_gr[[p2]], gr, margin = margin))

  dt$n12 <- n12
  dt$n2 <- n2
  dt$n23 <- n23
  dt$n123 <- n123

  dt$fold_change <- (n123 / n23) / (n12 / n2)
  cm <- rbind(c(n123, n12 - n123), c(n23 - n123, n2 - n12 - n23 + n123))
  dt$pval <- fisher.test(cm)$p.value
  return(dt)
}


#
#  Volcano plot for the observed and comparable control triplets
#

# for (margin in c(1000L, 2000L, 5000L))
for (margin in c(2000L))
{
  threeway_anno_list <- mclapply(seq_len(nrow(threeway_dt)), annotate_threeway, margin = margin)
  message("done")
  threeway_anno_dt <- do.call(rbind, threeway_anno_list)

  # ensure that subsequent anchors in a control triplet are at least 1 kb away after extending
  # by +/- margin) (the 1 kb is added so that a read overlapping two anchors will not contribute to
  # the counts) and also ensure that the distances are no more than 10 Mb away, and n123 >= 10
  dt <- threeway_anno_dt[group == "observed" | (n123 >= 10 & pmin(distance12, distance23) >= 2 * margin + 1000L & pmax(distance12, distance23) <= 10e6), ]

  dt[, padj := p.adjust(pval, method = "BH")]
  print(dt[pval < 0.01, ])

  print(threeway_anno_dt[p1 %in% c("P67", "P68") & p3 %in% c("P69", "P70"), ])
  print(threeway_anno_dt[p1 %in% c("P16", "P17") & p2 %in% c("P17", "P21") & p3 %in% c("P21", "P24"), ])
  print(threeway_anno_dt[p1 %in% c("P47", "P49") & p2 %in% c("P49", "P52") & p3 %in% c("P52", "P53"), ])

  pdf(paste0("analysis/plots/multiway_triplet_volcano_margin_", margin / 1e3, "kb.pdf"), width = 4, height = 3)

  p <- ggplot(dt, aes(log2(fold_change), -log10(pval), color = group)) +
    labs(x = expression(paste("Triplet cooperativity score (", log[2], ")")), y = "Triplet -log10 p-value") +
    # scale_x_continuous(labels = scales::comma) +
    # scale_y_continuous(labels = scales::comma) +
    geom_point(alpha = 0.6) +
    coord_cartesian(xlim = c(-1, 1) * max(abs(log2(dt$fold_change)))) +
    theme(legend.position = "bottom") +
    NULL
  print(p)

  dev.off()


  pdf(paste0("analysis/plots/multiway_triplet_stats_mirrorcontrol_ext100kb_margin_", margin / 1e3, "kb.pdf"), width = 2.25, height = 3.25)

  dt <- threeway_anno_dt[group != "control" & is.finite(log2(fold_change)), ]

  p <- ggplot(dt, aes(group, log2(fold_change))) +
    labs(x = NULL, y = expression(paste("Triplet cooperativity score (", log[2], ")"))) +
    geom_violin() +
    geom_boxplot(width = 0.2, outlier.shape = NA) +
    geom_jitter(position = position_jitter(seed = 4242), alpha = 0.4) +
    geom_signif(test = "wilcox.test", comparisons = list(c(1, 2)), margin_top = 0.1) +
    scale_x_discrete(labels = c(`observed` = "Observed\nlooping\ntriplets",
      `control_mirror` = "Control\ntriplets")) +
    NULL
  print(p)

  dev.off()
}


#
#  Export .bedGraph tracks for a 3-way interaction
#

anchors <- c("A16", "A22", "A25")
margin <- 2000L

out_dir <- "analysis/multiway_tracks"
if (!dir.exists(out_dir))
  dir.create(out_dir)

chrom_sizes <- fread(config$genome$chrom_sizes, header = F, col.names = c("chrom", "length"))

for (i in seq(1, length(anchors)))
  for (j in seq(i, length(anchors)))
  {
    p1 <- anchors[i]
    p2 <- anchors[j]
    this_anchors <- anchor_gr[anchor_dt$anchor_id %in% c(p1, p2)]
    stopifnot(length(unique(seqnames(this_anchors))) == 1)
    this_chrom <- seqnames(this_anchors)[1]
    this_chrom_length <- chrom_sizes$length[match(paste0("chr", this_chrom), chrom_sizes$chrom)]
    region <- GRanges(this_chrom, IRanges(1, this_chrom_length))

    dt <- interaction_profile(region, filter_interactions(mt_gr[[p1]], this_anchors, margin = margin))
    dt <- dt[value != 0, c("chrom", "start", "end", "value")]
    dt[, start := start - 1L]
    fout <- paste0("analysis/multiway_tracks/HiC_WT_L3_CNS_interactions_", p1,
      if (p1 == p2) "" else paste0("_", p2),
      "_margin_", margin / 1e3, "kb_bin_size_", bin_size / 1e3, "kb_normalized.bedGraph")
    write.table(dt, file = fout, quote = F, sep = '\t', row.names = F, col.names = F, append = T)
  }


#
#  Export .bedGraph tracks for a 4-way interaction
#

anchors <- c("A67", "A68", "A69", "A70")
margin <- 2000L

for (i in seq(1, length(anchors) - 1))
  for (j in seq(i + 1, length(anchors)))
  {
    p1 <- anchors[i]
    p2 <- anchors[j]
    this_anchors <- anchor_gr[anchor_dt$anchor_id %in% c(p1, p2)]
    stopifnot(length(unique(seqnames(this_anchors))) == 1)
    this_chrom <- seqnames(this_anchors)[1]
    this_chrom_length <- chrom_sizes$length[match(paste0("chr", this_chrom), chrom_sizes$chrom)]
    region <- GRanges(this_chrom, IRanges(1, this_chrom_length))

    dt <- interaction_profile(region, filter_interactions(mt_gr[[p1]], this_anchors, margin = margin),
      normalize = F)
    dt <- dt[value != 0, c("chrom", "start", "end", "value")]
    dt[, start := start - 1L]
    fout <- paste0("analysis/multiway_tracks/HiC_WT_L3_CNS_interactions_", p1,
      if (p1 == p2) "" else paste0("_", p2),
      "_margin_", margin / 1e3, "kb_bin_size_", bin_size / 1e3, "kb_normalized.bedGraph")
    write.table(dt, file = fout, quote = F, sep = '\t', row.names = F, col.names = F, append = T)
  }
