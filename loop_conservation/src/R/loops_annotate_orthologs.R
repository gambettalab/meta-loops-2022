options(warn = 1)

# Read config file
config <- yaml::yaml.load_file("config/config.yml")

library(data.table)
library(GenomicRanges)

max_TSS_distance <- 200L # for a peak to be considered TSS-proximal


#
#  Read the loop data
#

genome <- "D_mel"
other_genome <- "D_vir"

load(paste0("data/loops/loops_", genome, "_to_", other_genome, ".Rdata")) # match


#
#  Read the ortholog annotations fetched from OrthoDB
#

orthodb <- fread(paste0("data/loops/genes_loops_", genome, "_to_", other_genome, ".tsv"), header = T, sep = "\t")

# double-check that we are working on the same set of gene IDs
stopifnot(orthodb$gene_id %in% unlist(strsplit(match$lap$DHS_nearest_gene_id, ", ")))

# take only the orthologs within Metazoa (animals)
orthodb_filtered <- unique(orthodb[orthodb_level_name == "Metazoa",
  .(other_gene_id, other_gene_name, gene_n_orthologs = .N), by = .(gene_id)])


#
#  Read D. vir. gene coordinates 
#

gtf <- rtracklayer::import(config$genomes[[other_genome]]$gtf)
stopifnot(as.vector(strand(gtf)) %in% c("+", "-"))

# take only the exons of coding genes
coding_gtf <- gtf[gtf$type == "mRNA"]
coding_exons <- gtf[gtf$type == "exon" & gtf$gene_id %in% coding_gtf$gene_id]
# as we do not have transcript annotations, putative TSSes are defined as the 5' ends of all exons
tss <- resize(coding_exons, width = 1, fix = "start")


#
#  Extract the loop anchors, take the ones matched in the other species,
#  take the genes associated to the anchors, and complement these genes with orthologs 
#

# take the anchors of all loops
anchors <- match$lap[, list(
  DHS_TSS_proximal = DHS_TSS_proximal,
  DHS_distance_to_TSS = as.integer(strsplit(DHS_distance_to_TSS, ", ")[[1]]),
  gene_id = strsplit(DHS_nearest_gene_id, ", ")[[1]],
  gene_symbol = strsplit(DHS_nearest_gene_symbol, ", ")[[1]]
), by = c("loop_id", "anchor", "peak_id")]

message("Considering the ", length(unique(anchors$peak_id)), " anchors of ", length(unique(anchors$loop_id)), " loops")
dt <- unique(anchors[, c("peak_id", "DHS_TSS_proximal")])
message("Anchor split: ", sum(dt$DHS_TSS_proximal == 1), " promoter (TSS-proximal), ", sum(dt$DHS_TSS_proximal == 0), " intergenic (TSS-distal)")

# add matched anchors in the other species
anchors_matched <- merge(anchors,
  unique(match$lap_lifted[, c("loop_id", "anchor", "other_chrom", "other_start", "other_end", "other_loop_id", "other_anchor")]),
  by = c("loop_id", "anchor"), allow.cartesian = T, all.x = T, sort = F)
anchors_matched[, has_other_anchor := !is.na(other_chrom)]
anchors_matched[, other_anchor_looping := !is.na(other_loop_id)]

dt <- unique(anchors_matched[has_other_anchor == T, c("peak_id", "DHS_TSS_proximal")])
message()
message(nrow(dt), " loop anchors have a matching anchor in the other species")
message("Split of matched loop anchors: ", sum(dt$DHS_TSS_proximal == 1), " promoter (TSS-proximal), ", sum(dt$DHS_TSS_proximal == 0), " intergenic (TSS-distal)")

# add ortholog information
anchors_matched_orthologs <- merge(anchors_matched, orthodb_filtered, by = "gene_id", all.x = T, sort = F, allow.cartesian = TRUE)
# add ortholog coordinates
anchors_matched_orthologs <- merge(anchors_matched_orthologs, data.table(
  other_gene_id = tss$gene_id,
  other_gene_chrom = as.vector(seqnames(tss)),
  other_gene_tss = start(tss) - 1L,
  other_gene_strand = as.vector(strand(tss))
), by = "other_gene_id", all.x = T, sort = F)
anchors_matched_orthologs[, has_other_gene := !is.na(other_gene_id)]

dt <- unique(anchors_matched_orthologs[has_other_anchor & has_other_gene, c("peak_id", "DHS_TSS_proximal")])
message()
message(nrow(dt), " anchors have a matching anchor, and have an ortholog of anchor-associated gene(s)")
message("Split of matched loop anchors that have an ortholog: ", sum(dt$DHS_TSS_proximal == 1), " promoter (TSS-proximal), ", sum(dt$DHS_TSS_proximal == 0), " intergenic (TSS-distal)")


#
#  Do the matched anchors in the other species overlap orthologs of anchor-associated genes?
#

# calculate the distance only in cases where we have a matching anchor and an ortholog
sel <- with(anchors_matched_orthologs, has_other_anchor & !is.na(other_gene_chrom))
# in the remaining cases, store NA as the distance
anchors_matched_orthologs$other_anchor_gene_distance <- NA
anchors_matched_orthologs$other_anchor_gene_distance[sel] <- with(anchors_matched_orthologs[sel, ],
  distance(GRanges(other_gene_chrom, IRanges(other_gene_tss + 1L, other_gene_tss + 1L)),
    GRanges(other_chrom, IRanges(other_start + 1L, other_end))))

mymin <- function(x) { x <- x[!is.na(x)]; if (length(x) == 0) NA_integer_ else min(x) }
anchors_matched_orthologs_closest <- anchors_matched_orthologs[, list(other_anchor_gene_distance = mymin(other_anchor_gene_distance)),
  by = c("peak_id", "DHS_TSS_proximal", "DHS_distance_to_TSS", "gene_id", "gene_symbol", "has_other_anchor", "other_anchor_looping", "has_other_gene")]
anchors_matched_orthologs_closest <- merge(anchors_matched_orthologs_closest, anchors_matched_orthologs[is.finite(other_anchor_gene_distance), ], all.x = T, sort = F)

# sort according to the natural order of peak_id
up <- unique(anchors_matched_orthologs_closest$peak_id)
ups <- paste0("P", seq_along(up))
stopifnot(setequal(up, ups))
anchors_matched_orthologs_closest[, peak_id := factor(peak_id, ups)]
setkey(anchors_matched_orthologs_closest, peak_id)

dt <- unique(anchors_matched_orthologs_closest[other_anchor_gene_distance <= max_TSS_distance, c("peak_id", "DHS_TSS_proximal")])
message()
message(nrow(dt), " loop anchors are matching an ortholog (matching anchor in the other species overlaps an ortholog)")
message("Split of ortholog-matching loop anchors: ", sum(dt$DHS_TSS_proximal == 1), " promoter (TSS-proximal), ", sum(dt$DHS_TSS_proximal == 0), " intergenic (TSS-distal)")


#
#  Extract D. mel. anchor annotations for plotting
#

anchors_for_plot <- unique(anchors_matched_orthologs_closest[,
  c("peak_id", "DHS_TSS_proximal", "DHS_distance_to_TSS", "gene_id", "gene_symbol", "has_other_anchor", "other_chrom", "other_start", "other_end", "other_anchor_looping", "has_other_gene", "other_anchor_gene_distance", "other_gene_id", "other_gene_name")])

message()
message("Diagnostic table for resolving multi-matching cases:")
peak_ids_to_resolve <- anchors_for_plot[, list(.N), by = "peak_id"][N > 1, ]$peak_id
print(anchors_for_plot[peak_id %in% peak_ids_to_resolve, ])

# first, in D. mel. take the gene with the closest TSS to the anchor
anchors_for_plot <- anchors_for_plot[, list(
  DHS_distance_to_TSS = min(DHS_distance_to_TSS),
  gene_id = gene_id[which(DHS_distance_to_TSS == min(DHS_distance_to_TSS))],
  gene_symbol = gene_symbol[which(DHS_distance_to_TSS == min(DHS_distance_to_TSS))],
  other_anchor_gene_distance = other_anchor_gene_distance[which(DHS_distance_to_TSS == min(DHS_distance_to_TSS))],
  other_gene_id = other_gene_id[which(DHS_distance_to_TSS == min(DHS_distance_to_TSS))],
  other_gene_name = other_gene_name[which(DHS_distance_to_TSS == min(DHS_distance_to_TSS))]
), by = c("peak_id", "DHS_TSS_proximal", "has_other_anchor", "other_chrom", "other_start", "other_end", "other_anchor_looping", "has_other_gene")]

# second, in D. vir. take the candidate anchor with the smallest distance to the closest TSS
anchors_for_plot <- anchors_for_plot[, list(
  other_anchor_gene_distance = min(other_anchor_gene_distance),
  other_chrom = other_chrom[which(other_anchor_gene_distance == min(other_anchor_gene_distance))],
  other_start = other_start[which(other_anchor_gene_distance == min(other_anchor_gene_distance))],
  other_end = other_end[which(other_anchor_gene_distance == min(other_anchor_gene_distance))],
  other_anchor_looping = ifelse(length(unique(other_anchor_looping)) == 1, other_anchor_looping, other_anchor_looping[which(other_anchor_gene_distance == min(other_anchor_gene_distance))]),
  other_gene_id = other_gene_id[which(other_anchor_gene_distance == min(other_anchor_gene_distance))],
  other_gene_name = other_gene_name[which(other_anchor_gene_distance == min(other_anchor_gene_distance))]
), by = c("peak_id", "DHS_TSS_proximal", "DHS_distance_to_TSS", "gene_id", "gene_symbol", "has_other_anchor", "has_other_gene")]

# finally, concatenate everything that could not be resolved by anchor-to-closest-TSS distance
anchors_for_plot <- anchors_for_plot[, list(
  gene_id = paste(gene_id, collapse = ", "),
  gene_symbol = paste(gene_symbol, collapse = ", "),
  other_gene_id = paste(other_gene_id, collapse = ", "),
  other_gene_name = paste(other_gene_name, collapse = ", ")
), by = c("peak_id", "DHS_TSS_proximal", "DHS_distance_to_TSS", "has_other_anchor", "other_chrom", "other_start", "other_end", "other_anchor_looping", "has_other_gene", "other_anchor_gene_distance")]

# order by anchor-to-TSS distance in D. mel.
setkey(anchors_for_plot, DHS_distance_to_TSS)
anchors_for_plot[, peak_id := factor(peak_id, peak_id)]


#
#  Do the plotting
#

library(ggplot2)
library(patchwork)

pdf("analysis/plots/loops_D_mel_to_D_vir_anchor_orthologs.new.pdf", width = 4, height = 14)

p1 <- ggplot(anchors_for_plot, aes(x = DHS_distance_to_TSS / 1e3, y = peak_id)) +
  geom_segment(aes(xend = 0, yend = peak_id), color = "#e41a1c", alpha = 0.3) +
  geom_vline(xintercept = 0) +
  geom_text(aes(label = paste(gsub("-", "\u00ad", gene_symbol), " ")), hjust = 1, vjust = 0.5, size = 3, fontface = 3) +
  geom_point(color = "#e41a1c") +
  scale_x_reverse() +
  coord_cartesian(xlim = c(30, 0), clip = "off") +
  scale_y_discrete(limits = rev, position = "right") +
  xlab("Distance to TSS\nin D. mel. (kb)") +
  ylab(NULL) +
  theme_minimal() +
  NULL
# print(p1)

dt <- anchors_for_plot[, c("peak_id", "DHS_TSS_proximal", "has_other_anchor", "other_anchor_looping", "has_other_gene")]
# assert that that if the other anchor is looping, then it exists
stopifnot(with(dt, !other_anchor_looping | has_other_anchor))
dt[, has_other_anchor := 10L * as.integer(has_other_anchor) + as.integer(other_anchor_looping)]
dt[, other_anchor_looping := NULL]
dt[, has_other_gene := as.integer(has_other_gene)]
dt <- melt(dt, id.vars = "peak_id")

p2 <- ggplot(dt[value > 0, ], aes(x = variable, y = peak_id)) +
  geom_point(aes(color = as.factor(value))) +
  scale_color_manual(NULL, values = c(`1` = "#808285", `10` = "#0092CE", `11` = "#00447C"),
    labels = c(`1` = "yes", `10` = "candidate match in D. vir. is non\u00adlooping", `11` = "candidate match in D. vir. is looping"), guide = "none") +
  scale_x_discrete(labels = c("Anchor is TSS\u00adproximal in D. mel.", "Anchor has a candidate match\nin D. vir. (    looping,     non\u00adlooping)", "Gene has an ortholog in D. vir.")) +
  scale_y_discrete(limits = rev) +
  xlab(NULL) +
  ylab(NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  NULL
# print(p2)

p3 <- ggplot(anchors_for_plot, aes(x = pmin(other_anchor_gene_distance / 1e3, 30), y = peak_id)) +
  geom_segment(aes(xend = 0, yend = peak_id), color = "#e41a1c", alpha = 0.3) +
  geom_vline(xintercept = 0) +
  geom_text(aes(label = paste(" ", gsub("-", "\u00ad", gsub("Dvir\\\\", "", other_gene_name)))), hjust = 0, vjust = 0.5, size = 3, fontface = 3) +
  geom_point(data = anchors_for_plot[other_anchor_gene_distance / 1e3 <= 30, ], color = "#e41a1c") +
  coord_cartesian(xlim = c(0, 30), clip = "off") +
  scale_y_discrete(limits = rev, position = "right") +
  xlab("Distance to the\northolog TSS\nin D. vir. (kb)") +
  ylab(NULL) +
  theme_minimal() +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  NULL
# print(p3)

print(plot_spacer() + p1 + p2 + p3 + plot_spacer() + plot_layout(widths = c(1, 3, 3, 3, 1)))

dev.off()
