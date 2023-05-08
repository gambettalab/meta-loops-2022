options(warn = 1)

library(data.table)
library(GenomicFeatures)

max_TSS_distance <- 200L # for an anchor to be considered TSS-proximal

#
#  read the manually annotated loops
#

lap <- as.data.table(read.table("data/loops/loops_D_mel.tsv",
  header = T, sep = "\t", comment.char = "#", stringsAsFactors = F))

setnames(lap, "anchor_nearest_gene_symbol", "old_anchor_nearest_gene_symbol", skip_absent = TRUE)
setnames(lap, "anchor_nearest_gene_symbol.1", "old_anchor_nearest_gene_symbol.1", skip_absent = TRUE)
setnames(lap, "anchor_distance_to_TSS", "old_anchor_distance_to_TSS", skip_absent = TRUE)

#
#  ensure that the loop anchor boundaries are at DHS boundaries
#

# read the DHS, and take only the 10-12h Elav peaks
dhs_dt <- fread("data/DNase-seq_Reddington2020/DHS_clusters_dm6.tsv", header = T, sep = "\t")
dhs <- with(dhs_dt[`10-12h_Neuro_peak` > 0, ], GRanges(chr, IRanges(start + 1L, end)))
# convert chromosome names to FlyBase style (chr2L to 2L)
seqlevels(dhs) <- sub("^chr", "", seqlevels(dhs))

# find the DHS overlapping each anchor
gr <- with(lap, GRanges(anchor_chr, IRanges(anchor_start + 1L, anchor_end)))
ov_dt <- as.data.table(findOverlaps(gr, dhs))
ov_dt[, start := start(dhs)[subjectHits] - 1L]
ov_dt[, end := end(dhs)[subjectHits]]
ov_dt <- ov_dt[, list(start = min(start), end = max(end), count = .N), by = "queryHits"]

lap[, DHS_start := NA_integer_]
lap$DHS_start[ov_dt$queryHits] <- ov_dt$start
stopifnot(with(lap, all.equal(anchor_start, DHS_start)))
lap[, DHS_start := NULL]

lap[, DHS_end := NA_integer_]
lap$DHS_end[ov_dt$queryHits] <- ov_dt$end
stopifnot(with(lap, all.equal(anchor_end, DHS_end)))
lap[, DHS_end := NULL]

lap[, anchor_DHS_count := NA_integer_]
lap$anchor_DHS_count[ov_dt$queryHits] <- ov_dt$count

#
#  read transcript database
#

# read gene id to symbol mapping
gtf <- rtracklayer::import("data/genome/D.melanogaster/dm6/dmel_r6.36_FB2020_05/gtf/dmel-all-r6.36.gtf.gz")
gene_symbols <- unique(as.data.table(elementMetadata(gtf)[, c("gene_id", "gene_symbol")]))

# all FlyBase transcripts
txdb <- loadDb("data/genome/D.melanogaster/dm6/dmel_r6.36_FB2020_05/txdb/dmel-all-r6.36.sqlite")
tx <- transcripts(txdb, columns=c("gene_id", "tx_name", "TXTYPE"))
tx$gene_id <- sapply(tx$gene_id, paste)
tx$gene_symbol <- gene_symbols$gene_symbol[match(tx$gene_id, gene_symbols$gene_id)]

# take only the coding transcripts
coding_tx <- tx[tx$TXTYPE == "mRNA"]

# convert transcripts to TSSes
tss <- resize(coding_tx, width=1, fix='start')

#
#  find the nearest TSS for all loop anchors
#

annotate_nearest_gene <- function(lap, one_closest_non_overlapping = TRUE)
{
  gr <- with(lap, GRanges(anchor_chr, IRanges(anchor_start + 1L, anchor_end)))

  if (one_closest_non_overlapping)
    nearest_tss <- as.data.table(nearest(gr + max_TSS_distance, tss, select = "all"))
  else
  {
    # an option: take two closest TSSes instead of one for anchors not overlapping a TSS
    ov <- findOverlaps(gr + max_TSS_distance, tss)
    queryHits_not_overlaping <- setdiff(seq_along(gr), unique(queryHits(ov)))
    nearest_tss <- rbind(
      as.data.table(ov),
      data.table(
        queryHits = queryHits_not_overlaping,
        subjectHits = precede(gr[queryHits_not_overlaping], tss)
      ),
      data.table(
        queryHits = queryHits_not_overlaping,
        subjectHits = follow(gr[queryHits_not_overlaping], tss)
      )
    )
    setkey(nearest_tss, queryHits)
  }

  nearest_tss[, distance := distance(gr[queryHits], tss[subjectHits])]
  nearest_tss[, gene_id := tss[subjectHits]$gene_id]
  nearest_tss[, gene_symbol := tss[subjectHits]$gene_symbol]
  nearest_tss <- nearest_tss[, list(distance = min(distance)), by = list(queryHits, gene_id, gene_symbol)]

  nearest_tss_aggr <- nearest_tss[, list(
    anchor_TSS_proximal = as.integer(min(distance) <= max_TSS_distance),
    anchor_distance_to_TSS = paste(unique(distance), collapse=", "),
    anchor_nearest_gene_id = paste(unique(gene_id), collapse=", "),
    anchor_nearest_gene_symbol = paste(unique(gene_symbol), collapse=", ")
  ), by = queryHits]

  stopifnot(all.equal(nearest_tss_aggr$queryHits, seq_len(nrow(lap))))
  lap$anchor_TSS_proximal <- nearest_tss_aggr$anchor_TSS_proximal
  lap[, anchor_type := ifelse(anchor_TSS_proximal, "P", "I")]
  lap$anchor_distance_to_TSS <- nearest_tss_aggr$anchor_distance_to_TSS
  lap$anchor_nearest_gene_id <- nearest_tss_aggr$anchor_nearest_gene_id
  lap$anchor_nearest_gene_symbol <- nearest_tss_aggr$anchor_nearest_gene_symbol

  anchor_distance_fun <- function(anchor, anchor_summit)
  {
    return(max(anchor_summit[anchor == 2]) - min(anchor_summit[anchor == 1]))
  }
  anchor_distance_aggr <- lap[, list(loop_anchor_distance = anchor_distance_fun(anchor, anchor_summit)),
    by = list(loop_id)]
  lap$loop_anchor_distance <- with(anchor_distance_aggr, loop_anchor_distance[match(lap$loop_id, loop_id)])

  return(lap)
}

lap <- annotate_nearest_gene(lap)

# # add DHS start and end coordinates collapsed per (loop, anchor)
# lap_aggr <- lap[, list(
#   anchor_start_collapsed = min(anchor_start),
#   anchor_end_collapsed = max(anchor_end)
# ), by = list(loop_id, anchor_id)]
# lap$anchor_start_collapsed <- with(lap_aggr, anchor_start_collapsed[
#   match(paste(lap$loop_id, lap$anchor_id), paste(loop_id, anchor_id))])
# lap$anchor_end_collapsed <- with(lap_aggr, anchor_end_collapsed[
#   match(paste(lap$loop_id, lap$anchor_id), paste(loop_id, anchor_id))])

#
#  save the loops with extra annotations added
#

write.table(lap, file = "data/loops/loops_D_mel_annotated.tsv",
  quote = F, sep = "\t", row.names = F, col.names = T)

genes_proximal <- data.frame(gene_id =
  unique(do.call(c, strsplit(lap$anchor_nearest_gene_id[lap$anchor_TSS_proximal == 1], ", "))))
write.table(genes_proximal, file = "data/loops/genes_proximal_loops_D_mel.tsv",
  quote = F, sep = "\t", row.names = F, col.names = T)

genes_distal <- data.frame(gene_id =
  unique(do.call(c, strsplit(lap$anchor_nearest_gene_id[lap$anchor_TSS_proximal == 0], ", "))))
write.table(genes_distal, file = "data/loops/genes_distal_loops_D_mel.tsv",
  quote = F, sep = "\t", row.names = F, col.names = T)

#
#  save the gene universe for completeness
#

genes_universe <- data.frame(gene_id = unique(coding_tx$gene_id))
write.table(genes_universe, file = "data/loops/genes_universe_D_mel.tsv",
  quote = F, sep = "\t", row.names = F, col.names = T)

coding_tx_anno <- coding_tx
names(coding_tx_anno) <- with(coding_tx_anno, paste(gene_id, tx_name, sep = "_"))
rtracklayer::export.bed(coding_tx_anno, "data/loops/genes_universe_D_mel.bed")

#
#  statistics for the manuscript
#

message(nrow(genes_proximal), " genes attributed to promoter-proximal loop anchors (found within +/- ", max_TSS_distance, " bp of a promoter anchor).")
message(nrow(genes_distal), " genes attributed as the closest to promoter-distal loop anchors")
message(nrow(genes_universe), " genes in the universe")

#
#  save the loops with extra annotations, taking two closest TSSes instead of one for anchors not overlapping a TSS
#

lap_2TSS <- annotate_nearest_gene(lap, one_closest_non_overlapping = FALSE)

write.table(lap_2TSS, file = "data/loops/loops_D_mel_annotated_2TSS.tsv",
  quote = F, sep = "\t", row.names = F, col.names = T)

genes_proximal <- data.frame(gene_id =
  unique(do.call(c, strsplit(lap_2TSS$anchor_nearest_gene_id[lap_2TSS$anchor_TSS_proximal == 1], ", "))))
write.table(genes_proximal, file = "data/loops/genes_proximal_loops_D_mel_2TSS.tsv",
  quote = F, sep = "\t", row.names = F, col.names = T)

genes_distal <- data.frame(gene_id =
  unique(do.call(c, strsplit(lap_2TSS$anchor_nearest_gene_id[lap_2TSS$anchor_TSS_proximal == 0], ", "))))
write.table(genes_distal, file = "data/loops/genes_distal_loops_D_mel_2TSS.tsv",
  quote = F, sep = "\t", row.names = F, col.names = T)
