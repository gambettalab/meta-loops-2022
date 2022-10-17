library(data.table)
library(GenomicFeatures)

max_TSS_distance <- 200L # for a peak to be considered TSS-proximal

#
#  read the manually annotated loops
#

lap <- as.data.table(read.table("data/loops/loops_D_mel.tsv",
  header = T, sep = "\t", comment.char = "#", stringsAsFactors = F))

#
#  overlap with the previous version of manual annotations, and update old IDs
#

lap_old <- as.data.table(read.table("data/loops/loops_D_mel_old.tsv",
  header = T, sep = "\t", comment.char = "#", stringsAsFactors = F))

gr <- with(lap, GRanges(sub("^chr", "", anchor_chr), IRanges(anchor_start + 1L, anchor_end)))
gr_old <- with(lap_old, GRanges(sub("^chr", "", chrom), IRanges(DHS_start + 1L, DHS_end)))
ov <- findOverlaps(gr, gr_old)

# add old_loop_id
dt <- data.table(
  loop_id = lap$loop_id[queryHits(ov)],
  anchor = lap$anchor[queryHits(ov)],
  old_loop_id = lap_old$loop_id[subjectHits(ov)],
  old_anchor = lap_old$anchor[subjectHits(ov)]
)

dtm <- dcast(dt, loop_id + old_loop_id ~ paste0(anchor, old_anchor))[A1A1 * A2A2 > 0 | A1A2 * A2A1 > 0, ]
dtm <- dtm[, list(loop_id, old_loop_id)]
lap[, old_loop_id := NULL]
lap <- merge(lap, dtm, by = "loop_id", all.x = TRUE, sort = FALSE)

# add old_peak_id
dt <- data.table(
  loop_id = lap$loop_id[queryHits(ov)],
  peak_id = lap$peak_id[queryHits(ov)],
  anchor = lap$anchor[queryHits(ov)],
  old_loop_id = lap_old$loop_id[subjectHits(ov)],
  old_peak_id = lap_old$peak_id[subjectHits(ov)],
  old_anchor = lap_old$anchor[subjectHits(ov)]
)

dtm <- unique(dt[, list(loop_id, old_loop_id, peak_id, old_peak_id, anchor, old_anchor)])
lap[, old_peak_id := NULL]
lap[, old_anchor := NULL]
lap <- merge(lap, dtm, by = c("loop_id", "peak_id", "anchor", "old_loop_id"), all.x = TRUE, sort = FALSE)

# add old_aggregated_id
dtm <- unique(lap_old[, list(old_loop_id = loop_id, old_aggregated_id = aggregated_id)])
lap[, old_aggregated_id := NULL]
lap <- merge(lap, dtm, by = "old_loop_id", all.x = TRUE, sort = FALSE)

# add aggregated_id
dtm <- lap[, list(loop_id, aggregated_id = paste(unique(loop_id), collapse = "_")), by = "old_aggregated_id"]
dtm <- unique(dtm[, list(loop_id, aggregated_id)])
lap[, aggregated_id := NULL]
lap <- merge(lap, dtm, by = "loop_id", all.x = TRUE, sort = FALSE)

# aggregate new loop L56 together with L55 and L57
stopifnot(sum(lap$aggregated_id == "L55_L57") == 4)
lap$aggregated_id[lap$aggregated_id == "L55_L57"] <- "L55_L56_L57"
lap$aggregated_id[lap$loop_id == "L56"] <- "L55_L56_L57"

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
#  find the nearest TSS for all loop peaks
#

annotate_nearest_gene <- function(lap, one_closest_non_overlapping = TRUE)
{
  # convert chromosome names to FlyBase style (chr2L to 2L)
  gr <- GRanges(sub("^chr", "", lap$anchor_chr), IRanges(lap$anchor_start + 1L, lap$anchor_end))

  if (one_closest_non_overlapping)
    nearest_tss <- as.data.table(nearest(gr + max_TSS_distance, tss, select = "all"))
  else
  {
    # an option: take two closest TSSes instead of one for peaks not overlapping a TSS
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
    DHS_TSS_proximal = as.integer(min(distance) <= max_TSS_distance),
    DHS_distance_to_TSS = paste(unique(distance), collapse=", "),
    DHS_nearest_gene_id = paste(unique(gene_id), collapse=", "),
    DHS_nearest_gene_symbol = paste(unique(gene_symbol), collapse=", ")
  ), by = queryHits]

  stopifnot(all.equal(nearest_tss_aggr$queryHits, seq_len(nrow(lap))))
  lap$DHS_TSS_proximal <- nearest_tss_aggr$DHS_TSS_proximal
  lap$DHS_distance_to_TSS <- nearest_tss_aggr$DHS_distance_to_TSS
  lap$DHS_nearest_gene_id <- nearest_tss_aggr$DHS_nearest_gene_id
  lap$DHS_nearest_gene_symbol <- nearest_tss_aggr$DHS_nearest_gene_symbol

  # peak_aggr <- lap[, list(
  #   peak_TSS_proximal = max(DHS_TSS_proximal),
  #   peak_nearest_gene_symbol = paste(sort(unique(DHS_nearest_gene_symbol)), collapse=", ")
  # ), by = peak_id]
  # lap$peak_TSS_proximal <- with(peak_aggr, peak_TSS_proximal[match(lap$peak_id, peak_id)])
  # lap$peak_nearest_gene_symbol <- with(peak_aggr, peak_nearest_gene_symbol[match(lap$peak_id, peak_id)])

  DHS_distance_fun <- function(anchor, anchor_summit)
  {
    return(max(anchor_summit[anchor == "A2"]) - min(anchor_summit[anchor == "A1"]))
  }
  DHS_distance_aggr <- lap[, list(loop_DHS_distance = DHS_distance_fun(anchor, anchor_summit)),
    by = list(loop_id)]
  lap$loop_DHS_distance <- with(DHS_distance_aggr, loop_DHS_distance[match(lap$loop_id, loop_id)])

  return(lap)
}

lap <- annotate_nearest_gene(lap)

# # add DHS start and end coordinates collapsed per (loop, peak)
# lap_aggr <- lap[, list(
#   anchor_start_collapsed = min(anchor_start),
#   anchor_end_collapsed = max(anchor_end)
# ), by = list(loop_id, peak_id)]
# lap$anchor_start_collapsed <- with(lap_aggr, anchor_start_collapsed[
#   match(paste(lap$loop_id, lap$peak_id), paste(loop_id, peak_id))])
# lap$anchor_end_collapsed <- with(lap_aggr, anchor_end_collapsed[
#   match(paste(lap$loop_id, lap$peak_id), paste(loop_id, peak_id))])

#
#  save the loops with extra annotations added
#

write.table(lap, file = "data/loops/loops_D_mel_annotated.tsv",
  quote = F, sep = "\t", row.names = F, col.names = T)

genes_proximal <- data.frame(gene_id =
  unique(do.call(c, strsplit(lap$DHS_nearest_gene_id[lap$DHS_TSS_proximal == 1], ", "))))
write.table(genes_proximal, file = "data/loops/genes_proximal_loops_D_mel.tsv",
  quote = F, sep = "\t", row.names = F, col.names = T)

genes_distal <- data.frame(gene_id =
  unique(do.call(c, strsplit(lap$DHS_nearest_gene_id[lap$DHS_TSS_proximal == 0], ", "))))
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
