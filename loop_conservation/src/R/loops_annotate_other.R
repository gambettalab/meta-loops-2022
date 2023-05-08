options(warn = 1)

# Read config file
config <- yaml::yaml.load_file("config/config.yml")

library(data.table)
library(GenomicFeatures)

max_TSS_distance <- 200L # for an anchor to be considered TSS-proximal


#
#  Read the loop data
#

other_genome <- "D_vir"

lap_other <- as.data.table(read.table(paste0("data/loops/loops_", other_genome, ".tsv"),
  header = T, sep = "\t", quote = "", comment.char = "#", stringsAsFactors = F))


#
#  Align the loop anchor boundaries to DHS boundaries
#

# read the DHS
dhs_dt <- fread("data/Dvir_DNase-seq/dvir_25-28_IDR_0.05_in_GSE120751_assembly_cleaned.bed", header = F, sep = "\t")
names(dhs_dt)[1] <- "chrom"
dhs_dt[, start := pmin(V2, V3)] # note: 1-based inclusive, hence no +1L
dhs_dt[, end := pmax(V2, V3)]
dhs <- GRanges(dhs_dt)

# find the DHS overlapping each anchor
gr <- with(lap_other, GRanges(anchor_chr, IRanges(anchor_start + 1L, anchor_end)))
ov_dt <- as.data.table(findOverlaps(gr, dhs))
ov_dt[, start := start(dhs)[subjectHits] - 1L]
ov_dt[, end := end(dhs)[subjectHits]]
ov_dt <- ov_dt[, list(start = min(start), end = max(end), count = .N), by = "queryHits"]

# replace the start with the start of leftmost overlapping DHS
setnames(lap_other, "anchor_start", "anchor_start.MCG")
lap_other[, anchor_start := NA_integer_]
lap_other$anchor_start[ov_dt$queryHits] <- ov_dt$start

# replace the end with the end of rightmost overlapping DHS
setnames(lap_other, "anchor_end", "anchor_end.MCG")
lap_other[, anchor_end := NA_integer_]
lap_other$anchor_end[ov_dt$queryHits] <- ov_dt$end

lap_other[, anchor_DHS_count := NA_integer_]
lap_other$anchor_DHS_count[ov_dt$queryHits] <- ov_dt$count

lap_other[, anchor_midpoint := as.integer((anchor_start + anchor_end) / 2)]


#
#  Read D. vir. transcript database
#

# read gene id to symbol mapping
gtf <- rtracklayer::import(config$genomes[[other_genome]]$gtf)
gene_symbols <- unique(as.data.table(elementMetadata(gtf)[, c("gene_id", "gene_symbol")]))

# all Dvir_HiC transcripts
txdb <- loadDb("data/genome/Drosophila_Renschler2019/txdb/GSE120751_Dvir_HiC.sqlite")
tx <- transcripts(txdb, columns=c("gene_id", "tx_name", "TXTYPE"))
tx$gene_id <- sapply(tx$gene_id, paste)
tx$gene_symbol <- gene_symbols$gene_symbol[match(tx$gene_id, gene_symbols$gene_id)]

# take only the coding transcripts
coding_tx <- tx[tx$TXTYPE == "mRNA"]

# convert transcripts to TSSes
tss <- resize(coding_tx, width=1, fix='start')


#
#  Find the nearest TSS for all loop anchors
#

annotate_nearest_gene <- function(lap)
{
  gr <- with(lap, GRanges(anchor_chr, IRanges(anchor_start + 1L, anchor_end)))
  nearest_gene <- as.data.table(nearest(gr, tss, select = "all"))
  nearest_gene[, distance := distance(gr[queryHits], tss[subjectHits])]
  nearest_gene[, gene_id := tss[subjectHits]$gene_id]
  nearest_gene[, gene_symbol := tss[subjectHits]$gene_symbol]
  nearest_gene <- nearest_gene[, list(distance = min(distance)), by = list(queryHits, gene_id, gene_symbol)]

  nearest_gene_aggr <- nearest_gene[, list(
    anchor_TSS_proximal = as.integer(min(distance) <= max_TSS_distance),
    anchor_distance_to_TSS = paste(unique(distance), collapse=", "),
    anchor_nearest_gene_id = paste(unique(gene_id), collapse=", "),
    anchor_nearest_gene_symbol = paste(unique(gene_symbol), collapse=", ")
  ), by = queryHits]

  stopifnot(all.equal(nearest_gene_aggr$queryHits, seq_len(nrow(lap))))
  lap$anchor_TSS_proximal <- nearest_gene_aggr$anchor_TSS_proximal
  lap[, anchor_type := ifelse(anchor_TSS_proximal, "P", "I")]
  lap$anchor_distance_to_TSS <- nearest_gene_aggr$anchor_distance_to_TSS
  lap$anchor_nearest_gene_id <- nearest_gene_aggr$anchor_nearest_gene_id
  lap$anchor_nearest_gene_symbol <- nearest_gene_aggr$anchor_nearest_gene_symbol

  anchor_distance_fun <- function(anchor, anchor_midpoint)
  {
    return(max(anchor_midpoint[anchor == 2]) - min(anchor_midpoint[anchor == 1]))
  }
  anchor_distance_aggr <- lap[, list(loop_anchor_distance = anchor_distance_fun(anchor, anchor_midpoint)),
    by = list(loop_id)]
  lap$loop_anchor_distance <- with(anchor_distance_aggr, loop_anchor_distance[match(lap$loop_id, loop_id)])

  return(lap)
}

lap_other_annotated <- annotate_nearest_gene(lap_other)


#
#  Save the loops with extra annotations added
#

write.table(lap_other_annotated, file = paste0("data/loops/loops_", other_genome, "_annotated.tsv"),
  quote = F, sep = "\t", row.names = F, col.names = T)
