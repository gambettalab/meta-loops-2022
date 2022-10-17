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
#  Find the nearest TSS for all loop peaks
#

annotate_nearest_gene <- function(lap)
{
  gr <- with(lap, GRanges(chrom, IRanges(start + 1L, end)))
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
  lap$anchor_distance_to_TSS <- nearest_gene_aggr$anchor_distance_to_TSS
  lap$anchor_nearest_gene_id <- nearest_gene_aggr$anchor_nearest_gene_id
  lap$anchor_nearest_gene_symbol <- nearest_gene_aggr$anchor_nearest_gene_symbol

  anchor_distance_fun <- function(anchor, midpoint)
  {
    return(max(midpoint[anchor == "A2"]) - min(midpoint[anchor == "A1"]))
  }
  anchor_distance_aggr <- lap[, list(loop_anchor_distance = anchor_distance_fun(anchor, midpoint)),
    by = list(loop_id)]
  lap$loop_anchor_distance <- with(anchor_distance_aggr, loop_anchor_distance[match(lap$loop_id, loop_id)])

  return(lap)
}

lap_other_annotated <- annotate_nearest_gene(match$lap_other)


#
#  Number the peaks
#

peaks <- unique(lap_other_annotated[, c("chrom", "start", "end")])
peaks[, peak_id := paste0("D_vir_P", .I)]

# check that the peaks are non-overlapping (this is not required by design, but was the case here)
gr <- GRanges(peaks)
stopifnot(length(findOverlaps(gr, drop.self = T)) == 0)

lap_other_annotated <- merge(lap_other_annotated, peaks, by = c("chrom", "start", "end"), all.x = T, sort = F)


#
#  Save the loops with extra annotations added
#

write.table(lap_other_annotated, file = paste0("data/loops/loops_", other_genome, "_annotated.tsv"),
  quote = F, sep = "\t", row.names = F, col.names = T)
