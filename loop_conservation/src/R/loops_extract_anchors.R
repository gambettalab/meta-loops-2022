library(data.table)
library(GenomicRanges)
library(MapColoring)
library(RColorBrewer)

for (genome in c("D_mel", "D_vir"))
{
  #
  #  read the loops
  #

  if (genome == "D_mel")
  {
    lap <- fread(paste0("data/loops/loops_", genome, "_annotated.tsv"), header = T, sep = "\t")

    #
    #  reformat the loops to the wide format; save for use in Juicebox
    #

    dt <- dcast(lap, loop_id + aggregated_id + old_loop_id + old_aggregated_id ~ anchor,
      value.var = c("anchor_chr", "anchor_start", "anchor_end", "DHS_TSS_proximal"))
    dt[, color := "85,107,47"]
    setnames(dt, "anchor_chr_A1", "chr1")
    setnames(dt, "anchor_chr_A2", "chr2")
    setnames(dt, "anchor_start_A1", "x1")
    setnames(dt, "anchor_end_A1", "x2")
    setnames(dt, "anchor_start_A2", "y1")
    setnames(dt, "anchor_end_A2", "y2")
    setcolorder(dt, c("chr1", "x1", "x2", "chr2", "y1", "y2", "color"))

    write.table(dt, file = paste0("data/loops/long_range_loops_", genome, ".tsv"),
      quote = F, sep = "\t", row.names = F, col.names = T)
  }
  else
  {
    dt <- fread(paste0("data/loops/long_range_loops_", genome, ".tsv"), header = T, sep = "\t")
    dt[, DHS_TSS_proximal_A1 := 2]
    dt[, DHS_TSS_proximal_A2 := 2]
  }

  #
  #  calculate the adjacency matrix
  #

  maxgap_thr <- 5e3

  gr1 <- with(dt, GRanges(chr1, IRanges(x1 + 1L, x2)))
  gr2 <- with(dt, GRanges(chr2, IRanges(y1 + 1L, y2)))

  incidence <- function(gr1, gr2)
    sapply(seq_along(gr2), function(i) overlapsAny(gr1, gr2[i], maxgap = maxgap_thr))

  im <- incidence(gr1, gr1)
  im <- pmax(im, incidence(gr1, gr2))
  im <- pmax(im, incidence(gr2, gr1))
  im <- pmax(im, incidence(gr2, gr2))
  mode(im) <- "logical"

  #
  #  determine the coloring
  #

  coloring <- getColoring(im)
  N <- max(coloring)
  itemRgb = apply(col2rgb(brewer.pal(N, "Set2")), 2, paste, collapse = ",")

  #
  #  save the anchors in BED format with itemRgb column
  #

  bed_dt <- with(dt, data.table(
    chrom = c(chr1, chr2),
    start = c(x1, y1),
    end = c(x2, y2),
    name = c(paste0(loop_id, "_A1"), paste0(loop_id, "_A2")),
    score = rep(paste0(DHS_TSS_proximal_A1, DHS_TSS_proximal_A2), 2),
    strand = c(rep("+", nrow(dt)), rep("-", nrow(dt))),
    thickStart = c(x1, y1),
    thickEnd = c(x2, y2),
    rgb = rep(itemRgb[coloring], 2)
  ))
  setkey(bed_dt, chrom, start)

  fout <- paste0("data/loops/loops_", genome, "_anchors.bed")
  writeLines("track itemRgb=On", fout)
  write.table(bed_dt, file = fout, quote = F, sep = '\t', row.names = F, col.names = F, append = T)
}
