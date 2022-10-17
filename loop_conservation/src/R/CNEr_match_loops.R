options(warn = 1)

# Read config file
config <- yaml::yaml.load_file("config/config.yml")

params <- list(D_vir = "40_50")

# take CNEs closer that < segment_margin from the anchors
segment_margin <- 1e3

library(CNEr)
library(data.table)
library(GenomicRanges)
library(ggplot2)

#
#  Extract the loops in both species
#

extract_loops <- function(genome)
{
  if (genome == "D_mel")
    loops <- fread(paste0("data/loops/loops_", genome, "_annotated.tsv"), header = T, sep = "\t")
  else
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
  }

  return(loops)
}

#
#  Lift the anchors to the other species
#

lift_anchors <- function(lap, cne)
{
  # prepare the chain of conserved elements (CEs) for overlapping
  cne_dt <- as.data.table(cne)
  cne_dt[, first.midpoint := (first.start + first.end) / 2]
  cne_dt[, second.midpoint := (second.start + second.end) / 2]

  # determine the chain orientation
  cne_dt[, chain.first.ori := unique(head(first.midpoint, -1) < tail(first.midpoint, -1)), by = "chain_id"]
  cne_dt[, chain.second.ori := unique(head(second.midpoint, -1) < tail(second.midpoint, -1)), by = "chain_id"]

  cne_dt[, chain.first.start := ifelse(chain.first.ori, first.end, first.start)]
  cne_dt[, chain.first.end := ifelse(chain.first.ori, c(tail(first.start, -1), NA), c(tail(first.end, -1), NA)), by = "chain_id"]

  cne_dt[, chain.second.start := ifelse(chain.second.ori, second.end, second.start)]
  cne_dt[, chain.second.end := ifelse(chain.second.ori, c(tail(second.start, -1), NA), c(tail(second.end, -1), NA)), by = "chain_id"]

  # concatenate the ranges of CEs and the ranges of intervening chain segments
  concatenated_chain_gr <- c(CNEr::first(cne),
    with(cne_dt[!is.na(chaining_distance), ], GRanges(first.seqnames, IRanges(pmin(chain.first.start, chain.first.end), pmax(chain.first.start, chain.first.end)))))
  concatenated_other_chain_gr <- c(CNEr::second(cne),
    with(cne_dt[!is.na(chaining_distance), ], GRanges(second.seqnames, IRanges(pmin(chain.second.start, chain.second.end), pmax(chain.second.start, chain.second.end)))))
  concatenated_chain_id <- c(elementMetadata(cne)$chain_id, cne_dt[!is.na(chaining_distance), ]$chain_id)

  # do the overlap with source genome anchors
  if ("anchor_chr" %in% names(lap))
    lap_gr <- with(lap, GRanges(sub("^chr", "", anchor_chr), IRanges(anchor_start + 1L, anchor_end)))
  else
    lap_gr <- with(lap, GRanges(chrom, IRanges(start + 1L, end)))
  ov <- findOverlaps(lap_gr + segment_margin, concatenated_chain_gr)

  # data table for anchors lifted to the other species
  lap_lifted <- data.table()
  for (qh in unique(queryHits(ov)))
  {
    ov_subset <- ov[queryHits(ov) == qh]

    gr <- concatenated_chain_gr[subjectHits(ov_subset)]
    other_gr <- concatenated_other_chain_gr[subjectHits(ov_subset)]
    chain_id <- concatenated_chain_id[subjectHits(ov_subset)]

    lap_lifted <- rbind(lap_lifted, data.table(
      loop_id = lap$loop_id[qh],
      anchor = lap$anchor[qh],
      chrom = as.character(seqnames(gr)),
      start = start(gr) - 1L,
      end = end(gr),
      chain_id = chain_id,
      other_chrom = as.character(seqnames(other_gr)),
      other_start = start(other_gr) - 1L,
      other_end = end(other_gr)
    ))
  }

  # join the lifted anchor matches for the same alignment chain
  lap_lifted <- lap_lifted[, list(
    chrom = chrom[1],
    start = min(start),
    end = max(end),
    other_chrom = other_chrom[1],
    other_start = min(other_start),
    other_end = max(other_end)
  ), by = c("loop_id", "anchor", "chain_id")]

  # prepare a helper data table with CE chain ranges, to use for annotation
  chain_dt <- cne_dt[, list(
    chain_range = paste0(unique(first.seqnames), ":", min(first.start), "-", max(first.end)),
    other_chain_range = paste0(unique(second.seqnames), ":", min(second.start), "-", max(second.end)),
    other_chain_strand = paste(unique(second.strand), collapse = ",")
  ), by = "chain_id"]

  lap_lifted <- merge(lap_lifted, chain_dt, by = "chain_id", all.x = TRUE, sort = FALSE)

  return(lap_lifted)
}

#
#  Overlap the lifted anchors with the actual anchors in the other species
#

overlap_lifted_anchors <- function(lap_lifted, lap_other)
{
  # overlap the anchors from lap_lifted with the actual coordinates observed in the other species
  lap_lifted_gr <- with(lap_lifted, GRanges(other_chrom, IRanges(other_start + 1L, other_end)))
  if ("anchor_chr" %in% names(lap_other))
    lap_other_gr <- with(lap_other, GRanges(sub("^chr", "", anchor_chr), IRanges(anchor_start + 1L, anchor_end)))
  else
    lap_other_gr <- with(lap_other, GRanges(chrom, IRanges(start + 1L, end)))

  # prevent the warnings about no sequence levels in common
  ov <- suppressWarnings(findOverlaps(lap_lifted_gr + segment_margin, lap_other_gr))

  dt <- data.table(
    ind = queryHits(ov),
    other_loop_id = lap_other$loop_id[subjectHits(ov)],
    other_anchor = lap_other$anchor[subjectHits(ov)]
  )

  lap_lifted[, ind := .I]
  lap_lifted <- merge(lap_lifted, dt, by = "ind", all.x = TRUE, sort = FALSE)
  lap_lifted[, ind := NULL]

  return(lap_lifted)
}

#
#  Convert loop anchor data table to a loop data table with one row per loop
#

lap_to_loops <- function(lap)
{
  lap_A1 <- lap[anchor == "A1", ]
  lap_A1_sel <- lap_A1[, list(loop_id, aggregated_id,
    A1_midpoint = if ("anchor_summit" %in% names(.SD)) anchor_summit else midpoint)]
  if ("type" %in% names(lap))
    lap_A1_sel$A1_type = lap_A1$type
  if ("DHS_nearest_gene_symbol" %in% names(lap))
    lap_A1_sel$A1_DHS_nearest_gene_symbol = lap_A1$DHS_nearest_gene_symbol

  lap_A2 <- lap[anchor == "A2", ]
  lap_A2_sel <- lap_A2[, list(loop_id, aggregated_id,
    A2_midpoint = if ("anchor_summit" %in% names(.SD)) anchor_summit else midpoint)]
  if ("type" %in% names(lap))
    lap_A2_sel$A2_type = lap_A2$type
  if ("DHS_nearest_gene_symbol" %in% names(lap))
    lap_A2_sel$A2_DHS_nearest_gene_symbol = lap_A2$DHS_nearest_gene_symbol

  loops <- merge(lap_A1_sel, lap_A2_sel, sort = F)
  stopifnot(length(unique(loops$loop_id)) == nrow(loops))

  loops[, distance := abs(A2_midpoint - A1_midpoint)]
  loops[, A1_midpoint := NULL]
  loops[, A2_midpoint := NULL]

  return(loops)
}

#
#  Lift the source genome loops to the other species, taking all possible combinations of lifted anchors
#  and overlapping them against annotated loops in the other species
#

lift_loops <- function(match)
{
  # convert loop anchor data table to loop data table
  match$matches <- lap_to_loops(match$lap)
  match$matches$n_candidate <- NA
  match$matches$n_observed <- NA

  # data table for anchors of matched (observed) loops in the other species
  match$candidates <- data.table()

  # consider each source genome loop separately
  for (i in seq_len(nrow(match$matches)))
  {
    # focus on a given loop in the source genome
    li <- match$matches$loop_id[i]

    # take the loop anchor lifts, i.e. candidates for both anchors after lifting them to the other species
    lal1 <- match$lap_lifted[loop_id == li & anchor == "A1", ]
    lal2 <- match$lap_lifted[loop_id == li & anchor == "A2", ]
    match$matches$n_candidate[i] <- length(unique(lal1$chain_id)) * length(unique(lal2$chain_id))
    match$matches$n_observed[i] <- 0L

    # candidate loops (dt) are obtained as the cartesian product
    # of loop anchor lifts for anchor 1 (lal1) and anchor 2 (lal2)
    dt <- merge(lal1, lal2, by = "loop_id", allow.cartesian = TRUE)
    dt[, anchor.x := NULL]
    dt[, anchor.y := NULL]

    # consider a candidate loop, i.e. a row in dt, with candidate anchors defined by *.x and *.y
    # we have a match if other_loop_id is the same for the two lifted anchors,
    # and these lifted anchors correspond to two different anchors of the same loop in the other species
    dt[, other_loop_id := ifelse(other_loop_id.x == other_loop_id.y & other_anchor.x != other_anchor.y,
      other_loop_id.x, NA)]

    # change column names to denote the anchors by A1 and A2 instead of *.x and *.y
    setnames(dt, "chain_id.x", "A1_chain_id")
    setnames(dt, c("chrom.x", "start.x", "end.x"), c("A1_chrom", "A1_start", "A1_end"))
    setnames(dt, c("other_chrom.x", "other_start.x", "other_end.x"), c("A1_other_chrom", "A1_other_start", "A1_other_end"))
    setnames(dt, "chain_range.x", "A1_chain_range")
    setnames(dt, "other_chain_range.x", "A1_other_chain_range")
    setnames(dt, "other_chain_strand.x", "A1_other_chain_strand")
    setnames(dt, "other_loop_id.x", "A1_other_loop_id")
    setnames(dt, "other_anchor.x", "A1_other_anchor")
    setnames(dt, "exon_fraction.x", "A1_exon_fraction", skip_absent=TRUE)
    setnames(dt, "other_exon_fraction.x", "A1_other_exon_fraction", skip_absent=TRUE)

    setnames(dt, "chain_id.y", "A2_chain_id")
    setnames(dt, c("chrom.y", "start.y", "end.y"), c("A2_chrom", "A2_start", "A2_end"))
    setnames(dt, c("other_chrom.y", "other_start.y", "other_end.y"), c("A2_other_chrom", "A2_other_start", "A2_other_end"))
    setnames(dt, "chain_range.y", "A2_chain_range")
    setnames(dt, "other_chain_range.y", "A2_other_chain_range")
    setnames(dt, "other_chain_strand.y", "A2_other_chain_strand")
    setnames(dt, "other_loop_id.y", "A2_other_loop_id")
    setnames(dt, "other_anchor.y", "A2_other_anchor")
    setnames(dt, "exon_fraction.y", "A2_exon_fraction", skip_absent=TRUE)
    setnames(dt, "other_exon_fraction.y", "A2_other_exon_fraction", skip_absent=TRUE)

    # add some previously gathered annotations
    dt <- merge(match$matches[i], dt, by = "loop_id", all.y = TRUE, sort = FALSE)
    dt[, n_observed := NULL]

    match$matches$n_observed[i] <- length(unique(na.exclude(dt$other_loop_id)))
    match$candidates <- rbind(match$candidates, dt)
  }

  return(match)
}

#
#  Annotate the matched loops
#

annotate_matched_loops <- function(match)
{
  # put together loop_id and a matched other_loop_id
  dt <- unique(match$candidates[!is.na(other_loop_id), list(loop_id, other_loop_id)])
  stopifnot(!anyDuplicated(dt$loop_id))

  match$matches <- merge(match$matches, dt, by = "loop_id", all.x = TRUE, sort = FALSE)

  lap_other_A1 <- match$lap_other[anchor == "A1", list(loop_id, aggregated_id,
    A1_type = type,
    A1_chrom = if ("anchor_chr" %in% names(.SD)) anchor_chr else chrom,
    A1_start = if ("anchor_start" %in% names(.SD)) anchor_start else start,
    A1_end = if ("anchor_end" %in% names(.SD)) anchor_end else end,
    A1_midpoint = if ("anchor_summit" %in% names(.SD)) anchor_summit else midpoint)]

  lap_other_A2 <- match$lap_other[anchor == "A2", list(loop_id, aggregated_id,
    A2_type = type,
    A2_chrom = if ("anchor_chr" %in% names(.SD)) anchor_chr else chrom,
    A2_start = if ("anchor_start" %in% names(.SD)) anchor_start else start,
    A2_end = if ("anchor_end" %in% names(.SD)) anchor_end else end,
    A2_midpoint = if ("anchor_summit" %in% names(.SD)) anchor_summit else midpoint)]

  lap_other_dist <- merge(lap_other_A1, lap_other_A2, sort = FALSE)
  stopifnot(length(unique(lap_other_dist$loop_id)) == nrow(lap_other_dist))

  lap_other_dist[, distance := abs(A2_midpoint - A1_midpoint)]
  lap_other_dist[, A1_midpoint := NULL]
  lap_other_dist[, A2_midpoint := NULL]

  m <- match(match$matches$other_loop_id, lap_other_dist$loop_id)
  match$matches$other_distance <- lap_other_dist$distance[m]
  match$matches$A1_other_chrom <- lap_other_dist$A1_chrom[m]
  match$matches$A1_other_start <- lap_other_dist$A1_start[m]
  match$matches$A1_other_end <- lap_other_dist$A1_end[m]
  match$matches$A2_other_chrom <- lap_other_dist$A2_chrom[m]
  match$matches$A2_other_start <- lap_other_dist$A2_start[m]
  match$matches$A2_other_end <- lap_other_dist$A2_end[m]

  # add information on CE chain ranges
  dt <- match$candidates[!is.na(other_loop_id), list(
    A1_other_anchor = paste(unique(A1_other_anchor), collapse = ","),
    A1_chain_id = paste(unique(A1_chain_id), collapse = ","),
    A1_chain_range = paste(unique(A1_chain_range), collapse = ","),
    A1_other_chain_range = paste(unique(A1_other_chain_range), collapse = ","),
    A1_other_chain_strand = paste(unique(A1_other_chain_strand), collapse = ","),
    A2_other_anchor = paste(unique(A2_other_anchor), collapse = ","),
    A2_chain_id = paste(unique(A2_chain_id), collapse = ","),
    A2_chain_range = paste(unique(A2_chain_range), collapse = ","),
    A2_other_chain_range = paste(unique(A2_other_chain_range), collapse = ","),
    A2_other_chain_strand = paste(unique(A2_other_chain_strand), collapse = ",")
  ), by = "loop_id"]
  match$matches <- merge(match$matches, dt, by = "loop_id", all.x = TRUE, sort = FALSE)

  # Add column with the critical value for exon overlap fraction for chains of CEs used to lift the anchors.
  # I.e. when considering chains of CEs and the fraction of total alignment width that is covered by exons,
  # we could exclude some chains due to excessively high value of exon overlap fraction.
  # At what threshold value would such a limitation cause discarding this loop match between the species?
  if ("A1_exon_fraction" %in% names(match$candidates))
  {
    dt <- match$candidates[!is.na(other_loop_id), list(
      A1_exon_fraction = max(A1_exon_fraction, -Inf),
      A1_other_exon_fraction = max(A1_other_exon_fraction, -Inf),
      A2_exon_fraction = max(A2_exon_fraction, -Inf),
      A2_other_exon_fraction = max(A2_other_exon_fraction, -Inf),
      critical_exon_fraction = min(pmax(A1_exon_fraction, A1_other_exon_fraction, A2_exon_fraction, A2_other_exon_fraction), Inf)), by = "loop_id"]
    match$matches <- merge(match$matches, dt, by = "loop_id", all.x = TRUE, sort = FALSE)
  }

  # Add information on loop orientation in the other species.
  # Possible orientations are as follows:
  #   A1> ... A2> or <A2 ... <A1  original
  #   <A1 ... A2> or <A2 ... A1>  single-inverted
  #   A1> ... <A2 or A2> ... <A1  single-inverted
  #   <A1 ... <A2 or A2> ... A1>  double-inverted
  # Note that all the loops have original orientation in the source genome.

  # here we use the invariant that A2 comes after A1 in the genomic order in both species (checked above)
  other_ori <- with(match$matches, paste(A1_other_anchor, A1_other_chain_strand, A2_other_anchor, A2_other_chain_strand))
  match$matches$other_orientation <- factor(NA, c("original", "single_inv", "double_inv"))
  match$matches$other_orientation[other_ori %in% c("A1 + A2 +", "A2 - A1 -")] <- "original"
  match$matches$other_orientation[other_ori %in% c("A1 - A2 +", "A2 - A1 +", "A1 + A2 -", "A2 + A1 -")] <- "single_inv"
  match$matches$other_orientation[other_ori %in% c("A1 - A2 -", "A2 + A1 +")] <- "double_inv"

  return(match)
}

#
#  Match the source genome anchors to the other species
#

match_loops <- function(genome, other_genome, lap, lap_other, cne, cne_exon_fraction = NULL)
{
  # lift the anchors to the other species
  lap_lifted <- lift_anchors(lap, cne)
  # overlap the lifted anchors with the actual anchors in the other species
  lap_lifted <- overlap_lifted_anchors(lap_lifted, lap_other)
  # add the information on exon overlap fraction for chains of CEs used
  if (!is.null(cne_exon_fraction))
    lap_lifted <- merge(lap_lifted, cne_exon_fraction, by = "chain_id", all.x = TRUE, sort = FALSE)

  # update the original loop annotations to add number of candidates (counted as the number of distinct CE chains) per anchor
  dt <- lap_lifted[, list(n_candidate = length(unique(chain_id))), by = c("loop_id", "anchor")]
  lap <- merge(lap, dt, by = c("loop_id", "anchor"), all.x = TRUE, sort = FALSE)
  lap$n_candidate[is.na(lap$n_candidate)] <- 0L

  match <- list(
    genome = genome,
    other_genome = other_genome,
    lap = lap,
    lap_other = lap_other,
    lap_lifted = lap_lifted
  )

  match <- lift_loops(match)
  match <- annotate_matched_loops(match)
  return(match)
}

#
#  Print some summary statistics
#

print_statistics <- function(match)
{
  message("Out of ", nrow(match$lap), " ", match$genome, " loop anchors, ", sum(match$lap$n_candidate > 0), " have a candidate in ", match$other_genome)
  message("On average, a ", match$genome, " loop anchor has ", mean(match$lap$n_candidate), " candidates in ", match$other_genome)
  message("Out of ", nrow(match$matches), " ", match$genome, " loops, ",
    sum(match$matches$n_candidate > 0), " have a candidate in ", match$other_genome, " (i.e. both anchors matched), and ",
    sum(match$matches$n_observed > 0), " match a loop observed in ", match$other_genome)

  # confirm that the lifted loops are cis-chromosomal (after all, this is an important finding!)
  if (isTRUE(all.equal(match$candidates$A1_other_chrom, match$candidates$A2_other_chrom)))
    message("All the candidate loops in ", match$other_genome, " are cis-chromosomal.")

  message("Orientations of matched loops in ", match$other_genome, ":")
  print(table(match$matches[!is.na(other_loop_id), "other_orientation"], useNA = "ifany"))
  message()

  message("Out of ", length(unique(match$matches$aggregated_id)), " aggregated ", match$genome, " loops, ",
    length(unique(match$matches[n_candidate > 0]$aggregated_id)), " have a candidate in ", match$other_genome, " (i.e. both anchors matched), and ",
    length(unique(match$matches[n_observed > 0]$aggregated_id)), " have a match observed in ", match$other_genome)
  message()
}

#
#  Save loop anchor matches, lifted loops (candidates) nad matched loops (matches) to separate files
#

save_loops <- function(match)
{
  write.table(match$lap, file = paste0("data/loops/loops_", match$genome, "_to_", match$other_genome, "_anchors.tsv"),
    quote = F, sep = "\t", row.names = F, col.names = T)

  write.table(match$candidates, file = paste0("data/loops/loops_", match$genome, "_to_", match$other_genome, "_candidates.tsv"),
    quote = F, sep = "\t", row.names = F, col.names = T)

  write.table(match$matches, file = paste0("data/loops/loops_", match$genome, "_to_", match$other_genome, "_matches.tsv"),
    quote = F, sep = "\t", row.names = F, col.names = T)

  save(match, file = paste0("data/loops/loops_", match$genome, "_to_", match$other_genome, ".Rdata"))
}

#
#  Calculate the fraction of CEs overlapped by exons
#

calculate_cne_exon_fraction <- function(cne, exons, exons_other)
{
  cne_exon_fraction <- data.table(
    chain_id = elementMetadata(cne)$chain_id,
    first_width = width(first(cne)),
    second_width = width(second(cne))
  )

  overlap_width <- function(query, subject)
  {
    hits <- findOverlaps(query, subject)
    overlaps <- pintersect(query[queryHits(hits)], subject[subjectHits(hits)])
    result <- rep(0L, length(query))
    dt <- data.table(ind = queryHits(hits), width = width(overlaps))
    dt <- dt[, list(width = sum(width)), by = "ind"]
    result[dt$ind] <- dt$width
    return(result)
  }

  cne_exon_fraction$first_exon_overlap_width <- overlap_width(first(cne), GRanges(exons))
  cne_exon_fraction$second_exon_overlap_width <- overlap_width(second(cne), GRanges(exons_other))

  cne_exon_fraction <- cne_exon_fraction[, list(
    exon_fraction = sum(first_exon_overlap_width) / sum(first_width),
    other_exon_fraction = sum(second_exon_overlap_width) / sum(second_width)
  ), by = "chain_id"]

  return(cne_exon_fraction)
}

#
#  Main loop
#

match_and_annotate_loops <- function(genome, other_genome)
{
  # Read CE (conserved element) annotations
  load(paste0("data/lastz/ceFinalAnnotated.", genome, ".", other_genome, ".Rdata")) # cneFinal, exons, exons_other
  cne_exon_fraction <- calculate_cne_exon_fraction(cneFinal, exons, exons_other)

  # Read loop annotations
  lap <- extract_loops(genome)
  lap_other <- extract_loops(other_genome)

  match <- match_loops(genome, other_genome, lap, lap_other, cneFinal, cne_exon_fraction)
  print_statistics(match)
  save_loops(match)
}


match_and_annotate_loops("D_mel", "D_vir")
match_and_annotate_loops("D_vir", "D_mel")