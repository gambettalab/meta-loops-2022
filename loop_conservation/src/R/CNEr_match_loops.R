options(warn = 1)

# Read config file
config <- yaml::yaml.load_file("config/config.yml")

library(data.table)
library(GenomicRanges)
library(ggplot2)

#
#  Extract the loops in both species
#

source("src/R/CNEr_match_loops_functions.R")

extract_loops <- function(genome)
{
  loops <- fread(paste0("data/loops/loops_", genome, "_annotated.tsv"), header = T, sep = "\t")
  return(loops)
}

#
#  Print some summary statistics
#

print_statistics <- function(match)
{
  message("Out of ", length(unique(match$lap$anchor_id)), " ", match$genome, " loop anchors, ", length(unique(match$lap[n_candidate > 0, ]$anchor_id)), " have a candidate in ", match$other_genome)
  message("On average, a ", match$genome, " loop anchor has ", mean(unique(match$lap[, c("anchor_id", "n_candidate")])$n_candidate), " candidates in ", match$other_genome)
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
