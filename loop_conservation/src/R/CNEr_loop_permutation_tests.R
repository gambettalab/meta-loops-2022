options(warn = 1)

# Read config file
config <- yaml::yaml.load_file("config/config.yml")

library(data.table)
library(GenomicRanges)

#
#  Read the loop data, and the synteny alignment chains
#

source("src/R/CNEr_match_loops.R")

genome <- "D_mel"
other_genome <- "D_vir"

load(paste0("data/loops/loops_", genome, "_to_", other_genome, ".Rdata")) # match
load(paste0("data/lastz/ceFinalAnnotated.", genome, ".", other_genome, ".Rdata")) # cneFinal

#
#  Compare against permuted loops (randomly matching anchors, within each chromosome separately)
#

permute_loops <- function(lap)
{
  permute_chrom <- function(dt)
  {
    perm <- sample.int(nrow(dt))
    dt$loop_id <- dt$loop_id[perm]
    dt$aggregated_id <- dt$aggregated_id[perm]
    dt$anchor <- dt$anchor[perm]
    dt[, peak_id := NA]
    dt[, partner_peak_id := NA]
    return(dt)
  }
  return(lap[, permute_chrom(.SD), by = "anchor_chr"])
}

#
#  Check if the permuted loops are less likely to be observed in the other species
#

calculate_stats_permuted <- function(match)
{
  message(".", appendLF = F)
  lap_permuted <- permute_loops(match$lap)
  match_permuted <- match_loops(match$genome, match$other_genome, lap_permuted, match$lap_other, cneFinal)
  stats_permuted <- data.table(
    n_candidate = sum(match_permuted$matches$n_candidate > 0),
    n_observed = sum(match_permuted$matches$n_observed > 0),
    n_candidate_aggregated = length(unique(match_permuted$matches[n_candidate > 0]$aggregated_id)),
    n_observed_aggregated = length(unique(match_permuted$matches[n_observed > 0]$aggregated_id))
  )
  return(stats_permuted)
}

stats_permuted <- as.data.table(do.call(rbind, mclapply(1:1000, function(i) calculate_stats_permuted(match))))

write.table(stats_permuted, file = paste0("data/loops/loops_D_mel_to_", other_genome, "_permuted.tsv"), quote = F, sep = "\t", row.names = F, col.names = T)

stats <- data.table(
  enrich_candidate = sum(match$matches$n_candidate > 0) / mean(stats_permuted$n_candidate),
  enrich_observed = sum(match$matches$n_observed > 0) / mean(stats_permuted$n_observed),
  enrich_candidate_aggregated = length(unique(match$matches[n_candidate > 0]$aggregated_id)) / mean(stats_permuted$n_candidate_aggregated),
  enrich_observed_aggregated = length(unique(match$matches[n_observed > 0]$aggregated_id)) / mean(stats_permuted$n_observed_aggregated)
)

write.table(stats, file = paste0("data/loops/loops_D_mel_to_", other_genome, "_permuted_stats.tsv"), quote = F, sep = "\t", row.names = F, col.names = T)
