options(warn = 1)

# Read config file
config <- yaml::yaml.load_file("config/config.yml")

library(data.table)
library(GenomicRanges)
library(ggplot2)

genome <- "D_mel"
other_genome <- "D_vir"

# read the "reverse" loop matching (from other_genome to genome)
load(paste0("data/loops/loops_", other_genome, "_to_", genome, ".Rdata")) # match
other_match <- match
# now mask the loaded variables with the "actual" loop matching (from genome to other_genome)
load(paste0("data/loops/loops_", genome, "_to_", other_genome, ".Rdata")) # match

#
#  Helper function to append frequency to each class name for labelling
#

count_levels <- function(x, labels = NULL)
{
  tab <- table(x, useNA = "ifany")
  vec <- names(tab)
  m <- match(x, vec)

  if (!is.null(labels))
  {
    mm <- match(vec, names(labels))
    vec[!is.na(mm)] <- labels[na.exclude(mm)]
  }

  lab <- paste0(vec, " (", format(tab, trim = TRUE, big.mark = ","), ")")
  return(factor(lab[m], lab, ordered = is.ordered(x)))
}

#
#  Identify matched *conserved* (looping) / matched *non-looping* (non-conserved) / *unmatched* loops
#

# annotate match$matches
match$matches[, match_status := factor(ifelse(n_observed > 0, "conserved",
  ifelse(n_candidate > 0, "non-looping", "unmatched")), c("conserved", "non-looping", "unmatched"))]
stopifnot(with(match$matches, !is.na(other_loop_id) == (match_status == "conserved")))

conserved_loop_ids <- match$matches[match_status == "conserved", ]$loop_id
nonlooping_loop_ids <- match$matches[match_status == "non-looping", ]$loop_id

# extract unique matched non-looping loops from the candidate matches, and annotate them
unique_annotated_nonlooping <- function(nonlooping)
{
  # calculate the distance of non-looping loops in D. vir.
  nonlooping[, other_distance := abs(as.integer((A1_other_start + A1_other_end) / 2) - as.integer((A2_other_start + A2_other_end) / 2))]

  # determine the orientation of non-looping loops in D. vir.
  other_ori <- with(nonlooping, paste(
    ifelse(A1_other_end < A2_other_start, "A1", "A2"), A1_other_chain_strand,
    ifelse(A1_other_end < A2_other_start, "A2", "A1"), A2_other_chain_strand))
  nonlooping$other_orientation <- factor(NA, c("original", "single_inv", "double_inv"))
  nonlooping$other_orientation[other_ori %in% c("A1 + A2 +", "A2 - A1 -")] <- "original"
  nonlooping$other_orientation[other_ori %in% c("A1 - A2 +", "A2 - A1 +", "A1 + A2 -", "A2 + A1 -")] <- "single_inv"
  nonlooping$other_orientation[other_ori %in% c("A1 - A2 -", "A2 + A1 +")] <- "double_inv"
  rm(other_ori)

  # take each pair of loop anchors once only
  nonlooping_unique <- unique(nonlooping[, c("loop_id", "distance",
    "A1_chrom", "A1_start", "A1_end",
    "A1_other_chrom", "A1_other_start", "A1_other_end", "A1_other_chain_strand",
    "A2_chrom", "A2_start", "A2_end",
    "A2_other_chrom", "A2_other_start", "A2_other_end", "A2_other_chain_strand",
    "other_distance", "other_orientation")])

  return(nonlooping_unique)
}

# extract the matched non-looping loops (matched from D. mel. to D. vir, but none of the matches is looping)
# note: if a loop has a conserved match, we disregard all its non-looping (non-conserved) matches
nonlooping <- match$candidates[!loop_id %in% conserved_loop_ids, ]
stopifnot(nonlooping$loop_id %in% nonlooping_loop_ids)
nonlooping_unique <- unique_annotated_nonlooping(nonlooping)

#
#  Identify backward matched *conserved* (looping) / backward matched *non-looping* (non-conserved)
#  / backward *unmatched* loops (D. vir. to D. mel.)
#

# annotate other_match$matches
other_match$matches[, match_status := factor(ifelse(n_observed > 0, "conserved",
  ifelse(n_candidate > 0, "non-looping", "unmatched")), c("conserved", "non-looping", "unmatched"))]
stopifnot(with(other_match$matches, !is.na(other_loop_id) == (match_status == "conserved")))

other_conserved_loop_ids <- other_match$matches[match_status == "conserved", ]$loop_id
other_nonlooping_loop_ids <- other_match$matches[match_status == "non-looping", ]$loop_id

# check that the loop matching in both directions yields the same result
stopifnot(setequal(other_match$matches[match_status == "conserved", ]$other_loop_id, conserved_loop_ids))
stopifnot(setequal(match$matches[match_status == "conserved", ]$other_loop_id, other_conserved_loop_ids))

# extract the backward matched non-looping loops
# (matched from D. vir. to D. mel, but none of the matches is looping)
other_nonlooping <- other_match$candidates[!loop_id %in% other_conserved_loop_ids, ]
stopifnot(other_nonlooping$loop_id %in% other_nonlooping_loop_ids)
other_nonlooping_unique <- unique_annotated_nonlooping(other_nonlooping)

#
#  Save the conserved loops for use in Juicebox
#

save_conserved_loops <- function(dt, file)
{
  dt[, color := "51,204,0"]
  setnames(dt, "A1_other_chrom", "chr1")
  setnames(dt, "A1_other_start", "x1")
  setnames(dt, "A1_other_end", "x2")
  setnames(dt, "A2_other_chrom", "chr2")
  setnames(dt, "A2_other_start", "y1")
  setnames(dt, "A2_other_end", "y2")
  setcolorder(dt, c("chr1", "x1", "x2", "chr2", "y1", "y2", "color"))

  write.table(dt, file = file, quote = F, sep = "\t", row.names = F, col.names = T)
}

save_conserved_loops(match$matches[match_status == "conserved", ], paste0("data/loops/long_range_loops_", genome, "_matched_conserved_in_", other_genome, ".tsv"))
save_conserved_loops(other_match$matches[match_status == "conserved", ], paste0("data/loops/long_range_loops_", other_genome, "_matched_conserved_in_", genome, ".tsv"))

#
#  Save the matched non-looping loops
#

write.table(nonlooping_unique, file = paste0("data/loops/loops_", genome, "_matched_nonlooping_in_", other_genome, ".tsv"), quote = F, sep = "\t", row.names = F, col.names = T)
write.table(other_nonlooping_unique, file = paste0("data/loops/loops_", other_genome, "_matched_nonlooping_in_", genome, ".tsv"), quote = F, sep = "\t", row.names = F, col.names = T)

# for use in Juicebox
save_nonlooping_loops <- function(dt, file)
{
  dt[, color := "51,51,51"]
  setnames(dt, "A1_other_chrom", "chr1")
  setnames(dt, "A1_other_start", "x1")
  setnames(dt, "A1_other_end", "x2")
  setnames(dt, "A2_other_chrom", "chr2")
  setnames(dt, "A2_other_start", "y1")
  setnames(dt, "A2_other_end", "y2")
  setcolorder(dt, c("chr1", "x1", "x2", "chr2", "y1", "y2", "color"))

  write.table(dt, file = file, quote = F, sep = "\t", row.names = F, col.names = T)
}

save_nonlooping_loops(nonlooping_unique, paste0("data/loops/long_range_loops_", genome, "_matched_nonlooping_in_", other_genome, ".tsv"))
save_nonlooping_loops(other_nonlooping_unique, paste0("data/loops/long_range_loops_", other_genome, "_matched_nonlooping_in_", genome, ".tsv"))

#
#  Violin plot: loop anchor width
#

dt <- match$lap
dt[, species := "D. melanogaster"]
setnames(dt, c("anchor_chr", "anchor_start", "anchor_end", "anchor_summit"),
  c("chrom", "start", "end", "midpoint"))

dt_other <- match$lap_other
dt_other[, species := "D. virilis"]

dt_combined <- rbind(dt, dt_other, fill = TRUE)

# take each loop anchor once only
dt_combined_unique <- unique(dt_combined[, c("chrom", "start", "end", "species")])

pdf("analysis/plots/loops_D_mel_to_D_vir_anchor_width.pdf", width = 4, height = 2)

p <- ggplot(dt_combined_unique,
  aes(x = (end - start) / 1e3, y = count_levels(species))) +
  geom_violin() +
  scale_x_log10() +
  xlab("Anchor width (kb)") +
  ylab(NULL) +
  NULL
print(p)

dev.off()

message("Numbers of unique loop anchors:")
print(dt_combined_unique[, list(.N), by = "species"])

#
#  Violin plots: loop distance
#

dt2_combined <- dt_combined[, list(
  distance = midpoint[which(anchor == "A2")] - midpoint[which(anchor == "A1")]
), by = c("loop_id", "species")]
# check that A2 comes after A1 in the genomic order in both species
stopifnot(dt2_combined$distance > 0)

# add information on loop status (conserved/non-looping/unmatched)
dt2_combined <- merge(dt2_combined,
  rbind(match$matches[, c("loop_id", "match_status")], other_match$matches[, c("loop_id", "match_status")]),
  by = "loop_id", all.x = TRUE, sort = FALSE)
stopifnot(!is.na(dt2_combined$match_status))
setkey(dt2_combined, species, match_status)

pdf("analysis/plots/loops_D_mel_to_D_vir_distance_violinplot.pdf", width = 4, height = 2)

p <- ggplot(dt2_combined, aes(x = distance / 1e6, y = count_levels(species))) +
  geom_violin() +
  xlab("Loop distance (Mb)") +
  ylab(NULL) +
  NULL
print(p)

dev.off()

message("Numbers of unique loops:")
print(dt2_combined[, list(.N), by = c("species", "match_status")])

pdf("analysis/plots/loops_D_mel_to_D_vir_distance_violinplot_by_matched.pdf", width = 4, height = 3)

labels <- c(`non-looping` = "non\uadlooping in D. vir.")
p <- ggplot(dt2_combined, aes(x = distance / 1e6, y = count_levels(match_status, labels = labels), fill = match_status)) +
  facet_grid(rows = vars(species)) +
  geom_violin() +
  geom_jitter(alpha = 0.4) +
  scale_y_discrete(limits = rev) +
  scale_fill_manual(values = c("#e41a1c", "#ff7f00", "darkgray"), guide = "none") +
  theme(legend.position = "bottom") +
  xlab("Loop distance (Mb)") +
  ylab(NULL) +
  NULL
print(p)

dev.off()

message("Does the distance of D. melanogaster loops differ between conserved and non-looping loops? p = ",
  with(dt2_combined[species == "D. melanogaster"], wilcox.test(distance[match_status == "conserved"], distance[match_status == "non-looping"]))$p.value,
  " (Wilcoxon test)")
message("Does the distance of D. virilis loops differ between conserved and non-looping loops? p = ",
  with(dt2_combined[species == "D. virilis"], wilcox.test(distance[match_status == "conserved"], distance[match_status == "non-looping"]))$p.value,
  " (Wilcoxon test)")

#
#  Plot distances and orientations of matched loops between the species
#

dt <- match$matches[match_status == "conserved", ]
levels(dt$other_orientation) <- c(levels(dt$other_orientation), "complex")
dt$other_orientation[is.na(dt$other_orientation)] <- "complex"
labels <- c(original = "same as in D. mel.", `single_inv` = "one anchor inverted", double_inv = "both anchors inverted")

pdf("analysis/plots/loops_D_mel_to_D_vir_conserved.pdf", width = 5, height = 2.5)

p <- ggplot(dt, aes(x = distance / 1e6, y = other_distance / 1e6, shape = count_levels(other_orientation, labels), color = count_levels(other_orientation, labels))) +
  geom_point(alpha = 0.6) +
  scale_shape_discrete(name = "Orientation in D. vir.\nof conserved loops") +
  scale_color_brewer(palette = "Set1", name = "Orientation in D. vir.\nof conserved loops") +
  xlab("Distance in D. mel. (Mb)") +
  ylab("Distance in D. vir. (Mb)") +
  coord_fixed(xlim = c(0, NA), ylim = c(0, NA)) +
  guides(color = guide_legend(ncol = 1), fill = guide_legend(ncol = 1))
print(p)

dev.off()

dt <- rbind(dt,
  data.table(nonlooping_unique, match_status = "nonlooping_Dmel"),
  with(other_nonlooping_unique, data.table(distance = other_distance, other_distance = distance, other_orientation = other_orientation, match_status = "nonlooping_Dvir")),
  fill = T)
match_labels <- c(
  conserved = "conserved loops",
  nonlooping_Dmel = "non\uadlooping matches\nof D. mel. loops in D. vir.",
  nonlooping_Dvir = "non\uadlooping matches\nof D. vir. loops in D. mel.")

pdf("analysis/plots/loops_D_mel_to_D_vir_conserved_and_nonlooping.pdf", width = 8, height = 3)

p <- ggplot(dt, aes(x = distance / 1e6, y = other_distance / 1e6, shape = count_levels(other_orientation, labels), color = count_levels(other_orientation, labels))) +
  facet_grid(cols = vars(count_levels(match_status, labels = match_labels))) +
  geom_point(alpha = 0.6) +
  scale_shape_discrete(name = "Orientation in D. vir.") +
  scale_color_brewer(palette = "Set1", name = "Orientation in D. vir.") +
  xlab("Distance in D. mel. (Mb)") +
  ylab("Distance in D. vir. (Mb)") +
  coord_fixed(xlim = c(0, NA), ylim = c(0, NA)) +
  guides(color = guide_legend(ncol = 1), fill = guide_legend(ncol = 1))
print(p)

dev.off()

#
#  Distributions of the number of anchor candidate matches and loop candidate matches
#

dt <- unique(match$lap[, c("peak_id", "anchor_chr", "anchor_start", "anchor_end", "n_candidate")])
stopifnot(!anyDuplicated(dt$peak_id))

pdf("analysis/plots/loops_D_mel_to_D_vir_anchor_matches.pdf", width = 2.5, height = 2.5)

p <- ggplot(dt, aes(x = n_candidate)) +
  geom_histogram(binwidth = 1) +
  xlab("Number of candidate matches\nof D. mel. anchors in D. vir.") +
  ylab("Count") +
  NULL
print(p)

dev.off()

dt2 <- copy(match$matches)
dt2[, matched := loop_id %in% conserved_loop_ids]

pdf("analysis/plots/loops_D_mel_to_D_vir_loop_matches.pdf", width = 3, height = 2.5)

p <- ggplot(dt2, aes(x = n_candidate, fill = matched)) +
  geom_histogram(binwidth = 1) +
  scale_fill_manual(name = "Loop status", values = c("darkgray", "#e41a1c"),
    labels = c(`TRUE` = "matched", `FALSE` = "unmatched")) +
  xlab("Number of candidate matches\nof D. mel. loops in D. vir.") +
  ylab("Count") +
  NULL
print(p)

dev.off()

#
#  Permuted variants of the original D. mel. loops
#

# permuted loops are less likely to be observed in the other species
stats_permuted <- fread(paste0("data/loops/loops_D_mel_to_", other_genome, "_permuted.tsv"), header = T, sep = "\t")

pdf("analysis/plots/loops_D_mel_to_D_vir_permuted.pdf", width = 3, height = 2.5)

p <- ggplot(stats_permuted, aes(x = n_observed)) +
  geom_histogram(binwidth = 1) +
  geom_vline(data = data.table(n_observed = mean(stats_permuted$n_observed)), aes(xintercept = n_observed), lty = 2) +
  geom_vline(data = data.table(n_observed = n_obs), aes(xintercept = n_observed), color = "#ff7f00") +
  geom_label(x = n_obs, y = Inf, color = "#ff7f00", hjust = 0.5, vjust = 0,
    label = "Number of actual D. mel. loops\nconserved in D. vir.") +
  coord_cartesian(xlim = c(0, nrow(match$matches)), clip = "off") +
  ggtitle("\n") +
  theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "points")) +
  xlab("Number of permuted D. mel. loops\nconserved in D. vir.") +
  ylab("Count") +
  NULL
print(p)

dev.off()

message("Actual loops have ", sum(match$matches$n_observed > 0) / mean(stats_permuted$n_observed), " times more matched loops than randomly permuted loops")

#
#  Exon ratio of anchors of matched loops
#

dt <- with(match$matches[match_status == "conserved", ], rbind(
  data.table(
    species = "D. melanogaster",
    chain_range = c(A1_chain_range, A2_chain_range),
    exon_fraction = c(A1_exon_fraction, A2_exon_fraction)
  ),
  data.table(
    species = "D. virilis",
    chain_range = c(A1_other_chain_range, A2_other_chain_range),
    exon_fraction = c(A1_other_exon_fraction, A2_other_exon_fraction)
  )
))

# take each loop anchor once only
dt <- unique(dt)

message("Median fraction of aligned nucleotides at a loop anchor overlapped by gene exons: ", median(dt$exon_fraction))
message("Anchors having the fraction of aligned nucleotides overlapped by gene exons > 0.6:")
print(dt[exon_fraction > 0.6, ])

pdf("analysis/plots/loops_D_mel_to_D_vir_exon_ratio.pdf", width = 4.5, height = 2)

p <- ggplot(dt, aes(x = exon_fraction, y = count_levels(species))) +
  geom_violin() +
  geom_boxplot(width = 0.2) +
  xlab("Fraction of aligned nucleotides at a loop anchor\noverlapped by gene exons") +
  ylab(NULL) +
  NULL
print(p)

dev.off()
