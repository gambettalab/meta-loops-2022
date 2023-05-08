# Read config file
config <- yaml::yaml.load_file("config/config.yml")

library(data.table)
library(GenomicRanges)
library(parallel)

datasets <- c("HiC_WT_L3_CNS_Rep1", "HiC_WT_L3_CNS_Rep2", "HiC_WT_L3_CNS_Rep4")
anchor_margin <- 5000L # could be further reduced in the downstream analysis
options(mc.cores = 32)


#
#  Read the loop data
#

genome <- "D_mel"
other_genome <- "D_vir"

load(paste0("data/loops/loops_", genome, "_to_", other_genome, ".Rdata")) # match


#
#  Extract the anchors
#

anchor_dt <- unique(match$lap[, c("anchor_id", "anchor_chr", "anchor_start", "anchor_end", "anchor_summit")])
setkey(anchor_dt, anchor_chr, anchor_start)
anchor_gr <- with(anchor_dt, GRanges(anchor_chr, IRanges(anchor_start + 1L, anchor_end)))

# check that the anchors are non-overlapping (this is not required by design, but was the case here)
stopifnot(length(findOverlaps(anchor_gr, drop.self = T)) == 0)


#
#  Extract multi-way interactions for all the anchors
#

out_dir <- "data/HiC/multiway"
if (!dir.exists(out_dir))
  dir.create(out_dir)

extract_fun <- function(dataset, anchor_id)
{
  m <- match(anchor_id, anchor_dt$anchor_id)
  cmd <- paste0("src/py/filter_pairs_to_tsv.py data/HiC/bam/all_reads/", config$genome$symbol, "_", dataset, ".merge.fixmate.rmdup.pairs.gz ", out_dir, "/", config$genome$symbol, "_", dataset, "_", anchor_id, ".tsv.gz ", out_dir, "/", config$genome$symbol, "_", dataset, "_", anchor_id, ".stats chr", anchor_dt$anchor_chr[m], " ", anchor_dt$anchor_start[m] - anchor_margin, " ", anchor_dt$anchor_end[m] + anchor_margin)
  print(cmd)
  system(cmd)
}

for (dataset in datasets)
  mclapply(anchor_dt$anchor_id, extract_fun, dataset = dataset)

system(paste0("touch ", out_dir, "/anchors.txt"))
