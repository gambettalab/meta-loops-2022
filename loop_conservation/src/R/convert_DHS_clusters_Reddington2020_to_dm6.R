options(warn = 1)

library(data.table)
library(GenomicRanges)
library(rtracklayer)

dir <- "data/DNase-seq_Reddington2020"
chain_file <- "data/genome/D.melanogaster/dm3/liftOver/dm3ToDm6.over.chain.gz"

#
#  read the 59,864 DHS clusters
#

# 'tr' to convert the end-of-line characters from Mac to Unix
# 'sed' to avoid the curse of CAD (non-printable characters and other gibberish)
# 'grep' to remove comment lines
dt <- fread(cmd = paste0("cat ", dir, "/1-s2.0-S1534580720307991-mmc3.csv | tr '\r' '\n' | sed -e 's/\013/_/g' -e 's/\x89\xf6\xd5/-/g' -e \"s/\x89\xdb/'/g\" | grep -av '^#\\|^\"#'"), sep = ",")
# 'chr', 'start' and 'end' give the genomic coordinates (dm3) of each DHS cluster in BED format.
# Each DHS cluster is given a unique name ('uname') made from the 'chr' and the 'summit' position.
# The 'summit' (1-indexed) is determined during DHS cluster creation as the position with maximum DHS-summit density.

stopifnot(nrow(dt) == 59864)

#
#  convert the DHS coordinates to dm6
#

liftOver_dm3_to_dm6 <- function(dt) # data table with columns 'chr', 'start' and 'end' in BED format
{
  ftmp1 <- tempfile()
  ftmp2 <- tempfile()
  ftmp3 <- tempfile()

  bed_dm3 <- with(dt, data.table(chrom = chr, start = start, end = end, name = ".", score = seq_len(nrow(dt)), strand = "."))
  write.table(bed_dm3, file = ftmp1, quote = F, sep = '\t', row.names = F, col.names = F)
  system(paste("liftOver", ftmp1, chain_file, ftmp2, ftmp3))
  # stopifnot(file.info(ftmp3)$size == 0)
  bed_dm6 <- fread(ftmp2, header = F, sep = "\t")
  names(bed_dm6) <- c("chrom", "start", "end", "name", "score", "strand")
  unlink(c(ftmp1, ftmp2, ftmp3))

  # check that the rows were mapped uniquely (although some might have been lost)
  stopifnot(!anyDuplicated(bed_dm6$score))
  # print(summary(with(bed_dm6, end - start) - with(bed_dm3, end - start)))

  # create a new data table with the new coordinates
  dt_dm6 <- copy(dt)
  m <- match(bed_dm3$score, bed_dm6$score)
  dt_dm6$chr <- bed_dm6$chrom[m]
  dt_dm6$start <- bed_dm6$start[m]
  dt_dm6$end <- bed_dm6$end[m]
  return(dt_dm6)
}

liftOver_summit_dm3_to_dm6 <- function(dt) # data table with columns 'chr' and 'summit' (1-based)
{
  ftmp1 <- tempfile()
  ftmp2 <- tempfile()
  ftmp3 <- tempfile()

  bed_dm3 <- with(dt, data.table(chrom = chr, start = summit - 1L, end = summit, name = ".", score = seq_len(nrow(dt)), strand = "."))
  write.table(bed_dm3, file = ftmp1, quote = F, sep = '\t', row.names = F, col.names = F)
  system(paste("liftOver", ftmp1, chain_file, ftmp2, ftmp3))
  # stopifnot(file.info(ftmp3)$size == 0)
  bed_dm6 <- fread(ftmp2, header = F, sep = "\t")
  names(bed_dm6) <- c("chrom", "start", "end", "name", "score", "strand")
  unlink(c(ftmp1, ftmp2, ftmp3))

  # check that the rows were mapped uniquely (although some might have been lost)
  stopifnot(!anyDuplicated(bed_dm6$score))
  stopifnot(with(bed_dm6, end - start) == 1L)

  # create a new data table with the new coordinates
  dt_dm6 <- copy(dt)
  m <- match(bed_dm3$score, bed_dm6$score)
  dt_dm6$chr <- bed_dm6$chrom[m]
  dt_dm6$summit <- bed_dm6$end[m]
  return(dt_dm6)
}

dt_dm6 <- liftOver_dm3_to_dm6(dt)
dt_dm6$summit <- liftOver_summit_dm3_to_dm6(dt)$summit
dt_dm6[, uname := paste0(chr, ":", summit)]

# take only the DHS clusters that mapped successfully
dt_dm6 <- dt_dm6[!is.na(chr)]
setkey(dt_dm6, chr, start)
message("Out of ", nrow(dt), " DHS clusters, ", nrow(dt_dm6), " were mapped to dm6.")
# Note: There are 4 (unrelated) overlaps between these DHS clusters; all on chrM.

#
#  save the converted annotations
#

write.table(dt_dm6, file = paste0(dir, "/DHS_clusters_dm6.tsv"), quote = F, sep = '\t', row.names = F)
export.gff3(GRanges(dt_dm6), con = paste0(dir, "/DHS_clusters_dm6.gff"))
