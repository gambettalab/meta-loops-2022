library(GenomicFeatures)

suppressWarnings(dir.create("data/genome/D.melanogaster/dm6/dmel_r6.36_FB2020_05/txdb"))

#
#  note that chromosome names from FlyBase (2L, mitochondrion_genome, etc.) differ from UCSC (chr2L, chrM, etc.)
#

seqlengths <- read.table("data/genome/D.melanogaster/dm6/dmel_r6.36_FB2020_05/fasta/dmel-all-chromosome-r6.36.fasta.fai", sep="\t", stringsAsFactors=FALSE)

chromInfo <- Seqinfo(seqnames = seqlengths[[1]],
  seqlength=seqlengths[[2]],
  isCircular=rep(FALSE, nrow(seqlengths)),
  genome="dm6")

#
#  make and save the TxDb
#

txdb <- makeTxDbFromGFF("data/genome/D.melanogaster/dm6/dmel_r6.36_FB2020_05/gtf/dmel-all-r6.36.gtf.gz",
  chrominfo=chromInfo, organism="Drosophila melanogaster")

saveDb(txdb, file="data/genome/D.melanogaster/dm6/dmel_r6.36_FB2020_05/txdb/dmel-all-r6.36.sqlite")
