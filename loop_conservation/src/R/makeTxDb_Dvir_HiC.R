library(GenomicFeatures)

suppressWarnings(dir.create("data/genome/Drosophila_Renschler2019/txdb"))

#
#  note that chromosome names from FlyBase (2L, mitochondrion_genome, etc.) differ from UCSC (chr2L, chrM, etc.)
#

seqlengths <- read.table("data/genome/Drosophila_Renschler2019/GSE120751/suppl/GSE120751_Dvir_HiC.fa.fai", sep="\t", stringsAsFactors=FALSE)

chromInfo <- Seqinfo(seqnames = sub("^scaffold_", "hic_scaffold_", seqlengths[[1]]),
  seqlength=seqlengths[[2]],
  isCircular=rep(FALSE, nrow(seqlengths)),
  genome="Dvir_HiC")

#
#  make and save the TxDb
#

ftmp <- tempfile(fileext=".gtf")
system(paste0("zcat data/genome/Drosophila_Renschler2019/GSE120751/suppl/GSE120751_Dvir_HiC.gtf.gz | grep -Pv '\tstop_codon\t' > ", ftmp))

txdb <- makeTxDbFromGFF(ftmp,
  chrominfo=chromInfo, organism="Drosophila virilis")

unlink(ftmp)

saveDb(txdb, file="data/genome/Drosophila_Renschler2019/txdb/GSE120751_Dvir_HiC.sqlite")
