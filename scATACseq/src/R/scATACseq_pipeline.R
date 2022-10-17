library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(ggplot2)
library(patchwork)
set.seed(1234)

scATACseq.dir <- "data/scATACseq/WTL3CNS/outs"

counts <- Read10X_h5(filename = paste0(scATACseq.dir, "/filtered_peak_bc_matrix.h5"))
metadata <- read.csv(
  file = paste0(scATACseq.dir, "/singlecell.csv"),
  header = TRUE,
  row.names = 1
)

cns <- CreateSeuratObject(
  counts = counts,
  assay = 'peaks',
  project = 'ATAC',
  min.cells = 1,
  meta.data = metadata
)

fragment.path <- paste0(scATACseq.dir, "/fragments.tsv.gz")

cns <- SetFragments(
  object = cns,
  file = fragment.path
)

#
#  Computing QC Metrics
#

pdf("analysis/scATACseq_qc_metrics.pdf", width = 10, height = 6)

# Nucleosome banding pattern; by default, it is is calculated only on chr1 reads
cns <- NucleosomeSignal(object = cns, region = "2L-1-23513712")

cns$pct_reads_in_peaks <- cns$peak_region_fragments / cns$passed_filters * 100

# Note: the following does not work:
# cns$blacklist_ratio <- cns$blacklist_region_fragments / cns$peak_region_fragments
# Workaround:
# These regions were provided by the ENCODE consortium, and we encourage users to cite their paper if you use the regions in your analysis.
seqlevelsStyle(blacklist_dm6) <- 'Ensembl'

cns$blacklist_ratio <- FractionCountsInRegion(
  object = cns,
  assay = 'peaks',
  regions = blacklist_dm6,
  sep = c(":", "-")
)

# create GRanges object with TSS positions
gene.ranges <- genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
seqlevelsStyle(gene.ranges) <- 'Ensembl'
gene.ranges <- keepStandardChromosomes(gene.ranges, pruning.mode = 'coarse')

tss.ranges <- GRanges(
  seqnames = seqnames(gene.ranges),
  ranges = IRanges(start = start(gene.ranges), width = 2),
  strand = strand(gene.ranges)
)

cns <- TSSEnrichment(object = cns, tss.positions = tss.ranges)

p1 <- VlnPlot(cns, c('pct_reads_in_peaks', 'peak_region_fragments', "TSS.enrichment"), pt.size = 0.1)
p2 <- VlnPlot(cns, c('blacklist_ratio', 'nucleosome_signal'), pt.size = 0.1) & scale_y_log10()

print(p1 | p2)

cns$nucleosome_group <- ifelse(cns$nucleosome_signal > 5, 'High', 'Low')
cns$high.prf <- ifelse(cns$peak_region_fragments > 15000, 'High', 'Low')
cns$high.blacklist <- ifelse(cns$blacklist_ratio > 0.05, 'High', 'Low')
cns$high.tss <- ifelse(cns$TSS.enrichment > 2, 'High', 'Low')

print(FragmentHistogram(object = cns, group.by = 'nucleosome_group', region = "2L-1-2000000") + ggtitle("Nucleosome group"))
print(FragmentHistogram(object = cns, group.by = 'high.prf', region = "2L-1-2000000") + ggtitle("Peak region fragments"))
print(FragmentHistogram(object = cns, group.by = 'high.blacklist', region = "2L-1-2000000") + ggtitle("Blacklist ratio"))
print(FragmentHistogram(object = cns, group.by = 'high.tss', region = "2L-1-2000000") + ggtitle("TSS enrichment score"))

print(TSSPlot(cns, group.by = 'nucleosome_group') + ggtitle("Nucleosome group") + NoLegend())
print(TSSPlot(cns, group.by = 'high.prf') + ggtitle("Peak region fragments") + NoLegend())
print(TSSPlot(cns, group.by = 'high.blacklist') + ggtitle("Blacklist ratio") + NoLegend())
print(TSSPlot(cns, group.by = 'high.tss') + ggtitle("TSS enrichment score") + NoLegend())

dev.off()

# Finally we remove cells that are outliers for these QC metrics.

print(cns)

# "Step 4. Quality assessment of each single nucleus. For each single nucleus, we computed ATAC-seq quality metrics such as the fragment length distribution, transcription start site (TSS) enrichment, short-to-mononucleosomal reads ratio, total autosomal reads, and fraction of reads overlapping peaks. We removed nuclei with a) total reads outside the 5%–95% range (34578–226755) of all of the nuclei and b) TSS enrichment of <2.7 (5% tile) from further downstream analysis."

old_cns <- cns

cns <- subset(
  x = old_cns,
  subset =
    # peak_region_fragments > 3000 &
    peak_region_fragments >= quantile(old_cns$peak_region_fragments, 0.05) &
    peak_region_fragments <= quantile(old_cns$peak_region_fragments, 0.95)
    # pct_reads_in_peaks > 15 &
    # blacklist_ratio <= 0.05 &
    # nucleosome_signal <= 10 &
    # TSS.enrichment <= 2
)

print(cns)

#
#  Normalization and linear dimensional reduction
#

cns <- RunTFIDF(cns)
cns <- FindTopFeatures(cns, min.cutoff = 'q0')
cns <- RunSVD(
  object = cns,
  assay = 'peaks',
  reduction.key = 'LSI_',
  reduction.name = 'lsi'
)

pdf("analysis/scATACseq_normalization.pdf", width = 10, height = 6)
print(DepthCor(cns, n = 30))
dev.off()

# there is a very strong correlation between the first LSI component and the total number of counts for the cell, so we will perform downstream steps without this component.

cns <- RunUMAP(object = cns, reduction = 'lsi', dims = 2:30)
cns <- FindNeighbors(object = cns, reduction = 'lsi', dims = 2:30)
cns <- FindClusters(object = cns, verbose = FALSE, algorithm = 3)

# renumber the clusters from 1 (by default they are numbered from 0)
levels(Idents(cns)) <- seq_along(levels(Idents(cns)))
# save the clustering
write.table(data.frame(Barcode = names(Idents(cns)), Cluster = Idents(cns)),
  file = "data/scATACseq/filtered_tracks/clusters.tsv", quote = F, sep = '\t', row.names = F, col.names = F)

pdf("analysis/scATACseq_dim_reduction.pdf", width = 6, height = 6)
print(DimPlot(object = cns, label = TRUE) + NoLegend() + coord_fixed())
dev.off()
# print the colors used
# print(scales::hue_pal()(length(levels(Idents(cns)))))

#
#  Create a gene activity matrix
#

gene_symbol_annotations <- data.table::fread("zcat data/genome/D.melanogaster/Dm6/6.33/gene_names_annotation/fbgn_annotation_ID_fb_2020_02.tsv.gz", sep = "\t", col.names = c("gene_symbol", "organism_abbreviation", "primary_FBgn", "secondary_FBgns", "annotation_ID", "secondary_annotation_IDs"))
gene.ranges$gene_symbol <- gene_symbol_annotations$gene_symbol[match(gene.ranges$gene_id, gene_symbol_annotations$primary_FBgn)]

symbols <- c(
  "elav", # (differentiated neurons)
  "N", # (Notch) (neural progenitor cells, NPC)
  "repo", # (glia)
  # Neurons can be classified depending on the neurotransmitter they express:
  "VAChT", # for cholinergic
  "VGlut", # for glutamatergic
  "Gad1", # for GABAergic
  "twit", # for peptidergic
  "Vmat" # for monoaminergic neurons
)
sel.gene.ranges <- gene.ranges[gene.ranges$gene_symbol %in% symbols]

# Extend coordinates upstream to include the promoter
genebodyandpromoter.coords <- Extend(x = sel.gene.ranges, upstream = 2000, downstream = 0)

# create a gene by cell matrix
gene.activities <- FeatureMatrix(
  fragments = fragment.path,
  features = genebodyandpromoter.coords,
  cells = colnames(cns),
  chunk = 1 # formerly 20
)

# convert rownames from chromsomal coordinates into gene names
gene.key <- genebodyandpromoter.coords$gene_symbol
names(gene.key) <- GRangesToString(grange = genebodyandpromoter.coords)
rownames(gene.activities) <- gene.key[rownames(gene.activities)]

# add the gene activity matrix to the Seurat object as a new assay, and normalize it
cns[['pseudo']] <- CreateAssayObject(counts = gene.activities)
cns <- NormalizeData(
  object = cns,
  assay = 'pseudo',
  normalization.method = 'LogNormalize',
  scale.factor = median(cns$nCount_pseudo)
)

DefaultAssay(cns) <- 'pseudo'

pdf("analysis/scATACseq_gene_activity.pdf", width = 12, height = 8)
print(FeaturePlot(
  object = cns,
  features = symbols,
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3,
  coord.fixed = T
))
dev.off()
