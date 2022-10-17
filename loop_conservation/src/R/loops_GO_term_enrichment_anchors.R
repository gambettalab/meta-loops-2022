# ---
# title: "GO anchors"
# author: "Patrycja Rosa"
# date: "25 02 2022"
# output: html_document
# ---

library(ggplot2)
library(clusterProfiler)
library(org.Dm.eg.db)
library(stringr)


#read data

genes_proximal <- read.table("data/loops/genes_proximal_loops_D_mel.tsv", sep = "", quote = "\"'", header = TRUE)
genes_distal <- read.table("data/loops/genes_distal_loops_D_mel.tsv", sep = "", quote = "\"'", header = TRUE)
genes_distal_2TSS <- read.table("data/loops/genes_distal_loops_D_mel_2TSS.tsv", sep = "", quote = "\"'", header = TRUE)


# ALL GENES AS AN UNIVERSE

## Gene ontologies for proximal loops genes
### We are comparing proximal loops genes to the universe of all Drosophila melanogaster genes.

#GO over-representation analysis
proximal_ego <- enrichGO(gene = genes_proximal$gene_id,
                keyType = "ENSEMBL",
                OrgDb = org.Dm.eg.db,
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05)

xlim <- c(0, max(sapply(proximal_ego$GeneRatio, function(str) eval(parse(text=str))), 0.5))

pdf("analysis/plots/loops_GOtermenrichment_proximal.pdf", width = 6, height = 9)
p <- dotplot(proximal_ego, showCategory=20) +
  lims(x = xlim) +
  scale_color_gradient(name = "Adjusted\np\u00advalue", limits = c(0, 0.05)) +
  scale_y_discrete(labels=function(x) str_wrap(x, width=39))
print(p)
dev.off()


## Gene ontologies for distal loops genes
### We are comparing distal loops genes to the universe of all Drosophila melanogaster genes.

#GO over-representation analysis
distal_ego <- enrichGO(gene = genes_distal$gene_id,
                keyType = "ENSEMBL",
                OrgDb = org.Dm.eg.db,
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05)

pdf("analysis/plots/loops_GOtermenrichment_distal.pdf", width = 6, height = 9)
p <- dotplot(distal_ego, showCategory=20) +
  lims(x = xlim) +
  scale_color_gradient(name = "Adjusted\np\u00advalue", limits = c(0, 0.05)) +
  scale_y_discrete(labels=function(x) str_wrap(x, width=39))
print(p)
dev.off()


## Gene ontologies for distal 2TSS loops genes
### We are comparing distal 2TSS loops genes to the universe of all Drosophila melanogaster genes.

#GO over-representation analysis
distal_2TSS_ego <- enrichGO(gene = genes_distal_2TSS$gene_id,
                keyType = "ENSEMBL",
                OrgDb = org.Dm.eg.db,
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05)

pdf("analysis/plots/loops_GOtermenrichment_distal_2TSS.pdf", width = 6, height = 9)
p <- dotplot(distal_2TSS_ego, showCategory=20) +
  lims(x = xlim) +
  scale_color_gradient(name = "Adjusted\np\u00advalue", limits = c(0, 0.05)) +
  scale_y_discrete(labels=function(x) str_wrap(x, width=39))
print(p)
dev.off()
