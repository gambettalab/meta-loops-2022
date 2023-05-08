configfile: "config/config.yml"

GENOMES = ["D_mel", "D_vir"]

#
#  Default targets
#

rule all:
  input:
    "data/lastz/ceFinalAnnotated.D_mel.D_vir.Rdata",
    "data/loops/loops_D_mel_annotated.tsv",
    "data/loops/loops_D_vir_annotated.tsv",
    "data/loops/loops_D_mel_anchors.bed",
    "data/loops/loops_D_mel_to_D_vir.Rdata",
    "analysis/plots/loops_Dmel_to_Dvir.pdf",
    "analysis/plots/loops_D_mel_to_D_vir_conserved.pdf",
    "analysis/plots/loops_GOtermenrichment_proximal.pdf",
    "analysis/plots/loops_GOtermenrichment_distal.pdf",
    "analysis/plots/loops_D_mel_to_D_vir_anchor_orthologs.pdf",
    "analysis/plots/multiway_triplet_stats_mirrorcontrol_ext100kb_margin_2kb.pdf"

#
#  Extract soft-masked genomes for evolutionary analysis
#

rule softmask_fasta_D_vir:
  input:
    unmasked = "data/genome/Drosophila_Renschler2019/GSE120751/suppl/GSE120751_Dvir_HiC.fa",
    hardmasked = "data/genome/Drosophila_Renschler2019/GSE120751/suppl/GSE120751_Dvir_HiC_hardmasked.fa"
  output:
    "data/genome/Drosophila_Renschler2019/GSE120751/suppl/GSE120751_Dvir_HiC_softmasked.fa"
  conda:
    "../../env/biopython.yaml"
  shell:
    """
    src/py/softmask_fasta.py {input.unmasked} {input.hardmasked} {output}
    """

rule extract_softmasked_fasta_D_mel:
  input:
    config["genome"]["fasta"] # soft-masked FASTA from UCSC
  output:
    "data/lastz/fasta/D_mel.fa"
  params:
    chromosomes = config["genome"]["chromosomes"]
  conda:
    "../../env/samtools.yaml"
  shell:
    """
    samtools faidx {input} {params.chromosomes} | sed -e "s/^>chr/>/" > {output}
    """

rule extract_softmasked_fasta_D_vir:
  input:
    "data/genome/Drosophila_Renschler2019/GSE120751/suppl/GSE120751_Dvir_HiC_softmasked.fa"
  output:
    "data/lastz/fasta/D_vir.fa"
  params:
    chromosomes = config["genomes"]["D_vir"]["chromosomes"]
  conda:
    "../../env/samtools.yaml"
  shell:
    """
    samtools faidx {input} {params.chromosomes} > {output}
    """

#
#  Annotate loops using gene TSSes
#

rule makeTxDbFromFlyBase:
  input:
    "data/genome/D.melanogaster/dm6/dmel_r6.36_FB2020_05/gtf/dmel-all-r6.36.gtf.gz",
  output:
    "data/genome/D.melanogaster/dm6/dmel_r6.36_FB2020_05/txdb/dmel-all-r6.36.sqlite"
  conda:
    "../../env/r-genomicfeatures.yaml"
  shell:
    """
    Rscript src/R/makeTxDbFromFlyBase.R
    """

rule makeTxDb_Dvir_HiC:
  input:
    "data/genome/Drosophila_Renschler2019/GSE120751/suppl/GSE120751_Dvir_HiC.gtf.gz",
  output:
    "data/genome/Drosophila_Renschler2019/txdb/GSE120751_Dvir_HiC.sqlite"
  conda:
    "../../env/r-genomicfeatures.yaml"
  shell:
    """
    Rscript src/R/makeTxDb_Dvir_HiC.R
    """

rule loops_annotate:
  input:
    "data/loops/loops_D_mel.tsv",
    "data/genome/D.melanogaster/dm6/dmel_r6.36_FB2020_05/txdb/dmel-all-r6.36.sqlite"
  output:
    "data/loops/loops_D_mel_annotated.tsv",
    "data/loops/genes_proximal_loops_D_mel.tsv",
    "data/loops/genes_distal_loops_D_mel.tsv"
  conda:
    "../../env/r-genomicfeatures.yaml"
  shell:
    """
    Rscript src/R/loops_annotate.R
    """

rule loops_annotate_other:
  input:
    "data/loops/loops_D_vir.tsv",
    "data/genome/Drosophila_Renschler2019/txdb/GSE120751_Dvir_HiC.sqlite"
  output:
    "data/loops/loops_D_vir_annotated.tsv"
  conda:
    "../../env/r-genomicfeatures.yaml"
  shell:
    """
    Rscript src/R/loops_annotate_other.R
    """

rule loops_GO_term_enrichment_anchors:
  input:
    "data/loops/genes_proximal_loops_D_mel.tsv",
    "data/loops/genes_distal_loops_D_mel.tsv"
  output:
    "analysis/plots/loops_GOtermenrichment_proximal.pdf",
    "analysis/plots/loops_GOtermenrichment_distal.pdf"
  conda:
    "../../env/r-clusterprofiler-org.dm.eg.db.yaml"
  shell:
    """
    Rscript src/R/loops_GO_term_enrichment_anchors.R
    """

#
#  Extract BED file with loop anchors
#

rule loops_extract_anchors:
  input:
    expand("data/loops/loops_{genome}_annotated.tsv", genome = GENOMES)
  output:
    expand("data/loops/long_range_loops_{genome}.tsv", genome = GENOMES), # for Juicebox
    expand("data/loops/loops_{genome}_anchors.bed", genome = GENOMES)
  conda:
    "../../env/r-genomicranges-ggplot2.yaml"
  shell:
    """
    Rscript src/R/loops_extract_anchors.R
    """

#
#  Identification of conserved noncoding elements using CNEr
#

rule CNEr_make_chains:
  input:
    "data/lastz/fasta/D_mel.fa",
    "data/lastz/fasta/D_vir.fa"
  output:
    "data/lastz/axt/D_mel.D_vir.net.axt",
    "data/lastz/axt/D_vir.D_mel.net.axt"
  conda:
    "../../env/r-cner-genomicranges-ggplot2.yaml"
  shell:
    """
    Rscript src/R/CNEr_make_chains.R
    """

rule CNEr_detect_CNEs:
  input:
    "data/lastz/axt/D_mel.D_vir.net.axt",
    "data/lastz/axt/D_vir.D_mel.net.axt"
  output:
    "data/lastz/ceFinalAnnotated.D_mel.D_vir.Rdata",
    "data/lastz/cneFinalAnnotated.D_mel.D_vir.Rdata"
  conda:
    "../../env/r-cner-genomicranges-ggplot2.yaml"
  shell:
    """
    Rscript src/R/CNEr_detect_CNEs.R
    """

#
#  Match loops between D. melanogaster and the other species (D. virilis)
#

rule CNEr_match_loops:
  input:
    "data/lastz/ceFinalAnnotated.D_mel.D_vir.Rdata",
    "data/lastz/cneFinalAnnotated.D_mel.D_vir.Rdata",
    expand("data/loops/loops_{genome}_annotated.tsv", genome = GENOMES)
  output:
    "data/loops/loops_D_mel_to_D_vir.Rdata",
    "data/loops/loops_D_vir_to_D_mel.Rdata"
  conda:
    "../../env/r-cner-genomicranges-ggplot2.yaml"
  shell:
    """
    Rscript src/R/CNEr_match_loops.R
    """

rule plot_matched_loops:
  input:
    "data/loops/loops_D_mel_to_D_vir.Rdata"
  output:
    "analysis/plots/loops_Dmel_to_Dvir.pdf"
  conda:
    "../../env/r-cner-genomicranges-ggplot2.yaml"
  shell:
    """
    Rscript src/R/plot_matched_loops.R
    """

rule CNEr_loop_permutation_tests:
  input:
    "data/loops/loops_D_mel_to_D_vir.Rdata"
  output:
    "data/loops/loops_D_mel_to_D_vir_permuted.tsv"
  conda:
    "../../env/r-cner-genomicranges-ggplot2.yaml"
  shell:
    """
    Rscript src/R/CNEr_loop_permutation_tests.R
    """

rule plot_matched_loops_stats:
  input:
    "data/loops/loops_D_mel_to_D_vir.Rdata",
    "data/loops/loops_D_mel_to_D_vir_permuted.tsv",
    "data/loops/loops_D_mel_to_D_vir_shifted_stats.tsv"
  output:
    "analysis/plots/loops_D_mel_to_D_vir_conserved.pdf"
  conda:
    "../../env/r-cner-genomicranges-ggplot2.yaml"
  shell:
    """
    Rscript src/R/plot_matched_loops_stats.R
    """

#
#  Test if matched loop anchors overlap orthologs of anchor-associated genes
#

rule loops_fetch_orthologs:
  output:
    "data/loops/genes_loops_D_mel_to_D_vir.tsv",
    "data/loops/genes_loops_D_vir_to_D_mel.tsv"
  conda:
    "../../env/r-genomicfeatures.yaml"
  shell:
    """
    Rscript src/R/loops_fetch_orthologs.R
    """

rule plot_anchor_ortholog_looping:
  input:
    "data/loops/genes_loops_D_mel_to_D_vir.tsv",
    "data/loops/genes_loops_D_vir_to_D_mel.tsv"
  output:
    "analysis/plots/loops_D_mel_to_D_vir_anchor_orthologs.pdf"
  conda:
    "../../env/r-genomicranges-ggplot2.yaml"
  shell:
    """
    Rscript src/R/plot_anchor_ortholog_looping.R
    """

#
#  Plot the frequency of three-way interactions compared to their mirrored controls
#

rule process_multiway:
  input:
    expand("data/HiC/bam/all_reads/dm6_HiC_WT_L3_CNS_{replicate}.merge.fixmate.rmdup.pairs.gz",
      replicate = ['Rep1', 'Rep2', 'Rep4'])
  output:
    "data/HiC/multiway/anchors.txt"
  conda:
    "../../env/r-genomicfeatures.yaml"
  shell:
    """
    Rscript src/R/process_multiway.R
    """

rule plot_multiway:
  input:
    "data/HiC/multiway/anchors.txt"
  output:
    "analysis/plots/multiway_triplet_stats_mirrorcontrol_ext100kb_margin_2kb.pdf"
  conda:
    "../../env/r-genomicranges-ggplot2.yaml"
  shell:
    """
    Rscript src/R/plot_multiway.R
    """

#
#  Extract loop anchor FASTA sequences
#

rule loops_extract_anchor_D_mel_fasta:
  input:
    bed = "data/loops/{prefix}_D_mel_anchors_{margin}.bed",
    fasta_genome = config["genome"]["fasta"]
  output:
    "data/loops/{prefix}_D_mel_anchors_{margin}.fa"
  conda:
    "../../env/bedtools.yaml"
  shell:
    """
    bedtools getfasta -nameOnly -fi {input.fasta_genome} -bed {input.bed} -fo {output}
    """

rule loops_extract_anchor_D_vir_fasta:
  input:
    bed = "data/loops/{prefix}_D_vir_anchors_{margin}.bed",
    fasta_genome = "data/genome/Drosophila_Renschler2019/GSE120751/suppl/GSE120751_Dvir_HiC_softmasked.fa"
  output:
    "data/loops/{prefix}_D_vir_anchors_{margin}.fa"
  conda:
    "../../env/bedtools.yaml"
  shell:
    """
    bedtools getfasta -nameOnly -fi {input.fasta_genome} -bed {input.bed} -fo {output}
    """

rule loops_combine_anchors_fasta:
  input:
    "data/loops/loops_D_mel_anchors_{margin}.fa",
    "data/loops/loops_D_vir_anchors_{margin}.fa",
    "data/loops/loops_D_mel_matched_nonlooping_in_D_vir_anchors_{margin}.fa"
  output:
    "data/loops/loops_all_D_mel_D_vir_anchors_{margin}.fa"
  shell:
    """
    cat {input} > {output}
    """
