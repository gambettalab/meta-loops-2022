configfile: "config/config.yml"

#
#  Default targets
#

rule all:
  input:
    "data/lastz/ceFinalAnnotated.D_mel.D_vir.Rdata",
    "data/loops/loops_D_mel_annotated.tsv",
    "data/loops/loops_D_mel_anchors.bed",
    "data/loops/loops_D_mel_to_D_vir.Rdata"

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

rule loops_annotate:
  input:
    "data/loops/loops_D_mel.tsv"
  output:
    "data/loops/loops_D_mel_annotated.tsv"
  conda:
    "../../env/r-genomicfeatures.yaml"
  shell:
    """
    Rscript src/R/loops_annotate.R
    """

#
#  Extract BED file with loop anchors
#

rule loops_extract_anchors:
  input:
    "data/loops/loops_D_mel_annotated.tsv",
  output:
    "data/loops/long_range_loops_D_mel.tsv", # for Juicebox
    "data/loops/loops_D_mel_anchors.bed"
  conda:
    "../../env/r-genomicfeatures.yaml"
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
    "data/loops/loops_D_mel_annotated.tsv"
  output:
    "data/loops/loops_D_mel_to_D_vir.Rdata",
    "data/loops/loops_D_vir_to_D_mel.Rdata"
  conda:
    "../../env/r-cner-genomicranges-ggplot2.yaml"
  shell:
    """
    Rscript src/R/CNEr_match_loops.R
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
