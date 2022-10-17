configfile: "config/config.yml"

#
#  Default targets
#

rule all:
  input:
    "data/scATACseq/WTL3CNS/outs/web_summary.html",
    "data/scATACseq/tracks/bulk_coverage.bw",
    "data/scATACseq/tracks/bulk_peaks.bed".
    "data/scATACseq/tracks/clusters.tsv"

#
#  10x Genomics Cell Ranger pipeline
#

rule cellranger_atac:
  output:
    "data/scATACseq/WTL3CNS/outs/web_summary.html",
    "data/scATACseq/WTL3CNS/outs/possorted_bam.bam",
    "data/scATACseq/WTL3CNS/outs/possorted_bam.bam.bai",
    "data/scATACseq/WTL3CNS/outs/peaks.bed"
  shell:
    """
    cd data/scATACseq
    export PATH=/home/ajank/software/cellranger-atac-1.2.0:$PATH
    cellranger-atac count --id=WTL3CNS --fastqs=fastq-renamed --reference=../genome/D.melanogaster/dm6/cellranger-atac/BDGP6.28
    cd ../..
    """

#
#  Peaks, coverage tracks etc.
#

rule bulk_coverage:
  input:
    bam = "data/scATACseq/WTL3CNS/outs/possorted_bam.bam",
    bai = "data/scATACseq/WTL3CNS/outs/possorted_bam.bam.bai"
  output:
    "data/scATACseq/tracks/bulk_coverage.bw"
  threads:
    8
  conda:
    "../../env/deeptools.yaml"
  shell:
    """
    bamCoverage --numberOfProcessors {threads} -b {input.bam} -o {output} --normalizeUsing RPGC --effectiveGenomeSize 130000000
    """

rule bulk_peaks:
  input:
    "data/scATACseq/WTL3CNS/outs/peaks.bed"
  output:
    "data/scATACseq/tracks/bulk_peaks.bed"
  shell:
    """
    cp -p {input} {output}
    """

#
#  Signac analysis
#

rule signac_test:
  input:
    "data/scATACseq/WTL3CNS/outs/filtered_peak_bc_matrix.h5"
  output:
    "data/scATACseq/filtered_tracks/clusters.tsv"
  conda:
    "../../env/r-signac.yaml"
  shell:
    """
    Rscript src/R/scATACseq_pipeline.R
    """
