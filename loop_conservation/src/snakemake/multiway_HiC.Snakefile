configfile: "config/config.yml"

#
#  Default targets
#

rule all:
  input:
    lambda wildcards:
      rules.all_multiway.input

rule all_multiway:
  input:
    expand("data/HiC/bam/all_reads/dm6_HiC_WT_L3_CNS_{replicate}.merge.fixmate.rmdup.pairs.gz",
      replicate = ['Rep1', 'Rep2', 'Rep4'])

#
#  Hi-C read haplotype separation: remove duplicates, annotate haplotypes and separate them
#

rule picard_rmdup_dm6_HiC:
  input:
    "data/HiC/bam/merged_reads/dm6_HiC_{dataset}.merge.fixmate.bam"
  output:
    bam = "data/HiC/bam/all_reads/dm6_HiC_{dataset}.merge.fixmate.rmdup.bam",
    metrics = "data/HiC/bam/all_reads/dm6_HiC_{dataset}.merge.fixmate.rmdup.metrics.txt"
  log:
    "data/HiC/bam/all_reads/dm6_HiC_{dataset}.merge.fixmate.rmdup.log"
  params:
    tmpdir = config["paths"]["tmpdir"]
  conda:
    "../../env/picard-slim.yaml"
  shell:
    """
    picard -Xmx32g MarkDuplicates -INPUT {input} -OUTPUT {output.bam} -REMOVE_DUPLICATES true -READ_NAME_REGEX '[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*' -ASSUME_SORT_ORDER queryname -VALIDATION_STRINGENCY SILENT -METRICS_FILE {output.metrics} --MAX_RECORDS_IN_RAM 5000000 -TMP_DIR {params.tmpdir} &> {log}
    """

#
#  Extract multi-way interactions
#

rule pairtools_parse_multiway:
  input:
    "data/HiC/bam/all_reads/{dataset}.merge.fixmate.rmdup.bam"
  output:
    "data/HiC/bam/all_reads/{dataset}.merge.fixmate.rmdup.pairs.gz"
  params:
    genome = config["genome"]["symbol"],
    fasta_genome = config["genome"]["fasta"]
  conda:
    "../../env/pairtools.yaml"
  shell:
    """
    pairtools parse {input} \
      --assembly {params.genome} \
      --chroms-path {params.fasta_genome}.fai \
      --walks-policy all \
      --add-columns pos5,pos3 \
      --drop-sam \
      --no-flip \
      --output {output}
    """
