
# Meta-loop calling


## Dependencies

* [cooler](https://github.com/open2c/cooler) (v0.8.11).

* [R](https://www.R-project.org) (v4.2.0) with the following packages:

  * [optparse](https://CRAN.R-project.org/package=optparse) (v1.7.1).
  * [data.table](https://CRAN.R-project.org/package=data.table) (v1.14.2).
  * [igraph](https://CRAN.R-project.org/package=igraph) (v1.3.5).
  * [mlpack](https://CRAN.R-project.org/package=mlpack) (v4.0.0).
  * [EBImage](http://bioconductor.org/packages/EBImage) (v4.38.0).
  * [rhdf5](https://bioconductor.org/packages/rhdf5/) (v2.40.0).

All computations were done on Linux (CentOS 7). Version number used for the analysis are specified above (parenthesis).


## Usage

### Data preparation

#### Drosophila melanogaster

Download Hi-C data in bedpe format for all WT third instar larval CNS replicates from [Kaushal et al. 2021](http://dx.doi.org/10.1038/s41467-021-21366-2) (Gene Expression Omnibus, accession code [GSE146750](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146750)):

```
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4405nnn/GSM4405399/suppl/GSM4405399_brain_WT_1.bedpe.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4405nnn/GSM4405400/suppl/GSM4405400_brain_WT_2.bedpe.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4405nnn/GSM4405401/suppl/GSM4405401_brain_WT_4.bedpe.gz
```

Generate cooler files binned at 1kb (ignoring chrM and chrY)

```
zgrep  -v -e chrM -e chrY GSM4405399_brain_WT_1.bedpe.gz | awk -F '\t' 'BEGIN{OFS="\t"}{for(i=1;i<=$7;i++)print $1,int(($2+$3)/2),$4,int(($5+$6)/2)}' | cooler cload pairs --chrom1 1 --pos1 2 --chrom2 3 --pos2 4 --assembly dm6 data/dm6.chrom.sizes.tsv:1000 - GSM4405399_brain_WT_1.cool
zgrep  -v -e chrM -e chrY GSM4405400_brain_WT_2.bedpe.gz | awk -F '\t' 'BEGIN{OFS="\t"}{for(i=1;i<=$7;i++)print $1,int(($2+$3)/2),$4,int(($5+$6)/2)}' | cooler cload pairs --chrom1 1 --pos1 2 --chrom2 3 --pos2 4 --assembly dm6 data/dm6.chrom.sizes.tsv:1000 - GSM4405400_brain_WT_2.cool
zgrep  -v -e chrM -e chrY GSM4405401_brain_WT_4.bedpe.gz | awk -F '\t' 'BEGIN{OFS="\t"}{for(i=1;i<=$7;i++)print $1,int(($2+$3)/2),$4,int(($5+$6)/2)}' | cooler cload pairs --chrom1 1 --pos1 2 --chrom2 3 --pos2 4 --assembly dm6 data/dm6.chrom.sizes.tsv:1000 - GSM4405401_brain_WT_4.cool
```

Merge replicates

```
cooler merge Kaushal2021_WT_larval_CNS.cool GSM4405399_brain_WT_1.cool GSM4405400_brain_WT_2.cool GSM4405401_brain_WT_4.cool
```

Create multi-resolution cooler file (normalized with iterative correction)

```
cooler zoomify --balance --out Kaushal2021_WT_larval_CNS.mcool Kaushal2021_WT_larval_CNS.cool
```

Clean temporary files

```
rm GSM4405399_brain_WT_1.bedpe.gz
rm GSM4405400_brain_WT_2.bedpe.gz
rm GSM4405401_brain_WT_4.bedpe.gz
rm GSM4405399_brain_WT_1.cool
rm GSM4405400_brain_WT_2.cool
rm GSM4405401_brain_WT_4.cool
rm Kaushal2021_WT_larval_CNS.cool
```


#### Drosophila virilis

Download Hi-C data (pooled replicates) in mcool format for drosophila virilis adult CNS (Gene Expression Omnibus, accession code [GSE214704](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE214704)) and save as `Dvir_WT_adult_CNS_HiC.mcool`.


### Meta-loop calling

Use `meta_loops.R` to call meta-loops on drosophila melanogaster WT third instar larval CNS replicates from Kaushal et al. 2021

```
./meta_loops.R --output=Dmel-meta-loops.tsv --resolution=2000 --chrs=chr2L,chr2R,chr3L,chr3R,chr4,chrX Kaushal2021_WT_larval_CNS.mcool
```

The list of meta-loops will be saved in the file `Dmel-meta-loops.tsv`. This file is in tab separated format with header in first row, one row per meta-loop and columns chr1, start1, end1 (first anchor), chr2, start2, end2 (second anchor). 


For drosophila virilis adult CNS

```
./meta_loops.R --output=Dvir-meta-loops.tsv --resolution=4000 --chrs=2L,2R,3L,3R,4,X Dvir_WT_adult_CNS_HiC.mcool
```

The list of meta-loops will be saved in the file `Dvir-meta-loops.tsv`.

