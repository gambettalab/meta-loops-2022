#!/bin/bash

# Run the pipeline locally, using 32 cores
snakemake --snakefile src/snakemake/loops.Snakefile --use-conda -j 32 -k
