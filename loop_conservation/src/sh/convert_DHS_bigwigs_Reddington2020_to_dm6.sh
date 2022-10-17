#!/bin/bash

dir="data/DNase-seq_Reddington2020"

mkdir -p $dir/wigs_dm6

for f in `cd $dir/wigs; ls *.bigwig`; do
  echo $f
  CrossMap.py bigwig data/genome/D.melanogaster/dm3/liftOver/dm3ToDm6.over.chain.gz $dir/wigs/$f $dir/wigs_dm6/${f%.bigwig}
done
