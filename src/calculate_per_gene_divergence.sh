#!/bin/bash

vcf=$1
genes=$2
out=$3

< $vcf zcat | 
bedtools intersect -a $genes -b - -c |
awk '{print $0"\t"$3-$2}' |
bedtools groupby -i - -g 4 -c 5,6 -o sum |
awk '{print $0"\t"$2/$3}' \
> $out