#!/bin/bash

VCF=$1
GENES=$2
OUT=$3

< $VCF zcat | 
bedtools intersect -a $GENES -b - -c |
awk '{print $0"\t"$3-$2}' |
bedtools groupby -i - -g 4 -c 5,6 -o sum |
awk '{print $0"\t"$2/$3}' \
> $OUT