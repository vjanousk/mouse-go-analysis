#!/bin/bash

DIVERGENCE=$1
GO_DB=$2
MINNUMGENES=$3
OUT=$4

sort -k1,1 $DIVERGENCE |
join -t $'\t' -11 -22 - <( <$GO_DB sort -k2,2 ) |
sort -k5,7 |
bedtools groupby -i - -g 5-7 -c 4 -o mean,count |
awk -v numg=$MINNUMGENES -F $'\t' '$5 > numg' \
> $OUT