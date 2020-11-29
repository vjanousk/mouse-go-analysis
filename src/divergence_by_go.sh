#!/bin/bash

divergence=$1
go_db=$2
minnumgenes=$3
out=$4

sort -k1,1 $divergence |
join -t $'\t' -11 -22 - <( <$go_db sort -t $'\t' -k2,2 ) |
sort -k5,7 |
bedtools groupby -i - -g 5-7 -c 4 -o mean,count |
awk -v numg=$minnumgenes -F $'\t' '$5 > numg' \
> $out