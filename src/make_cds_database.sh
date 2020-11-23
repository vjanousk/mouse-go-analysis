#!/bin/bash

# Prepare CDS database (based on MGI genes)
GENES=$1 # data/00-source-data/MGI.gff3.gz
CDS_DB=$2

TAB=$'\t'

< $GENES zcat | 
awk '$1=="X" && $2=="NCBI" && $3=="CDS"' | 
sed -E "s/^(.+${TAB}).+gene_id=([^;]+).+$/\1\2/" | 
awk '{ print $1"\t"$4"\t"$5"\t"$9 }' \
> $CDS_DB
