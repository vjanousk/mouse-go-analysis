#!/bin/bash

# Prepare genes to gene ontology association table
GO2GENES=$1
GOTERMS=$2
GO_DB=$3

< $GO2GENES zcat | 
grep -v "^\!" |
cut -f2,5 |
sort -k2,2 |
join -t $'\t' -12 -22 - <( < $GOTERMS zcat | sort -k2,2 ) \
> $WD/go2genes.txt