#!/bin/bash

# Prepare genes to gene ontology association table
go2genes=$1
goterms=$2
go_db=$3

< $go2genes zcat | 
grep -v "^\!" |
cut -f2,5 |
sort -k2,2 |
join -t $'\t' -12 -22 - <( < $goterms zcat | sort -k2,2 ) \
> $go_db