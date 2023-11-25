#!/bin/bash

# Script to get reduced VCF with divergent variants

qual=$1
dp=$2
file=$3 # Full path to the source VCF file (mgp.v5.snps.dbSNP142.clean.vcf.gz)
annotation=$4
output=$5 # Path to the final VCF file

# Prepare annotation file
bcftools filter \
    -i'QUAL>'$qual' && INFO/DP>'$dp \
    $file | 
vcftools --vcf - \
    --hardy \
    --stdout |
tail -n+2 |
cut -f1-3 |
tr "/" " " |
awk '$3==1 && $5==1' |
cut -f1-2 \
> $annotation

# Retrieve divergent positions
vcftools --gzvcf $file \
    --positions $annotation \
    --recode \
    --stdout \
    > $output

bgzip $output

# Delete the log file produced by vcftools
log=out.log
if test -f "$log"; then
    rm $log
fi

