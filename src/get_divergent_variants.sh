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
cut -f1-2 |
awk '{ print $0,"\t1" }' \
> $annotation

bgzip $annotation
tabix -s1 -b2 -e2 $annotation.gz

# Annotate original VCF file with the annotation file and select only divergent sites
header=hdr.txt
echo -e '##INFO=<ID=DIV,Number=1,Type=Integer,Description="Divergent SNP">' >> $header

bcftools annotate -a $annotation.gz \
    -h $header \
    -c CHROM,POS,INFO/DIV \
    $file | 
bcftools filter \
    -i 'DIV=1' \
    > $output

bgzip $output

# Delete the log file produced by vcftools
log=out.log
if test -f "$log"; then
    rm $log
fi

# Delete header file
if test -f "$header"; then
    rm $header
fi
