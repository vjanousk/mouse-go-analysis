#!/bin/bash

# Script to get reduced VCF with divergent variants

QUAL=$1
DP=$2
FILE=$3 # Full path to the source VCF file (mgp.v5.snps.dbSNP142.clean.vcf.gz)
ANNOTATION=$4
OUTPUT=$5 # Path to the final VCF file

# Prepare annotation file
bcftools filter \
-i'QUAL>'$QUAL' && INFO/DP>'$DP \
$FILE | 
vcftools --vcf - \
--hardy \
--stdout |
tail -n+2 |
cut -f1-3 |
tr "/" " " |
awk '$3==1 && $5==1' |
cut -f1-2 |
awk '{ print $0,"\t1" }' \
> $ANNOTATION

bgzip $ANNOTATION
tabix -s1 -b2 -e2 $ANNOTATION.gz

# Annotate original VCF file with the annotation file and select only divergent sites
HEADER=hdr.txt
echo -e '##INFO=<ID=DIV,Number=1,Type=Integer,Description="Divergent SNP">' >> $HEADER

bcftools annotate -a $ANNOTATION.gz \
-h $HEADER \
-c CHROM,POS,INFO/DIV \
$FILE | 
bcftools filter \
-i 'DIV=1' \
> $OUTPUT

bgzip $OUTPUT

# Delete the log file produced by vcftools
LOG=out.log
if test -f "$LOG"; then
    rm $LOG
fi

# Delete header file
if test -f "$HEADER"; then
    rm $HEADER
fi
