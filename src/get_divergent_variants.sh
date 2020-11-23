# Script to get reduced VCF with divergent variants

QUAL=$1
DP=$2
WD=$3 # Working directory (where the output should be stored)
FILE=$4 # Full path to the source VCF file (mgp.v5.snps.dbSNP142.clean.vcf.gz)
OUTFILENAME=$5 # Name of the final VCF file
OUTPUT=$WD/$OUTFILENAME # Path to the final VCF file

SNPLIST=$WD/snplist.tab

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
> $SNPLIST

bgzip $SNPLIST
tabix -s1 -b2 -e2 $SNPLIST.gz

# Annotate original VCF file with the annotation file and select only divergent sites
HEADER=$WD/hdr.txt
echo -e '##INFO=<ID=DIV,Number=1,Type=Integer,Description="Divergent SNP">' >> $HEADER

bcftools annotate -a $SNPLIST.gz \
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
