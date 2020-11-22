WD=$1
VCF=$2
GENES=$3
OUTFILE=$4
VCFPATH=$WD/$VCF
GENESPATH=$WD/$GENES
OUT=$WD/$OUTFILE

< $VCFPATH zcat | 
bedtools intersect -a $GENESPATH -b - -c |
awk '{print $0"\t"$3-$2}' |
bedtools groupby -i - -g 4 -c 5,6 -o sum |
awk '{print $0"\t"$2/$3}' \
> $OUT