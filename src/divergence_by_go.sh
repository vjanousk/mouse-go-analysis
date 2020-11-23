GO2GENES=$1
GOTERMS=$2
DIVERGENCE=$3
OUTFILE=$4
WD=$5
MINNUMGENES=$6
OUT=$WD/$OUTFILE

< $GO2GENES zcat | 
grep -v "^\!" |
cut -f2,5 |
sort -k2,2 |
join -t $'\t' -12 -22 - <( < $GOTERMS zcat | sort -k2,2 ) \
> $WD/go2genes.txt

sort -k1,1 $DIVERGENCE |
join -t $'\t' -11 -22 - <( <$WD/go2genes.txt sort -k2,2 ) |
sort -k5,7 |
bedtools groupby -i - -g 5-7 -c 4 -o mean,count |
awk -v numg=$MINNUMGENES -F $'\t' '$5 > numg' \
> $OUT





