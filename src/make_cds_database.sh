GENES=$1 # data/00-source-data/MGI.gff3.gz
WD=$2
OUTFILE=$3
OUT=$WD/$OUTFILE

TAB=$'\t'

< $GENES zcat | 
awk '$1=="X" && $2=="NCBI" && $3=="CDS"' | 
sed -E "s/^(.+${TAB}).+gene_id=([^;]+).+$/\1\2/" | 
awk '{ print $1"\t"$4"\t"$5"\t"$9 }' \
> $OUT