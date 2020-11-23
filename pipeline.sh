#!/bin/bash

# Main pipeline to prepare GO enrichment analysis

QUALITY=$1
READDEPTH=$2
MINNUMGENES=$3
SOURCEVCF=$4
ANNOTATION=$5
DIVERGENCEVCF=$6
CDSDB=$7
DIVERGENCE=$8
GO2GENES=$9
GOTERMS=$10
GO_DB=$11
DIVBYGO=$12

if [ $# !=13 0 ]; then
    echo "Incorect number of command line arguments"
    exit
fi

echo "Positional argument 1 (=QUALITY): $QUALITY"
echo "Positional argument 2 (=READDEPTH): $READDEPTH"
echo "Positional argument 3 (=MINNUMGENES): $MINNUMGENES"
echo "Positional argument 4 (=SOURCEVCF): $SOURCEVCF"
echo "Positional argument 5 (=ANNOTATION): $ANNOTATION"
echo "Positional argument 6 (=DIVERGENCEVCF): $DIVERGENCEVCF"
echo "Positional argument 7 (=CDSDB): $CDSDB"
echo "Positional argument 8 (=DIVERGENCE): $DIVERGENCE"
echo "Positional argument 9 (=GO2GENES): $GO2GENES"
echo "Positional argument 10 (=GOTERMS): $GOTERMS"
echo "Positional argument 11 (=GO_DB): $GO_DB"
echo "Positional argument 12 (=DIVBYGO): $DIVBYGO"


# Prepare VCF file with divergent variants

echo "Get divergent variants..."

bash src/get_divergent_variants.sh \
$QUALITY \
$READDEPTH \
$SOURCEVCF \
$ANNOTATION \
$DIVERGENCEVCF

# Calculate per gene divergence

echo "Calculate per-gene divergece..."

bash src/calculate_per_gene_divergence.sh \
$DIVERGENCEVCF.gz \
$CDSDB \
$DIVERGENCE

# Divergence by GO

echo "Calculate divergece by GO..."

bash src/divergence_by_go.sh \
$DIVERGENCE \
$GO_DB \
$MINNUMGENES \
$DIVBYGO

# Run the R script

echo "Prepare ggplot bargraph..."

bash src/plot.R

echo "Done."