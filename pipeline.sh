#!/bin/bash

#set -e

# Main pipeline to prepare GO enrichment analysis

quality=$1
readdepth=$2
minnumgenes=$3
sourcevcf=$4
annotation=$5
divergencevcf=$6
cds_db=$7
divergence=$8
go_db=$9
div_go=$10

if [ $# != 11 ]; then
    echo "Incorect number of command line arguments"
    exit
fi

echo "Positional arguments:"
echo "---------------------"

echo "1.  quality: $quality"
echo "2.  readdepth: $readdepth"
echo "3.  minnumgenes: $minnumgenes"
echo "4.  sourcevcf): $sourcevcf"
echo "5.  annotation): $annotation"
echo "6.  divergencevcf): $divergencevcf"
echo "7.  cdsdb): $cds_db"
echo "8.  divergence): $divergence"
echo "9.  go_db): $go_db"
echo "10. divbygo): $div_go"

# Prepare VCF file with divergent variants

echo "Get divergent variants..."

src/get_divergent_variants.sh \
$quality \
$readdepth \
$sourcevcf \
$annotation \
$divergencevcf

# Calculate per gene divergence

echo "Calculate per-gene divergece..."

src/calculate_per_gene_divergence.sh \
$divergencevcf.gz \
$cds_db \
$divergence

# divergence by GO

echo "Calculate divergece by GO..."

src/divergence_by_go.sh \
$divergence \
$go_db \
$minnumgenes \
$div_go

# Run the R script

echo "Prepare ggplot bargraph..."

Rscript src/plot.R

echo "Done."