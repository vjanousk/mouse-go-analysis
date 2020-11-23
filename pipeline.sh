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

# Prepare VCF file with divergent variants

bash src/get_divergent_variants.sh \
$QUALITY \
$READDEPTH \
$SOURCEVCF \
$ANNOTATION \
$DIVERGENCEVCF

# Calculate per gene divergence

bash src/calculate_per_gene_divergence.sh \
$DIVERGENCEVCF.gz \
$CDSDB \
$DIVERGENCE

# Divergence by GO

bash src/divergence_by_go.sh \
$DIVERGENCE \
$GO_DB \
$MINNUMGENES \
$DIVBYGO