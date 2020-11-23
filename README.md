# Mouse Gene Ontology enrichment analysis pipeline

This repository represents a bioinformatic pipeline to carry Gene Ontology enrichment analysis for low and high divergence genes among two house mouse subspecies. The pipeline is used as an exercise at the course on UNIX and work with genomic data.

## Overview

SNP variants for two mouse strains (PWD/PhJ, WSB/EiJ) were downloaded from the Mouse Genome Project FTP site (ftp://ftp-mouse.sanger.ac.uk/). PWD/PhJ and WSB/EiJ represent Mus musculus musculus and Mus musculus domesticus subspecies, respectively.

The aim is to identify genes with high relative divergence between the two strains and carry Gene Ontology enrichment analysis for genes according to the divergence.

## Description of the pipeline

1. Selecting SNPs that are divergent between the two strains

Other criteria used for selection is the PHRED quality and read depth (DP). Divergent SNPs are identified using Fst function built in the `vcftools`. SNPs are considered to be divergent when Fst equals 1.

```bash
QUALITY=50
READDEPTH=10
SOURCEVCF=data/00-source-data/mgp.v5.snps.dbSNP142.clean.X.vcf.gz

WORKINGDIR=data/01-divergence
DIVERGENCEVCF=out-vars.vcf

bash src/get_divergent_variants.sh \
$QUALITY \
$READDEPTH \
$WORKINGDIR \
$SOURCEVCF \
$DIVERGENCEVCF
```

2. Prepare CDS database

`MGI.gff3.gz` represents a full report containing detailed information on genes, mRNAs, exons and CDS. For the divergence analysis only CDS are needed. CDS database is prepared in this step and `gff3` is converted to `bed` to work more easily with the CDS data.

```bash
SOURCEGENES=data/00-source-data/MGI.gff3.gz
WORKINGDIR=data/01-divergence
CDSDB=mgi-cds.bed

bash src/make_cds_database.sh \
$SOURCEGENES \
$WORKINGDIR \
$CDSDB
```

3. Calculate the per gene divergence

Once the list of divergent SNPs between the two strains and the CDS database are created, the divergence per gene can be calculated. Combination of `bedtools` tools and `awk` commands is used to find SNPs overlapping CDS parts of the genes and calculate sums and relative divergence by genes.

```bash
DIVERGENCE=divergence.bed

bash src/calculate_per_gene_divergence.sh \
$WORKINGDIR \
$DIVERGENCEVCF.gz \
$CDSDB \
$DIVERGENCE
```

4. Calculate the average relative divergence by Gene Ontology category

Per-gene relative divergences are used to calculate the average relative divergence for individual GO terms. Combinatino of the built-in UNIX `join` and `sort` commands is used along with `groupby` that is part of the `bedtools` tools suite. Dataset representing association between genes and GO terms provided by Mouse Genome Informatics (http://www.informatics.jax.org) and Gene Ontology Consortium (http://geneontology.org) is joined to dataset on with gene relative divergences. The average for every GO term is then calculated omitting low prevalence GO terms.

```bash
GO2GENES=data/00-source-data/gene_association.mgi.gz
GOTERMS=data/00-source-data/go_terms.mgi.gz
DIVERGENCE=data/01-divergence/divergence.bed
DIVBYGO=divergence_by_go.txt
WORKINGDIR=data/02-go
MINNUMGENES=9

bash src/divergence_by_go.sh \
$GO2GENES \
$GOTERMS \
$DIVERGENCE \
$DIVBYGO \
$WORKINGDIR \
$MINNUMGENES
```