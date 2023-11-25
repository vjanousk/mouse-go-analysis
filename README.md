# Mouse Gene Ontology enrichment analysis pipeline

This repository represents a bioinformatic pipeline to carry out gene-set enrichment analysis for genes having low and high genomic divergence between two house mouse subspecies. The pipeline is used as an exercise at the course on Unix and work with genomic data.

## Overview

SNP variants for two mouse strains (PWD/PhJ, WSB/EiJ) were downloaded from the Mouse Genome Project FTP site [Mouse Genome Project](https://www.sanger.ac.uk/data/mouse-genomes-project/). PWD/PhJ and WSB/EiJ represent *Mus musculus musculus* and *Mus musculus domesticus* subspecies, respectively.

The aim is to identify genes with high relative divergence between the two strains and carry Gene Ontology enrichment analysis for genes according to the divergence.

## Run the pipeline

#### Test & install required software

If `bedtools`, `vcftools` or `bcftools` are not installed, the script below will install these tools.

```bash
chmod +x install.sh

./install.sh
```

#### Prepare data directories

```bash
mkdir -p data/00-source-data data/01-divergence data/02-go
```

#### Define variables

Several types of variables defined. Filtering parameters provide thresholds on filtering quality and number of genes used at different stages of the pipeline.

```bash
# Filtering parameters
quality=50
readdepth=10
minnumgenes=30

# Working directories (removed from the git in .gitignore)
wd_source=data/00-source-data
wd_divergence=data/01-divergence
wd_go=data/02-go

# Source files
sourcevcf=$wd_source/mgp.v5.snps.dbSNP142.clean.X.vcf.gz
sourcegenes=$wd_source/MGI.gff3.gz
go2genes=$wd_source/gene_association.mgi.gz
goterms=$wd_source/go_terms.mgi.gz

# Helping datasets
cds_db=$wd_source/mgi-cds.bed
go_db=$wd_source/go2genes.txt

# Output files
annotation=$wd_divergence/annotation.tab
divergencevcf=$wd_divergence/out-vars.vcf
divergence=$wd_divergence/divergence.bed
div_go=$wd_go/divergence_by_go.txt
```

#### Make scripts executable

```bash
chmod +x src/make_cds_database.sh src/make_go_database.sh

chmod +x src/calculate_per_gene_divergence.sh src/divergence_by_go.sh src/get_divergent_variants.sh

chmod +x pipeline.sh
```

#### Prepare data

```bash
# Prepare CDS database
src/make_cds_database.sh $sourcegenes $cds_db

# Prepare GO database
src/make_go_database.sh $go2genes $goterms $go_db
```

#### Run the pipeline

```bash
./pipeline.sh \
$quality \
$readdepth \
$minnumgenes \
$sourcevcf \
$annotation \
$divergencevcf \
$cds_db \
$divergence \
$go_db \
$div_go
```

#### Resulting ggplot graph

![results](results/go-enrichment.jpg)

## Run the pipeline step-by-step

#### 1. Prepare CDS & GO databases

`MGI.gff3.gz` represents a full report containing detailed information on genes, mRNAs, exons and CDS. For the divergence analysis only CDS are needed. CDS database is prepared in this step and `gff3` is converted to `bed` to work more easily with the CDS data.

```bash
src/make_cds_database.sh $sourcegenes $cds_db
```

`go_terms.mgi.gz` and `gene_association.mgi.gz` represents GO terms and association between genes and GO terms IDs provided by Mouse Genome Informatics [Mouse Genome Informatics](http://www.informatics.jax.org) and Gene Ontology Consortium [Gene Ontology](http://geneontology.org). In the command below joined dataset of list of genes with GO term enrichment is prepared.

```bash
src/make_go_database.sh $go2genes $goterms $go_db
```

#### 2. Selecting SNPs that are divergent between the two strains

Other criteria used for selection is the PHRED quality and read depth (DP). Divergent SNPs are identified using Fst function built in the `vcftools`. SNPs are considered to be divergent when Fst equals 1.

```bash
src/get_divergent_variants.sh \
$quality \
$readdepth \
$sourcevcf \
$annotation \
$divergencevcf
```

#### 3. Calculate the per gene divergence

Once the list of divergent SNPs between the two strains and the CDS database are created, the divergence per gene can be calculated. Combination of `bedtools` tools and `awk` commands is used to find SNPs overlapping CDS parts of the genes and calculate sums and relative divergence by genes.

```bash
src/calculate_per_gene_divergence.sh \
$divergencevcf.gz \
$cds_db \
$divergence
```

#### 4. Calculate the average relative divergence by Gene Ontology category

Per-gene relative divergences are used to calculate the average relative divergence for individual GO terms. Combinatino of the built-in Unix `join` and `sort` commands is used along with `groupby` that is part of the `bedtools` tools suite. GO dataset is joined to dataset on with gene relative divergences. The average for every GO term is then calculated omitting low prevalence GO terms.

```bash
src/divergence_by_go.sh \
$divergence \
$go_db \
$minnumgenes \
$div_go
```

#### 5. Prepare a barplot showing results of the GO enrichment analysis

To plot the results of the GO enrichment analysis `Rscript` is used. Library `ggplot2` is the most suitable tool to provide fast and efficient plot.

```bash
Rscript src/plot.R
```
