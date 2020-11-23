# Pipeline to 

QUALITY=50
READDEPTH=10
SOURCEVCF=data/00-source-data/mgp.v5.snps.dbSNP142.clean.X.vcf.gz
SOURCEGENES=data/00-source-data/MGI.gff3.gz

# Prepare VCF file with divergent variants
WORKINGDIR=data/01-divergence
DIVERGENCEVCF=out-vars.vcf

bash src/get_divergent_variants.sh \
$QUALITY \
$READDEPTH \
$WORKINGDIR \
$SOURCEVCF \
$DIVERGENCEVCF

# Prepare CDS database (based on MGI genes)
WORKINGDIR=data/01-divergence
CDSDB=mgi-cds.bed

bash src/make_cds_database.sh \
$SOURCEGENES \
$WORKINGDIR \
$CDSDB

# Calculate per gene divergence
DIVERGENCE=divergence.bed

bash src/calculate_per_gene_divergence.sh \
$WORKINGDIR \
$DIVERGENCEVCF.gz \
$CDSDB \
$DIVERGENCE

# Divergence by GO
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
