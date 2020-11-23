# Create directory where to compile software
SW=.sw
if test -d "$SW"; then
    mkdir $SW
fi

cd .sw

# Install bedtools if not installed
COMMAND=bedtools

if ! type "$COMMAND" > /dev/null
then 
echo "$COMMAND will be installed"
wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz
tar -zxvf bedtools-2.29.1.tar.gz
cd bedtools2
make
fi

# Install bcftools if not installed
COMMAND=bcftools

if ! type "$COMMAND" > /dev/null
then 
echo "$COMMAND will be installed"
git clone git://github.com/samtools/bcftools.git
cd bcftools
make
fi

# Install vcftools if not installe
COMMAND=vcftools

if ! type "$COMMAND" > /dev/null
then 
echo "$COMMAND will be installed"
git clone https://github.com/vcftools/vcftools.git
cd vcftools
make
fi