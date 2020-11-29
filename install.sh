# Create directory where to compile software
sw="sw"

[ ! -d "$sw" ] && mkdir "$sw"

cd sw


# Install bedtools if not installed
command=bedtools

if ! type "$command" > /dev/null
then 
echo "$command is being installed"
wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz
tar -zxvf bedtools-2.29.1.tar.gz
cd bedtools2
make
else
echo "$command is already installed"
fi

# Install bcftools if not installed
command=bcftools

if ! type "$command" > /dev/null
then 
echo "$command is being installed"
git clone git://github.com/samtools/bcftools.git
cd bcftools
make
else
echo "$command is already installed"
fi

# Install vcftools if not installe
command=vcftools

if ! type "$command" > /dev/null
then 
echo "$command is being installed"
git clone https://github.com/vcftools/vcftools.git
cd vcftools
make
else
echo "$command is already installed"
fi