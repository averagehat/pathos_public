# if needed . . . 
# export https_proxy=http://127.0.0.1:3128
# export http_proxy=http://127.0.0.1:3128 

# e.g. bash install.sh $PWD/../miniconda

set -x
set -e

if [ $# -ne 1 ]; then
  echo 1>&2 "Usage: $0 <INSTALL_BIN_DIRECTORY>"
  exit 3
fi

# CONDA install
BIN=$1
export PATH=$BIN:$PATH


# conda requirements
conda config --add channels bioconda
conda install -y --file conda-requirements-python.txt 
conda install -y --file conda-requirements-software.txt 

# python packages not present on conda
pip install -r requirements.txt 

# == PriceSeqFilter
PRICESOURCE=PriceSource140408
if [ ! -f ${PRICESOURCE}.tar.gz ]; then
    wget https://derisilab.ucsf.edu/software/price/${PRICESOURCE}.tar.gz 
fi
tar -xvf ${PRICESOURCE}.tar.gz 
cd PriceSource140408 && make && ln -s $PWD/PriceSeqFilter $BIN/PriceSeqFilter && cd ..


# === cd-hit-dup ONLY (not the other cd-hit tools) 
git clone https://github.com/weizhongli/cdhit
cd cdhit/cd-hit-auxtools/ && make && ln -s $PWD/cd-hit-dup  $BIN/cd-hit-dup


# == ncbi BLAST 
BLAST_VER1=2.2.30
BLAST_VER=ncbi-blast-${BLAST_VER1}+ 
BLAST_TGZ=$BLAST_VER-x64-linux.tar.gz
BLAST_URL=https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/$BLAST_VER1/$BLAST_TGZ
if [ ! -f ${BLAST_TGZ}.tar.gz ]; then
   wget $BLAST_URL
fi
tar -xvf $BLAST_TGZ
# wget $BLAST_URL -O- | tar xzvf -
ln -s $PWD/$BLAST_VER/bin/blastdbcmd $BIN/blastdbcmd
ln -s $PWD/$BLAST_VER/bin/blastn $BIN/blastn
ln -s $PWD/$BLAST_VER/bin/blastx $BIN/blastx
ln -s $PWD/$BLAST_VER/bin/blastp $BIN/blastp
ln -s $PWD/$BLAST_VER/bin/update_blastdb.pl $BIN/update_blastdb.pl
ln -s $PWD/$BLAST_VER/bin/makeblastdb $BIN/makeblastdb
ln -s $PWD/$BLAST_VER/bin/blastdbcheck $BIN/blastdbcheck
