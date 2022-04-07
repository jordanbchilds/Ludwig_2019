#!/bin/bash

# Additional packages are used in this pipeline which are already installed as modules on the Newcastle University HPC, and are loaded as modules using SLURM. To run the pipeline on a more general Linux system, these need to be downloaded and installed:

mkdir ./software/
cd ./software/

# export ./software/bin/ to path 

  ## Samtools ##
wget https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2
tar -xvf samtools-1.12.tar.bz2
rm samtools-1.12.tar.bz2

prefix="`pwd`"
cd ./samtools-1.12
./configure --prefix=$prefix

make
make install

# symlink to ../bin/

export PATH=`pwd`/bin/:$PATH;

#  ## fastq ##
#
#chmod +x ./${EXECUTABLE}
#
#
#  ## Mutserve ##
#
#set -e
#
#NAME="Mutserve"
#VERSION="v2.0.0-rc12"
#GITHUB_USER="seppinho"
#GITHUB_REPO="mutserve"
#EXECUTABLE="mutserve"
#ZIP="mutserve.zip"
#
#INSTALLER_URL=https://github.com/${GITHUB_USER}/${GITHUB_REPO}/releases/download/${VERSION}/${ZIP}
#
#echo "Installing ${NAME} ${VERSION}..."
#
#echo "Downloading ${NAME} from ${INSTALLER_URL}..."
#curl -fL ${INSTALLER_URL} -o ${ZIP}
#
## execute installer
#unzip ./${ZIP}
#
## change mod for executables
#chmod +x ./${EXECUTABLE}
#
## remove installer
#rm ./${ZIP}
#
#echo ""
#GREEN='\033[0;32m'
#NC='\033[0m'
#echo -e "${GREEN}${NAME} ${VERSION} installation completed. ${NC}"
#echo ""
