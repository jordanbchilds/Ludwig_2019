#!/bin/bash

# Additional packages are used in this pipeline which are already installed as modules on the Newcastle University HPC, and are loaded as modules using SLURM. To run the pipeline on a more general Linux system, these need to be downloaded and installed:

mkdir ./software/
# export ./software/bin/ to path 


#### Samtools ####

if [ -f "software/samtools-1.12/samtools" ]; then 
  echo "samtools is installed";
else
  echo "samtools is not installed. Downloading...";
  wget https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2
  tar -xvf samtools-1.12.tar.bz2
  rm samtools-1.12.tar.bz2
  
  prefix="`pwd`"
  cd ./software/samtools-1.12
  ./configure --prefix=$prefix
  
  make
  make install
  # symlink to ../bin/
  
  export PATH=`pwd`/bin/:$PATH;
  cd ../../
fi


#### fastqc ####

if [ -f software/bin/fastqc/fastqc ]; then
  echo "samtools is already installed"
else
  echo "fastqc is not installed" 
  cd ./software/bin
  wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
  unzip fastqc_v0.11.9.zip 
  rm fastqc_v0.11.9.zip 
  mv FastQC/* .
  
  chmod 755 fastqc
  #chmod +x ./${EXECUTABLE} 
  export PATH=`pwd`/bin/:$PATH;
  
  cd ../../
fi


#### multiqc ####

#if [-f software/m]


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
