#!/bin/bash

# Additional packages are used in this pipeline which are already installed as modules on the Newcastle University HPC, and are loaded as modules using SLURM. To run the pipeline on a more general Linux system, these need to be downloaded and installed:

mkdir ./software/
# export ./software/bin/ to path 


#### Samtools ####

if [ -f "software/samtools-1.12/samtools" ]; then 
  echo "samtools is installed";
else
  echo "samtools is not installed. Downloading...";
  cd software/
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
  cd ../../
fi


#### fastqc ####

if [ -f software/bin/fastqc ]; then
  echo "samtools is already installed"
else
  echo "fastqc is not installed" 
  cd ./software
  wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
  unzip fastqc_v0.11.9.zip 
  rm fastqc_v0.11.9.zip 
  chmod 755 FastQC/fastqc
  #chmod +x ./${EXECUTABLE} 
  mv FastQC/* bin/
  export PATH=`pwd`/bin/:$PATH; 
  cd ../
fi


#### multiqc ####

if [ -f software/MultiQC/setup.py ]; then
  echo "multiqc is already installed"
else
  echo "multiqc is not installed. Downloading..."
  cd software
  git clone https://github.com/ewels/MultiQC.git
  python MultiQC/setup.py install
  export PATH=`pwd`/bin/:$PATH;
  cd ../
fi


  ## Bowtie2 ##

if [ -f software/bin/bowtie2]; then
  echo "bowtie2 is already installed"
else
  echo "bowtie2 is not installed" 
  cd ./software
  wget "https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.2/bowtie2-2.3.4.2-linux-x86_64.zip"
  unzip "bowtie2-2.3.4.2-linux-x86_64.zip"
  rm "bowtie2-2.3.4.2-linux-x86_64.zip"
  mv bowtie2-2.3.4.2-linux-x86_64/* bin/
  #chmod +x ./${EXECUTABLE} 
  export PATH=`pwd`/bin/:$PATH; 
  cd ../
fi


  ## Mutserve ##

if [ -f software/bin/mutserve ]; then
  echo "Mutserve is already installed"
else
  echo "Mutserve is not installed. Installing..."
  cd software/bin/
  set -e
  
  NAME="Mutserve"
  VERSION="v2.0.0-rc12"
  GITHUB_USER="seppinho"
  GITHUB_REPO="mutserve"
  EXECUTABLE="mutserve"
  ZIP="mutserve.zip"
  
  INSTALLER_URL=https://github.com/${GITHUB_USER}/${GITHUB_REPO}/releases/download/${VERSION}/${ZIP}
  
  echo "Installing ${NAME} ${VERSION}..."
  
  echo "Downloading ${NAME} from ${INSTALLER_URL}..."
  curl -fL ${INSTALLER_URL} -o ${ZIP}
  
  # execute installer
  unzip ./${ZIP}
  
  # change mod for executables
  chmod +x ./${EXECUTABLE}
  
  # remove installer
  rm ./${ZIP}
  echo ""
  GREEN='\033[0;32m'
  NC='\033[0m'
  echo -e "${GREEN}${NAME} ${VERSION} installation completed. ${NC}"
  echo ""
  cd ../../
fi
