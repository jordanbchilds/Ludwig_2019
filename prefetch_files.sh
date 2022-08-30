#!/bin/bash
#
#SBATCH --chdir=/nobackup/proj/clsclmr/Ludwig_2019
#SBATCH -p defq
#SBATCH -A clsclmr
#SBATCH -t 48:00:00
#SBATCH -c 8
#

  ## Parse Arguments ##
 
 POSITIONAL_ARGS=()
 
while [[ $# -gt 0 ]]; do
  case $1 in
 -g|--group-name)
    g="$2"
    echo "Group name: ${g}"
    shift # past argument
    shift # past value 
  ;;

  -h|--help)
    echo "-g|--group-name       base name for data/group_(group_name)_SRRs.txt, eg. SRP149534"
    shift # past argument
    exit 1
  ;;

  -*|--*)
    echo "Unknown option ${1}. Valid aruguments:"
    echo "-g|--group-name       base name for data/group_(group_name)_SRRs.txt, eg. SRP149534"
    exit 1
  ;;

  *)
    POSITIONAL_ARGS+=("$1") # save positional arg
    shift # past argument
  ;;

  esac
done


    ## load modules
module load Python/3.8.6-GCCcore-10.2.0;
module load SAMtools/1.12-GCC-10.2.0;
module load parallel/20200522-GCCcore-10.2.0


    ## check if sratools is installed
if [ -f "software/sratoolkit.2.11.0-ubuntu64/bin/fastq-load" ]; then 
  echo "SRA-tools is installed";
else
  echo "SRA-tools is not installed. Please see the README.md document to install and configure.";
  exit 1
fi

# Export to shell PATH variable
export PATH=`pwd`/software/bin:$PATH:`pwd`/software/sratoolkit.2.11.0-ubuntu64/bin/;


# create directories
mkdir fastq
mkdir nuc
mkdir sra
mkdir data

# To download the samples, you might be tempted to use fastq-dump from sra-tools.
# However, this is slow and unable to resume from broken connection.
# Better to run sra-tools prefetch first, then fastq-dump on result


     ## Use custom python script to download metadata for .sra downloads ##

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115218
# https://doi.org/10.1016/j.cell.2019.01.022


# superseries (with all Ludwig data, including RNA-seq, scATAC, scRNA, difference cell lines etc...): 
#gse='GSE115218';
# subseries with bulk ATAC-seq TF1 cells:
gse='GSE115208';

# parse.py gets list of SRA sequence names from GSE series of "Human lineage tracing enabled by mitochondrial mutations and single cell genomics"

# install GEOparse python3 module. pip3 install _ won't install if module is already installed.
pip3 install GEOparse;

# run parse.py
cd data/
python3 ../parse.py $gse;
cd ../


    ## Prefetch .sra files (depreciated) ##
## Use $GSE to select runs:
#echo "Prefetching all SRRs in ${gse}_sra.txt ...";
#prefetch --option-file "data/${gse}_sra.txt";
#echo "Done";


    ## Prefetch and validate sra files ##

# read bulk ATAC-seq from TF1 cells into array 

# Get SRR numbers from metadata table SraRunTable.txt: downloaded from https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA474183&o=acc_s%3Aa
# Select which group of runs to analyse with a keyword in categories.txt 
# each line of categories.txt contains a string that distinguishes it as a sugcategory in metadata file: eg. RNA-seq, SRP149534 (one of 14 subseries), MiSeq 

# interesting keywords for subsets of data (see categories.txt):
# 1. Grouped sequencing type: RNA-Seq vs ATAC-Seq vs OTHER (mitosc-seq)
# 2. Platform: NextSeq 500 vs Illumina MiSeq
# 3. Subseries (cell type, sequening type, same library preparation?): 14 Subseries of GSE115218: (SRP149534 SRP149535 SRP149536 SRP149537 SRP149538 SRP149539 SRP149540 SRP149541 SRP149542 SRP149545 SRP156531 SRP156532 SRP168762 SRP168821)

# The bulk-ATAC-seq of TF1 cells all have the subseries ID: SRP149534
#readarray -t types < categories.txt
types=("$g")
echo "Group names: $types"

for j in ${types[@]}
do
  echo "Group name: $j"
  grep $j data/SraRunTable.txt | cut -d ',' -f 1 > data/group_${j}_SRRs.txt; 
  # table is comma separated but has some commas within strings in one "cell": delete everything before "SRX", then everything after the first comma to get the SRX number for prefetch
  grep $j data/SraRunTable.txt | sed s/.*SRX/SRX/g | sed s/\,.*//g > data/group_${j}_SRXs.txt;  
 
  # prefetch using SRX numbers TODO don't prefetch if fastq file present
  prefetch --option-file "data/group_${j}_SRXs.txt";
 

  # for ATAC-seq of TF1 cells:
  readarray -t rts < data/group_${j}_SRRs.txt 
  #readarray -t rts < data/individual.txt
  
  # validate each prefetched file and output any missing or incomplete to 'failed_to_prefetch.txt'
  for i in "${rts[@]}"
  do
  
    echo $i
    vdb-validate sra/sra/${i}.sra &> sra/${i}_validation.txt;
    
    if grep -q 'err' sra/${i}_validation.txt; 
    then
    echo ${i} >> failed_to_prefetch.txt
    prefetch ${i}
    fi
    
    if grep -q "could not be found" sra/${i}_validation.txt;
    then
    echo ${i} >> failed_to_prefetch.txt
    prefetch ${i}
    fi
    
  done
  
  
      ## Convert prefetched .sra files to fasta format ##
  
  # loop to find if .sra file has been dumped (converted to fastq.gz) and if not add file to list
  for i in "${rts[@]}";
  do
  if test -f "fastq/${i}_1.fastq.gz";
  then
    echo "${i}_1.fastq.gz file exists";
  else
    echo ${i} >> dump_list.txt
    echo "${i}.sra added to dump_list.txt";
  fi
  done
  
  cat dump_list.txt | parallel --jobs 8 "fastq-dump --split-files --gzip --outdir fastq/ sra/sra/{}.sra";
  rm dump_list.txt;
  
  # Delete large sra files after dumping
  #rm -rf sra;

done

# zip fastq files
find fastq/SRR*.fastq | parallel --jobs 8 "gzip -r {}"

 ## Download reference genome ##
if [ -f "nuc/hg38.fa" ]; then
  echo "hg38 reference genome already downloaded";
else
  echo "Downloading hg38 reference genome..."
  cd nuc/;
  wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz;
  gunzip hg38.fa.gz;
  cd ..;
fi

module purge;

