#!/bin/bash
#
#SBATCH --chdir=/nobackup/proj/clsclmr/Ludwig_2019
#SBATCH -p defq
#SBATCH -A clsclmr
#SBATCH -t 48:00:00
#SBATCH -c 8
#

    ## load modules
module load Python/3.8.6-GCCcore-10.2.0;
module load SAMtools/1.12-GCC-10.2.0;
module load parallel/20200522-GCCcore-10.2.0


    ## check if sratools is installed
if [ -f "sratoolkit.2.11.0-ubuntu64/README.md" ]; then 
  echo "SRA-tools is installed";
else
  echo "SRA-tools is not installed. Please see the README.md document to install and configure.";
  exit 1
fi

# Export to shell PATH variable
export PATH=$PATH:`pwd`/sratoolkit.2.11.0-ubuntu64/bin/;


# create directories
mkdir fastq;
mkdir nuc;
mkdir sra;


# To download the samples, you might be tempted to use fastq-dump from sra-tools.
# However, this is slow and unable to resume from broken connection.
# Better to run sra-tools prefetch first, then fastq-dump on result


     ## Use custom python script to download metadata ##

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
python3 parse.py $gse;



    ## Prefetch .sra files ##

#echo "Prefetching all SRRs in ${gse}_sra.txt ...";
#prefetch --option-file "${gse}_sra.txt";
#echo "Done";


    ## Validate	sra files ##
# read bulk ATAC-seq from TF1 cells into array TODO get SRR numbers automatically from GSE
readarray -t rts < group_SRP149534_SRRs.txt 
#readarray -t rts < individual.txt

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


#readarray -t rts < group_SRP149534_SRRs.txt;


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
rm -rf sra;

# zip fastq files
find fastq/SRR*.fastq | parallel --jobs 8 "gzip -r {}"

 ## Download reference genome ##
if [ -f "nuc/hg38.fa" ]; then
  echo "hg38 reference genome already downloaded";
else
  cd nuc/;
  wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz;
  gunzip hg38.fa.gz;
  cd ..;
fi

module purge;

