#!/bin/bash
#SBATCH --chdir=/nobackup/proj/clsclmr/Ludwig_2019
#SBATCH -p defq
#SBATCH -A clsclmr
#SBATCH -t 03:00:00
#SBATCH -c 8
#

## load modules
module load FastQC/0.11.8-Java-1.8.0_144;
module load parallel/20200522-GCCcore-10.2.0;
module load MultiQC/1.7-foss-2018b-Python-3.6.6;

mkdir fastQC_results;


  ## run fastqc for specific SRR*.fastq.gz files instead of whole group (forward and reverse reads). eg. SRR7245916 
#find fastq/SRR7245916*.fastq.gz | parallel --jobs 8 "fastqc --noextract --outdir fastQC_results/ {}" ;
#find fastq/SRR72458*.fastq.gz | parallel --jobs 8 "fastqc --noextract --outdir fastQC_results/ {}" ;


  ## See metadata for all from: https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=1&WebEnv=MCID_60b8a39352def33200839b51&o=acc_s%3Aa

# get file sizes for SRR*.fastq.gz files, and change file name to SRR number
#ls -la fastq/ > fastqgz_list.txt;
#sed -i 's|fastq\/||g' fastqgz_list.txt;
#sed -i 's|.fastq.gz||g' fastqgz_list.txt;

  # fastqgz_list.txt was combined with metadata in R to create file metadata_ls.csv (depreciated: replaced with SraRunTable.txt)


# interesting potential subsets of data (see categories.txt):
# 1. Grouped sequencing type: RNA-Seq vs ATAC-Seq vs OTHER (mitosc-seq)
# 2. Platform: NextSeq 500 vs Illumina MiSeq
# 3. Subseries (cell type, sequening type, same library preparation?): 14 Subseries of GSE115218: (SRP149534 SRP149535 SRP149536 SRP149537 SRP149538 SRP149539 SRP149540 SRP149541 SRP149542 SRP149545 SRP156531 SRP156532 SRP168762 SRP168821)

# each line of categories.txt contains a string that distinguishes it as a sugcategory in metadata file: eg. RNA-seq, SRP149534 (one of 14 subseries), MiSeq etc. 

# Read keyword for a group/s of SRRs
readarray -t types < categories.txt

echo "Group/s: ${types[@]}"

# Change to supported UTF-8 variable
locale;
export LANG=en_GB.utf8
export LC_ALL="en_GB.utf8"
#locale;


  ## Adaptor trimming ##
# either trimmomatic or bamclipper?


  ## Run multiqc ##

for j in ${types[@]}
do
  mkdir tmp_multiqc
  
  # extract all lines from data/SraRunTable.txt containing $j from categories.txt line
  #grep ${j} data/SraRunTable.txt | cut -d ',' -f 1 > data/group_${j}_SRRs.txt; 
    
  # read SRRs of group into array
  echo "Number of samples: $(wc -l data/group_${j}_SRRs.txt)"
  readarray -t group_SRRs < data/group_${j}_SRRs.txt
  
  # run fastqc on each SRR in group, move each SRR*_fastqc.zip file into tmp_multiqc/, then run multiqc for all in the group
  for i in ${group_SRRs[@]}
  do
    echo ${i}
    find fastq/${i}_1.fastq.gz | parallel --jobs 8 "fastqc --noextract --outdir fastQC_results/ {}" ;
    find fastq/${i}_2.fastq.gz | parallel --jobs 8 "fastqc --noextract --outdir fastQC_results/ {}" ;
    
    cp fastQC_results/${i}_1_fastqc.zip tmp_multiqc/${i}_1_fastqc.zip 
    cp fastQC_results/${i}_2_fastqc.zip tmp_multiqc/${i}_2_fastqc.zip 
  done

  # run multiqc
  multiqc tmp_multiqc/*_fastqc.zip -f -n group_${j}_multiQC_report -o multiQC/

  rm -r tmp_multiqc/
done

# Revert to system locale environment variables
export LANG=C.UTF-8 ;
export LC_ALL= ;
locale;

module purge
