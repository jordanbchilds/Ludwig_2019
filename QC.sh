#!/bin/bash
#SBATCH --workdir=/nobackup/proj/clsclmr/Ludwig_2019
#SBATCH -p short
#SBATCH -A clsclmr
#SBATCH -t 00:10:00
#SBATCH -c 8
#

## load modules
module load FastQC/0.11.8-Java-1.8.0_144;
#module load parallel/20200522-GCCcore-10.2.0;
module load MultiQC/1.7-foss-2018b-Python-3.6.6;
#
#
## run fastqc on each SRR*.fastq.gz file
#find fastq/SRR*.fastq.gz | parallel --jobs 8 "fastqc --noextract --outdir multiQC/ {}"
#


 ## Get metadata for all from: https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=1&WebEnv=MCID_60b8a39352def33200839b51&o=acc_s%3Aa

# get file sizes for SRR*.fastq.gz files, and change file name to SRR number
#ls -la fastq/ > fastqgz_list.txt;
#sed -i 's|fastq\/||g' fastqgz_list.txt;
#sed -i 's|.fastq.gz||g' fastqgz_list.txt;

 ## fastqgz_list.txt was combined with metadata in R to create file metadata_ls.csv 


# SRR7246683 temporarily excluded for bad file


# grep all lines with RNA-seq as sequencing technology, take SRR number and put into group_RNA_SRRs.txt
grep 'RNA-Seq' metadata_ls.csv | cut -d ',' -f 2 > group_RNA-seq_SRRs.txt;

mkdir tmp_multiqc;


# interesting subsets of data:
# 1. Grouped sequencing type: RNA-seq vs ATAC-seq vs OTHER (mitosc-seq)
# 2. Platform: NextSeq 500 vs Illumina MiSeq
# 3. Subseries (cell type, sequening type, same library preparation?): 14 Subseries of GSE115218


# read SRRs of group into array
readarray -t group_SRRs < group_RNA-seq_SRRs.txt

# mv each SRR*_fastqc.zip file into tmp_multiqc; run multiqc on all in group
for i in ${group_SRRs[@]}
  do
  #echo ${i}
  cp multiQC/${i}_fastqc.zip tmp_multiqc/${i}_fastqc.zip 
done

multiqc tmp_multiqc/*_fastqc.zip -n group_RNA-Seq_multiQC_report
rm -r tmp_multiqc/


#ls -la fastq/ | grep "sequencing type" | cut -f ? > seq-type-sra_list.txt
#grep "sequencing type" GSE115218.txt | cut " " -f 1 | multiqc > seq-type-multiqc



module purge
