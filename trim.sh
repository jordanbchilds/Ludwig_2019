#!/bin/bash
#SBATCH --workdir=/nobackup/proj/clsclmr/Ludwig_2019
#SBATCH -p defq
#SBATCH -A clsclmr
#SBATCH -t 12:00:00
#SBATCH -c 8
#

module load cutadapt/1.18-foss-2018b-Python-3.6.6;

readarray -t types < categories.txt;


# trim reads for low quality. REMEMBER: change hisat2 input to _trimmed.fastq for alignment
for i in ${types[@]}; do

## loop so SRRs are read


echo "trimming ${i}";
cutadapt -j 0 -q 30 -a TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -A GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -o fastq/${i}_2_trimmed.fastq.gz -p fastq/${i}_2_trimmed.fastq.gz fastq/${i}_1.fastq.gz fastq/${i}_2.fastq.gz;
# ^ NOT reverse complement of adaptors

done


## Split into nuclear sequences and mitochondrial sequences

