#!/bin/bash
#SBATCH --workdir=/nobackup/proj/clsclmr/Ludwig_2019
#SBATCH -p defq
#SBATCH -A clsclmr
#SBATCH -t 48:00:00
#SBATCH -c 8
#

module load parallel/20200522-GCCcore-10.2.0

#gzip each fastq file

#gzip -r fastq/SRR724588*.fastq
find fastq/SRR*.fastq | parallel --jobs 8 "gzip -r {}"
#find fastq/ -name *.fastq | parallel --jobs 8

module purge


