#!/bin/bash
#SBATCH --workdir=/nobackup/proj/clsclmr/Ludwig_2019
#SBATCH -p defq
#SBATCH -A clsclmr
#SBATCH -t 48:00:00
#SBATCH -c 8
#

# load modules
module load FastQC/0.11.8-Java-1.8.0_144
module load parallel/20200522-GCCcore-10.2.0


find fastq/SRR*.fastq.gz | parallel --job 8 "fastqc --noextract --outdir multiQC {}"

module purge
