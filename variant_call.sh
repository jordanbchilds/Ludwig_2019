#!/bin/bash
#SBATCH --chdir=/nobackup/proj/clsclmr/Ludwig_2019
#SBATCH -p defq
#SBATCH -A clsclmr
#SBATCH -t 12:00:00
#SBATCH -c 8
#

    ## load modules
module load Python/3.8.6-GCCcore-10.2.0;
module load SAMtools/1.12-GCC-10.2.0;
module load BCFtools/1.10.2-foss-2019b


# only for mitochondrial chromosome chrM
# MAF?

# read bulk ATAC-seq from TF1 cells into array
readarray -t rts < multiQC/group_SRP149534_SRRs.txt;

for rt in "${rts[@]}"; do 
echo $rt;
bcftools mpileup --threads 8 -r J01415.2 -f nuc/GCA_000001405.28_GRCh38.p13_genomic.fna bam/${rt}_sorted.bam | bcftools call --threads 8 -mv -Oz -o vcf/${rt}_mito.vcf.gz;
# -r chrM ???
#bcftools stats -F nuc/GCA_000001405.28_GRCh38.p13_genomic.fna -s - <study.vcf.gz> > <study.vcf.gz.stats>


# filter variants
#bcftools filter -O z -o <study_filtered..vcf.gz> -s LOWQUAL -i'%QUAL>10' <study.vcf.gz>
done 


module purge;
