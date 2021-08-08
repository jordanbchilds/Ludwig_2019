#!/bin/bash
#
#SBATCH --chdir=/nobackup/proj/clsclmr/Ludwig_2019
#SBATCH -p defq
#SBATCH -A clsclmr
#SBATCH -t 24:00:00
#SBATCH -c 8
#

module load SAMtools/1.12-GCC-10.2.0;

mkdir coverages/;


# ls -1 prints each file on a new line
ls -1 -d bam/* | grep -v "bai" > ls_bam_files.txt;


  ## SAMtools depth ##
# Calculate depth for all positions (-a), comment line of column names (-H), base (-q) and mapping quality (-Q) greater than 20: (based on default settings for mutserve variant caller), region chrM (-r), remove depth limit (-d 0) 
#samtools depth -f ls_bam_files.txt -a -H -q 20 -Q 20 -r chrM -d 0 -o coverages/depths_qfilt.txt

# no specified quality limits: defaults?
samtools depth -f ls_bam_files.txt -a -H -r chrM -d 0 -o coverages/depths.txt;



  ## SAMtools coverage ##
# calculate mean depth, SD, no reads aligned, mean base quality, mean mapping quality, proportion of bases with depth <1 (why 1??) and output tab separated file. 
# list of bam files (-b), min base quality (-q), min mapping quality (-Q), mitochondrial chromosome (-r chrM).

#echo "mean_coverage_qfilt";
samtools coverage -b ls_bam_files.txt -q 20 -Q 20 -r chrM -o coverages/mean_coverage_qfilt.txt;
#echo "mean_coverage";
samtools coverage -b ls_bam_files.txt -r chrM -o coverages/mean_coverage.txt;


## loop for individual bam files
#
## read list of bam files into array
#readarray -t bams < ls_bam_files.txt
#
#for i in "${bams[@]}";
#do
## with filters
#samtools coverage $i -q 20 -Q 20 -r chrM -o coverages/coverage_qfilt_${i}.txt;
## without filters
#samtools coverage $i -r chrM -o coverages/coverage_${i}.txt;
#done
#
#

rm ls_bam_files.txt;




##  ##Get qualimap ##
#if ! [-f "qualimap_v2.2.1/qualimap"]; then
#  wget https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.2.1.zip;
#  unzip qualimap_v2.2.1.zip;
#  rm qualimap_v2.2.1.zip;
#fi
#
#module load Java/11.0.2;

#
# ATACseqQC bioconductor package?
#  ## Run qualimap ##
# 'qualimap multiqc' needs a config file made up of two columns: sample_name path_to_bam
#qualimap_v2.2.1/qualimap multi-bamqc -d qualimap_config.txt -r -c;
#
#
## read bulk ATAC-seq from TF1 cells into array
#readarray -t rts < multiQC/group_SRP149534_SRRs.txt;
#touch qualimap_coverages.txt;
#
#for rt in "${rts[@]}";
#do
#echo ${rt} >> qualimap_coverages.txt;
#grep MT /nobackup/proj/clsclmr/Ludwig_2019/bam/${rt}_sorted_stats/genome_results.txt >> qualimap_coverages.txt;
#
#done
#
#
#
