#!/bin/bash
#
#SBATCH --chdir=/nobackup/proj/clsclmr/Ludwig_2019
#SBATCH -p defq
#SBATCH -A clsclmr
#SBATCH -t 24:00:00#SBATCH -c 8
#


module load SAMtools/1.12-GCC-10.2.0;


mkdir alignment_stats/;

bamdir="bam_hg38nodups"

# Extract bam file names (using bam files instead of group SRRs to get bam/prefix, and to avoid excluded samples)
ls -1 -d ${bamdir}/* | grep -v "bai" > ls_bam_files.txt  # ls -1 prints each file on a new line 


  ## SAMtools depth ##

# Calculate depth for all positions (-a), comment line of column names (-H), base (-q) and mapping quality (-Q) greater than 20: (based on default settings for mutserve variant caller), region chrM (-r), remove depth limit (-d 0). Reads with UNMAP, SECONDARY, QCFAIL, or DUP flags are excluded by default.
echo "Calculating depths of filtered reads"
#samtools depth -f ls_bam_files.txt -a -H -q 30 -Q 18 -r chrM -d 0 -o alignment_stats/depths_qfilt_30_${bamdir}.txt
samtools depth -f ls_bam_files.txt -a -H -q 20 -Q 18 -r chrM -d 0 -o alignment_stats/depths_qfilt_${bamdir}.txt
echo "calculating depths of reads (no base or mapping quality filters)"
#samtools depth -f ls_bam_files.txt -a -H -r chrM -d 0 -o alignment_stats/depths_${bamdir}.txt;


  ## SAMtools coverage ##

# calculate mean depth, no reads, mean base quality, mean mapping quality, breadth of coverage and output tab separated file. 
# list of bam files (-b), min base quality (-q), min mapping quality (-Q), mitochondrial chromosome (-r chrM).

echo "calculating mean coverage of filtered reads across all files";
#samtools coverage -b ls_bam_files.txt --excl-flags UNMAP,SECONDARY,QCFAIL,DUP -q 30 -Q 18 -r chrM -o alignment_stats/mean_coverage_qfilt_30_${bamdir}.txt;
samtools coverage -b ls_bam_files.txt --excl-flags UNMAP,SECONDARY,QCFAIL,DUP -q 20 -Q 18 -r chrM -o alignment_stats/mean_coverage_qfilt_${bamdir}.txt;
echo "calculating mean coverage of all reads, no base or mapping quality filters";
#samtools coverage -b ls_bam_files.txt -r chrM -o alignment_stats/mean_coverage_${bamdir}.txt;

# loop for individual bam files

# read list of bam files into array
readarray -t bams < ls_bam_files.txt;

#echo "SRRfile	rname	startpos	endpos	numreads	covbases	coverage	meandepth	meanbaseq	meanmapq" > alignment_stats/all_coverages_${bamdir}.txt;
echo "SRRfile	#rname	startpos	endpos	numreads	covbases	coverage	meandepth	meanbaseq	meanmapq" > alignment_stats/all_coverages_qfilt_${bamdir}.txt;
#echo "SRRfile	#rname	startpos	endpos	numreads	covbases	coverage	meandepth	meanbaseq	meanmapq" > alignment_stats/all_coverages_qfilt_30_${bamdir}.txt;
 

for i in "${bams[@]}";
do
 # remove bam/ prefix from i (bam/SRR*.bam) to get just SRR*.bam
 i_nodir=`echo ${i} | sed "s/${bamdir}\///"`
 i_nodir=${i_nodir/\.bam/}
 echo $i_nodir 
 # # with filters (baseq30)
  echo "${i_nodir}: coverage of filtered reads"
 # samtools coverage $i -q 30 -Q 18 -r chrM --excl-flags UNMAP,SECONDARY,QCFAIL,DUP -o alignment_stats/coverage_qfilt_${i_nodir}.txt;
 # echo "${i_nodir}	`grep chrM alignment_stats/coverage_qfilt_${i_nodir}.txt;`" >> alignment_stats/all_coverages_qfilt_30_${bamdir}.txt
 # rm alignment_stats/coverage_qfilt_${i_nodir}.txt; 
  
 # with jfilters
 samtools coverage $i -q 20 -Q 18 -r chrM --excl-flags UNMAP,SECONDARY,QCFAIL,DUP -o alignment_stats/coverage_qfilt_${i_nodir}.txt;
 echo "${i_nodir}	`grep chrM alignment_stats/coverage_qfilt_${i_nodir}.txt;`" >> alignment_stats/all_coverages_qfilt_${bamdir}.txt
 rm alignment_stats/coverage_qfilt_${i_nodir}.txt; 
 
# ## without filters
# echo "${i_nodir}: coverage of all reads"
# samtools coverage $i -r chrM -o alignment_stats/coverage_${i_nodir}.txt;
# echo "${i}	`grep chrM alignment_stats/coverage_${i_nodir}.txt;`" >> alignment_stats/all_coverages_${bamdir}.txt
# rm alignment_stats/coverage_${i_nodir}.txt;
 
done


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
#readarray -t rts < data/group_SRP149534_SRRs.txt;
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

module purge;
