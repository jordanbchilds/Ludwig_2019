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
module load Bowtie2/2.3.4.2-foss-2018b;
module load parallel/20200522-GCCcore-10.2.0

mkdir bam;

gse='GSE115208';

# read bulk ATAC-seq from TF1 cells into array
readarray -t rts < group_SRP149534_SRRs.txt;


  ## build indices ##

if test -f "nuc/parent_consensus.1.bt2"; then
  ref="nuc/parent_consensus"
  echo "Parent clone consensus sequence available and indexed: nuc/parent_consensus.1.bt2"
elif test -f "nuc/parent_consensus.fa"; then
  echo "parent clone consensus fasta available, indexing..."
  bowtie2-build --threads 8 nuc/parent_consensus.fa nuc/parent_consensus; 
  echo "bowtie2 referenece indices built"
  ref="nuc/parent_consensus"
 
elif test -f "nuc/hg38.1.bt2"; then
  echo "bowtie2-build reference indices already built";
  ref="nuc/hg38"
else
  echo "building reference indices..."
  bowtie2-build --threads 8 nuc/hg38.fa nuc/hg38;
  ref="nuc/hg38" 
fi


#locale;
#export LANG=en_GB.utf8
#export LC_ALL="en_GB.utf8" 
#locale;


  ## Align reads ##

for rt in "${rts[@]}"
do

 echo ${rt};
 #echo Number of reads: $(cat fastq/${rt}.fastq|wc -l)/4|bc

 if [ -f "bam/${rt}.bam.bai" ]; then
  echo "${rt} already aligned";
 else 
   
  echo "Aligning ${rt} to ${ref} reference genome...";

# bowtie2 parameters: 
# -p 8 cores, -1 forward and -2 reverse read, -x ref, --local alignment (soft-clipping allowed), very sensitive (-L 20: 20 bp substrings in multiseed, -i s,1,0.50: shorter intervals between seed substrings, -D 20 -R 3: see manual), -t: time to align in stout,  out? -X 2000???.

# samtools view parameters: 
# - - (input from stdin, AGAIN to output to stout), -h (header), eg. -F 0 (do not output alignments with FLAG integer), eg. -q 10 (skip alignments with MAPQ quality <10), --un-gz unmapped reads to zipped file, -u outputs uncompressed bam into pipe.

# samtools sort: normally sorts by leftmost coordinates (for indexing), -n sorts by QNAME (for fixmate - so mate pairs can be labelled) 

# samtools fixmate adds mate score (ms) tags (used by markdup to select best reads): -m ms tags.

# samtools markdup: -s print basic stats, -r will REMOVE duplicates, 

  bowtie2 -p 8 -1 fastq/${rt}_1.fastq.gz -2 fastq/${rt}_2.fastq.gz -x ${ref} --local --sensitive -t --un-gz fastq/${rt}_unmapped.fastq | samtools view --threads 8 - -h -u | samtools sort --threads 8 -n - -u | samtools fixmate --threads 8 -m -u - - | samtools sort --threads 8 - -u | samtools markdup --threads 8 -s -u - - > bam/${rt}.bam

# Duplicate commands using intermediate files instead of piping: (after sorting by QNAME)
#echo "collate and fixmate"
##samtools collate --threads 8 bam/${rt}_sorted.bam -O | 
#samtools fixmate --threads 8 -m bam/${rt}_sorted1.bam bam/${rt}_fixmated.bam
#rm bam/${rt}_sorted1.bam  
#echo "sort by coordinates"
#samtools sort --threads 8 bam/${rt}_fixmated.bam -o bam/${rt}_sorted2.bam
#rm bam/${rt}_fixmated.bam 
#echo "markdup"
#samtools markdup --threads 8 -s bam/${rt}_sorted2.bam bam/${rt}.bam

# TODO add removal of optical duplicates with -d _ : what distance for Nextseq?

  # index sorted bam files:
  echo "Index"; 
  samtools index -@ 8 bam/${rt}.bam ;
    
# create bam containing only mitochondrial aligned reads
#    samtools view --threads 8 -b -h bam/${rt}_sorted.bam chrM > bam/${rt}_sorted_chrM.bam;
 
    # get bam files headers: (-H header)
    #samtools view -H bam/${rt}_sorted.bam >> ../bam_summaries.txt 

 fi
done


#locale;
export LANG=C.UTF-8 ; 
export LC_ALL= ;
locale;


# Copy stats from slurm outfile (stout) to alignment stats
cp slurm-${SLURM_JOB_ID}.out alignment_stats/alignment_stdout.txt
echo 'SRR Overall_alignment_rate Total_reads_bowtie2 Total_reads_markdup Total_duplicates Estimated_unique_lib_size' >> alignment_stats/alignment_and_duplicate_summary.txt

for rt in "${rts[@]}"
do
 overall_alignment=`grep -A 20 "$rt" alignment_stats/alignment_stdout.txt | grep "overall" | cut -d "%" -f 1`;
 bowtie2_total_reads=`grep -A 20 "$rt" alignment_stats/alignment_stdout.txt | grep "reads; of these:" | cut -d " " -f 1`;
 markdup_total_reads=`grep -A 41 "$rt" alignment_stats/alignment_stdout.txt | grep "READ:" | cut -d " " -f 2`;
 total_duplicates=`grep -A 41 "$rt" alignment_stats/alignment_stdout.txt | grep "DUPLICATE TOTAL:" | cut -d " " -f 3`;
 unique_lib_size=`grep -A 41 "$rt" alignment_stats/alignment_stdout.txt | grep "ESTIMATED_LIBRARY_SIZE" | cut -d " " -f 2`; 
 echo "$rt $overall_alignment $bowtie2_total_reads $markdup_total_reads $total_duplicates $unique_lib_size" >> alignment_stats/alignment_and_duplicate_summary.txt 
done


module purge;
