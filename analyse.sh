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

gse='GSE115218';

# read bulk ATAC-seq from TF1 cells into array
readarray -t rts < group_SRP149534_SRRs.txt;


  ## build indices ##
if test -f "nuc/btref.1.bt2"; then
  echo "bowtie2-build reference indices";
else
  bowtie2-build --threads 8 nuc/hg38.fa nuc/btref;
fi

#locale;
#export LANG=en_GB.utf8
#export LC_ALL="en_GB.utf8" 
locale;


  ## Align reads ##

for rt in "${rts[@]}"
do

 echo ${rt};
 #echo Number of reads: $(cat fastq/${rt}.fastq|wc -l)/4|bc

 if [ -f "bam/${rt}_sorted.bam.bai" ]; then
    echo "${rt} already aligned";
 else 
   
    echo "Aligning ${rt} to whole genome...";

# bowtie2 parameters: -p 8 cores, -1 forward and -2 reverse read, -x ref, --local alignment (soft-clipping allowed), very sensitive (-L 20: 20 bp substrings in multiseed, -i s,1,0.50: shorter intervals between seed substrings, -D 20 -R 3: see manual), -t: time to align in stout,  out? -X 2000???.

# samtools view parameters: - (input from stdin), -h (header), eg. -F 0 (do not output alignments with FLAG integer), eg. -q 10 (skip alignments with MAPQ quality <10), --un-gz unmapped reads to zipped file, -u outputs uncompressed bam into pipe.

# samtools fixmate adds mate score (ms) tags (used by markdup to select best reads): -m ms tags.

# samtools sort: sorts by leftmost coordinates for indexing

# samtools markdup: -s print basic stats, -r REMOVE duplicates, 

    bowtie2 -p 8 -1 fastq/${rt}_1.fastq.gz -2 fastq/${rt}_2.fastq.gz -x nuc/btref --local --sensitive -t --un-gz fastq/${rt}_unmapped.fastq | \
samtools view --threads 8 - -h -u | \
samtools fixmate --threads 8 - -m -u | \
samtools sort --threads 8 - -u | \ 
samtools markdup --threads 8 - -s > bam/${rt}_sorted.bam ;
# ^add removal of optical duplicates with -d _ : what distance for Nextseq?

    # index sorted bam files: -b .bai
    echo "bam/${rt}_sorted.bam aligned. Indexing.."; 
    samtools index --threads 8 -b bam/${rt}_sorted.bam ;

    # create bam containing only mitochondrial aligned reads
#    samtools view --threads 8 -b -h bam/${rt}_sorted.bam MT > bam/${rt}_sorted_chrM.bam;
 
    # get bam files headers: (-H header, -? human readable)
    samtools view -H -? bam/${rt}_sorted.bam >> ../bam_summaries.txt 

 fi
done


locale;
export LANG=C.UTF-8 ; 
export LC_ALL= ;
locale;


# Copy text from slurm outfile to alignment stats
cp slurm-${SLURM_JOB_ID}.out alignment_stdout.txt
echo 'SRR Overall_alignment_rate' > alignment_summary.txt

for rt in "${rts[@]}"
do

 echo ${rt};
 overall=`grep -A 20 "$rt" alignment_stdout.txt | grep "overall" | cut -d "%" -f 1`;
 echo "$rt $overall" >> alignment_summary.txt;

done


module purge;
