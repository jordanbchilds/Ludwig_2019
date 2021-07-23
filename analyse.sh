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

echo "modules loaded";

# Export sratools bin to shell PATH variable
export PATH=$PATH:`pwd`/sratoolkit.2.11.0-ubuntu64/bin/;

gse='GSE115218';

# read bulk ATAC-seq from TF1 cells into array
readarray -t rts < multiQC/group_SRP149534_SRRs.txt;



    ## Convert prefetched .sra files to fasta format ##

# loop to find if .sra file has been dumped (converted to fastq.gz) and if not add file to list

for i in "${rts[@]}";
do
if test -f "fastq/${i}_1.fastq.gz";
then
  echo "${i}_1.fastq.gz file exists";
else
  echo ${i} >> dump_list.txt
  echo "${i}.sra added to dump_list.txt";
fi
done

cat dump_list.txt | parallel --jobs 8 "fastq-dump --split-files --gzip --outdir fastq/ sra/sra/{}.sra";
rm dump_list.txt;

#rm -rf sra;

 ## Download reference genome ##
#cd nuc/;
#wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
#gunzip hg38.fa.gz
#cd ..;

#echo "bowtie2-build reference indices";
#bowtie2-build --threads 8 nuc/hg38.fa nuc/btref;

locale;
export LANG=en_GB.utf8
export LC_ALL="en_GB.utf8" 
locale;

  ## Align reads ##

for rt in "${rts[@]}"
do

 echo ${rt};
 #echo Number of reads: $(cat fastq/${rt}.fastq|wc -l)/4|bc

 if [ -f "bam/${rt}_sorted.bam" ]; then
    echo "${rt} already aligned";
 else 
   
    echo "Aligning ${rt} to whole genome...";
    # bowtie2 parameters: -p 8 cores, forward and reverse read, ref, local alignment (soft-clipping allowed), very sensitive (-L 20: 20 bp substrings in multiseed, -i s,1,0.50: shorter intervals between seed substrings, -D 20 -R 3: see manual), -t: time to align in stout,  out? -X 2000???.
    # samtools view parameters:  first filter: - (input from stdin), -h (header), eg. -F 0 (do not output alignments with FLAG integer), eg. -q 10 (skip alignments with MAPQ quality <10), -u outputs uncompressed bam into pipe.
    bowtie2 -p 8 -1 fastq/${rt}_1.fastq.gz -2 fastq/${rt}_2.fastq.gz -x nuc/btref --local --sensitive -t --un-gz fastq/${rt}_unmapped.fastq | samtools view --threads 8 - -h -u | samtools sort --threads 8 - > bam/${rt}_sorted.bam ;
    echo "bam/${rt}_sorted.bam finished. Indexing..";
    
    # index sorted bam files
    samtools index -@ 8 bam/${rt}_sorted.bam ;

    # create bam containing only mitochondrial aligned reads
#    samtools view --threads 8 -b -h bam/${rt}_sorted.bam MT > bam/${rt}_sorted_chrM.bam;

 fi
done


# calculate coverage
#samtools coverage -b multiQC/group_SRP149534_SRRs.txt -r MT



locale;
export LANG=C.UTF-8 ;                                                                                                   export LC_ALL= ;
locale;



module purge;
