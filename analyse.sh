#!/bin/bash
#
#SBATCH --workdir=/nobackup/proj/clsclmr/Ludwig_2019
#SBATCH -p defq
#SBATCH -A clsclmr
#SBATCH -t 48:00:00
#SBATCH -c 8
#

    ## load modules
module load Python/3.8.6-GCCcore-10.2.0;
module load SAMtools/1.12-GCC-10.2.0; 
module load HISAT2/2.1.0-foss-2017b;
module load Bowtie2/2.3.4.2-foss-2018b;

#Export sratools bin to shell PATH variable
export PATH=$PATH:`pwd`/sratoolkit.2.11.0-ubuntu64/bin/;

gse='GSE115218';

#readarray -t rts < ExamplePath_sra.txt;
readarray -t rts < SRR_Acc_List.txt;


    ## Convert prefetched .sra files to fasta format ##

## loop to convert each prefetched .sra file to fasta format
#for i in "${rts[@]}"
#do
##echo "next"
#echo $i
#fastq-dump --outdir "fastq" "sra/sra/${i}"
#done

#rm -rf sra;



## Split into nuclear sequences and mitochondrial sequences
#python3 split_genome.py GCA_000001405.28_GRCh38.p13_genomic.fna;

## Build indices for reference sequences
#bowtie2-build --threads 8 mito/mito.fna mito/mito&
#bowtie2-build --threads 8 nuc/nuc.fna nuc/nuc;


#hisat2-build -p 8 mito/mito.fna mito/mito;
#echo histat2-build mitochondrial indices building stopped
#hisat2-build -p 8 nuc/nuc.fna nuc/nuc;
#echo histat2-build mitochondrial indices building stopped
 
# remove .fna files so bowtie2 recognises indexes properly when aligning
#rm mito/mito.fna;
#rm nuc/nuc.fna;


for rt in "${rts[@]}"
do

 echo ${rt};
 #echo Number of reads: $(cat fastq/${rt}.fastq|wc -l)/4|bc

 if [ -f "bam/${rt}_header.bam" ]; then
    echo "${rt} already aligned";
 else 
    
    echo "Aligning ${rt} to mitochondrial genome...";
    hisat2 -p 8 -x /nobackup/proj/clsclmr/Ludwig_2019/nuc/nuc -U /nobackup/proj/clsclmr/Ludwig_2019/fastq/${rt}.fastq --un /nobackup/proj/clsclmr/Ludwig_2019/fastq/${rt}_unmapped.fastq -S /nobackup/proj/clsclmr/Ludwig_2019/fastq/${rt}_tmp.sam;
    hisat2 -p 8 -x /nobackup/proj/clsclmr/Ludwig_2019/mito/mito -U /nobackup/proj/clsclmr/Ludwig_2019/fastq/${rt}_unmapped.fastq -S ${rt}_aligned_mito.sam;

    echo "Generating output files...";
    samtools view -Sb ${rt}_aligned_mito.sam -u| samtools view -h -f 0 -q 1 - >  ${rt}_unsorted.sam;
    samtools view -Sb ${rt}_unsorted.sam -u|samtools sort - bam/${rt}_header;
    samtools view -h bam/${rt}_header.bam > ${rt}_header.sam
    samtools index bam/${rt}_header.bam bam/${rt}_header.bam;
    rm ${rt}_unsorted.sam;
    rm ${rt}_header.sam;
    rm ${rt}_aligned_mito.sam;
 fi

done
module purge;
