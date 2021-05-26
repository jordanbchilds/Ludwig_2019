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
module load HISAT2/2.1.0-foss-2017b
module load Bowtie2/2.3.4.2-foss-2018b

#Export sratools bin to shell PATH variable

export PATH=$PATH:`pwd`/sratoolkit.2.11.0-ubuntu64/bin/;

gse='GSE115218';

#readarray -t rts < ExamplePath_sra.txt;
readarray -t rts < SRR_Acc_List.txt;


    ## Convert prefetched .sra files to fasta format ##

# loop to convert each prefetched .sra file to fasta format
#for i in "${rts[@]}"
#do
##echo "next"
#echo $i
#fastq-dump --outdir "fastq" "sra/sra/${i}"
#done

#rm -rf sra;



## Split into nuclear sequences and mitochondrial sequences
python3 split_genome.py GCA_000001405.28_GRCh38.p13_genomic.fna;

# Build indices for reference sequences
# The second line takes many hours to complete...
#bowtie2-build mito/mito.fna mito/mito&
#bowtie2-build nuc/nuc.fna nuc/nuc;


hisat2-build -p 8 mito/mito.fna mito/mito;
echo histat2-build mitochondrial indices building stopped
hisat2-build -p 8 nuc/nuc.fna nuc/nuc;
echo histat2-build mitochondrial indices building stopped
 

for rt in "${rts[@]}"
do

 echo ${rt};
 #echo Number of reads: $(cat fastq/${rt}.fastq|wc -l)/4|bc

 if [ -f "bam/${rt}_header.bam" ]; then
    echo "${rt} already aligned";
 else 
    echo "Aligning ${rt} to mitochondrial genome...";
    #bowtie2 -p 22 -D20 -R 10 -N 1 -L 20 -i C,1 -x mito/mito -U fastq/${rt}.fastq -S ${rt}_aligned_mito.sam
    #bowtie2 -p 22 --very-sensitive-local -x mito/mito -U fastq/${rt}.fastq -S ${rt}_aligned_mito.sam;
	bowtie2 -p 8 --very-sensitive -x ./nuc/nuc -U fastq/${rt}.fastq --un ./fastq/${rt}_unmapped.fastq -S fastq/${rt}_tmp.sam;
	bowtie2 -p 8 --very-sensitive -x ./mito/mito -U ./fastq/${rt}_unmapped.fastq -S ${rt}_aligned_mito.sam;

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

