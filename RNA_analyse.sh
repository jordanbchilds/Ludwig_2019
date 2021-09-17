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
module load HISAT2/2.1.0-foss-2017b;
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

rm -rf sra;

## Download reference genome ##
#wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GCA_000001405.28_GRCh38.p13_genomic.fna.gz ; 
#gunzip GCA_000001405.28_GRCh38.p13_genomic.fna.gz;


echo "hisat2-build reference indices";
hisat2-build -p 8 nuc/hg38.fa nuc/ref;



 ## Align reads ##

for rt in "${rts[@]}"
do

 echo ${rt};
 #echo Number of reads: $(cat fastq/${rt}.fastq|wc -l)/4|bc

 if [ -f "bam/${rt}_sorted.bai" ]; then
    echo "${rt} already aligned";
 else 
    
    echo "Aligning ${rt} to nuclear genome...";
    hisat2 -p 8 -x nuc/ref -1 fastq/${rt}_1.fastq.gz -2 fastq/${rt}_2.fastq.gz --un fastq/${rt}_unmapped.fastq -S sam/${rt}_aligned.sam | samtools view --threads 8 - -h -u | samtools sort --threads 8 - > bam/${rt}_sorted.bam ;
    echo "bam/${rt}_sorted.bam aligned. Indexing..";

    # index sorted bam files
    samtools index -@ 8 bam/${rt}_sorted.bam ;




#    echo "Generating output files...";
#  first filter: -u outputs uncompressed bam into pipe: -h (header), -f 0 (do not output alignments with 0 bits), -q 1 (skip alignments with MAPQ quality <1)    
    #samtools view -@ 8 -Sb sam/${rt}_aligned.sam -u| samtools view -@ 8 -h -f 1 -q 10 > sam/${rt}_unsorted.sam; 
    #samtools view -@ 8 -Sb sam/${rt}_unsorted.sam -u| samtools sort --threads 8 > bam/${rt}_sorted.bam;  # -u pipes bam of rt_unsorted.sam: sorts by reference index
    #samtools view -h bam/${rt}_sorted.bam > ${rt}_header.sam  # why
    samtools index -@ 8 bam/${rt}_sorted.bam bam/${rt}_sorted.bai;  # index sorted bam file
    rm sam/${rt}_unsorted.sam;
    #rm ${rt}_header.sam;
    #rm bam/${rt}_sorted.bam
    rm sam/${rt}_aligned.sam;
 fi
done



#  ##Get qualimap ##
#if ! [-f "qualimap_v2.2.1/qualimap"]; then
#  wget https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.2.1.zip;
#  unzip qualimap_v2.2.1.zip;
#  rm qualimap_v2.2.1.zip;
#fi
#
#module load Java/11.0.2;


# ATACseqQC bioconductor package?
  ## Run qualimap ##
# 'qualimap multiqc' needs a config file made up of two columns: sample_name path_to_bam
#qualimap_v2.2.1/qualimap multi-bamqc -d qualimap_config.txt -r ;

#  ## Get RNA-SeQC ##
#wget http://www.broadinstitute.org/cancer/cga/tools/rnaseqc/RNA-SeQC_v1.1.8.jar
#
#  ## Get gtf annotation file for Hg19 reference genome ##
#wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz



module purge;
