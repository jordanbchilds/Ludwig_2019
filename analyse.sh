



# Read sra accession names into shell array rts
#readarray -t rts < ExamplePath_sra.txt;
readarray -t rts < $gse\_sra.txt;


    ## Convert prefetched .sra files to fasta format ##

cd sra/sra;

# loop to convert each prefetched .sra file to fasta format
for i in "${rts[@]}"
do
fastq-dump --outdir ../../fastq ${i}.sra
done

rm -rf sra;



## Split into nuclear sequences and mitochondrial sequences
python3 split_genome.py GCA_000001405.28_GRCh38.p13_genomic.fna;

# Build indices for reference sequences
# The second line takes many hours to complete...
#bowtie2-build mito/mito.fna mito/mito&
#bowtie2-build nuc/nuc.fna nuc/nuc;

#
hisat2-build mito/mito.fna mito/mito;
hisat2-build nuc/nuc.fna nuc/nuc;


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
	#bowtie2 -p 22 --very-sensitive -x nuc/nuc -U fastq/${rt}.fastq --un fastq/${rt}_unmapped.fastq -S fastq/${rt}_tmp.sam;
	bowtie2 -p 22 --very-sensitive -x mito/mito -U fastq/${rt}_unmapped.fastq -S ${rt}_aligned_mito.sam;

    echo "Generating output files...";
    samtools view -Sb ${rt}_aligned_mito.sam -u| samtools view -h -f 0 -q 1 - >  ${rt}_unsorted.sam;
    samtools view -Sb ${rt}_unsorted.sam -u|samtools sort - bam/${rt}_header;
    #samtools view -h bam/${rt}_header.bam > ${rt}_header.sam
    samtools index bam/${rt}_header.bam bam/${rt}_header.bam;
    rm ${rt}_unsorted.sam;
    rm ${rt}_header.sam;
