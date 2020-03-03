# Make sure you have SRA toolkit installed.
# Note that version from ubuntu repos is too out of date:
# https://ncbi.github.io/sra-tools/install_config.html

# Don't forget to add the contents of sra-tools directory to your path
mkdir sra;
mkdir fastq;
mkdir bam;
mkdir pileup;
mkdir frames;
mkdir frames_examine;
mkdir reports;
mkdir nuc;
mkdir mito;

# To download the samples, you might be tempted to use fastq-dump from sra-tools.
# However, this is slow and unable to resume from broken connection.
# Better to run sra-tools prefetch first, then fastq-dump on result

# Use custom python script to download metadata
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115218
# https://doi.org/10.1016/j.cell.2019.01.022
gse='GSE115218';
python3 parse.py $gse;

# Next, update the prefetch download directory as required:
# However, command below doesn't seem to work...  Watch out, huge SRA files stored at ~/ncbi/public/sra
# Need to delete once have .fastq files
echo '/repository/user/main/public/rt = '"\"$(pwd)/sra\"" > $HOME/ncbi/user-settings.mkfg;
prefetch $(<$gse\_sra.txt) 
fastq-dump --outdir fastq $(<$gse\_sra.txt)
rm -rf sra;

readarray -t rts < $gse\_sra.txt;
#readarray -t rts < ExamplePath_sra.txt;

# Download reference human genome
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GCA_000001405.28_GRCh38.p13_genomic.fna.gz
gunzip GCA_000001405.28_GRCh38.p13_genomic.fna.gz;
grep '^>' GCA_000001405.28_GRCh38.p13_genomic.fna > seqnames.txt;

# Split into nuclear sequences and mitochondrial sequences
python3 split_genome.py GCA_000001405.28_GRCh38.p13_genomic.fna;

# Build indices for reference sequences
bowtie2-build mito/mito.fna mito/mito&
bowtie2-build nuc/nuc.fna nuc/nuc;

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
	bowtie2 -p 22 --very-sensitive-local -x nuc/nuc -U fastq/${rt}.fastq -un fastq/${rt}_unmapped.fastq;
	bowtie2 -p 22 --very-sensitive-local -x mito/mito -U fastq/${rt}_unmapped.fastq -S ${rt}_aligned_mito.sam;

    echo "Generating output files...";
    samtools view -Sb  ${rt}_aligned_mito.sam -u| samtools view -h -f 0 -q 1 - >  ${rt}_unsorted.sam;
    samtools view -Sb ${rt}_unsorted.sam -u|samtools sort - bam/${rt}_header;
    #samtools view -h bam/${rt}_header.bam > ${rt}_header.sam
    samtools index bam/${rt}_header.bam bam/${rt}_header.bai;
    rm ${rt}_unsorted.sam;
    rm ${rt}_header.sam;
    rm ${rt}_aligned_mito.sam;
    #samtools mpileup -a bam/${rt}_header.bam > pileup/${rt}.pileup -f mtDNA.fa
  fi

done
