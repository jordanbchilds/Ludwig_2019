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
module load SAMtools;  # which version

    ## Make sure you have SRA toolkit installed.
# if not see configure_sratoolkilt.sh
# (implement a check for SRA toolkit)
# Export to shell PATH variable
export PATH=$PATH:`pwd`/sratoolkit.2.11.0-ubuntu64/bin/;

# make directories
mkdir fastq;
mkdir bam;
mkdir pileup;
mkdir frames;
mkdir frames_examine;
mkdir reports;
mkdir nuc;
mkdir mito;
mkdir sra;


# To download the samples, you might be tempted to use fastq-dump from sra-tools.
# However, this is slow and unable to resume from broken connection.
# Better to run sra-tools prefetch first, then fastq-dump on result


     ## Use custom python script to download metadata ##

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115218
# https://doi.org/10.1016/j.cell.2019.01.022

gse='GSE115218';
# parse.py gets list of SRA sequence names from GSE series of "Human lineage tracing enabled by mitochondrial mutations and single cell genomics"

# install GEOparse python3 module. pip3 install _ won't install if module is already installed.
pip3 install GEOparse;

# run parse.py
python3 parse.py $gse;


    ## Prefetch .sra files ##
# Need to delete once have .fastq files.
#vdb-dump --info

prefetch --option-file $gse\_sra.txt;
# fasterq-dump????


    ## Download reference human genome ##

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GCA_000001405.28_GRCh38.p13_genomic.fna.gz ; 
gunzip GCA_000001405.28_GRCh38.p13_genomic.fna.gz;
rm GCA_000001405.28_GRCh38.p13_genomic.fna.gz;
grep '^>' GCA_000001405.28_GRCh38.p13_genomic.fna > seqnames.txt;



