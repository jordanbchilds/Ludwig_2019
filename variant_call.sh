#!/bin/bash
#SBATCH --chdir=/nobackup/proj/clsclmr/Ludwig_2019
#SBATCH -p defq
#SBATCH -A clsclmr
#SBATCH -t 24:00:00
#SBATCH -c 8
#

    ## load modules
module load Java/11.0.2; 
module load BCFtools/1.10.2-foss-2019b;
module load SAMtools/1.12-GCC-10.2.0; 

export PATH=`pwd`/software/bin/:$PATH

mkdir vcf_consensus-nodups/

  ## check if mutserve is installed
if [ -f "software/bin/mutserve" ]; then 
  echo "mutserve already installed";
else
  echo "installing mutserve...";
  mkdir mutserve/;
  cd mutserve/;
  source ../scripts/mutserve_installer.sh;
  cd ../;
fi

  ## Variant call
# read bulk ATAC-seq from TF1 cells into array
readarray -t rts < data/group_SRP149534_SRRs.txt;

if test -f "nuc/parent_consensus.fa"; then
  ref="nuc/parent_chrM_consensus.fa"
else
  ref="mutserve/rCRS.fasta"
fi


echo "Reference fasta file: ${ref}"

for rt in "${rts[@]}"
do

  echo ${rt};

  if [ -f "vcf_consensus-nodups/${rt}_annotated.txt" ]; then
    echo "${rt} already called";
  else 
    
    # default settings: min heteroplasmy level=0.01, mapping quality=20, base quality=20, alignment quality=30
    echo "Calling variants for ${rt}...";
    ./mutserve/mutserve call bamc/${rt}.bam --threads 8 --baseQ 20 --mapQ 18 --level 0.001 --reference ${ref} --output vcf_consensus-nodups/${rt}.vcf.gz ;

    gunzip vcf_consensus-nodups/${rt}.vcf.gz

    echo "left aligning with bcftools";
    bcftools norm vcf_consensus-nodups/${rt}.vcf -f ${ref} --multiallelics +snps -o vcf_consensus-nodups/${rt}_normalised.vcf -Ov
    
    echo "Annotating mutserve .txt output";
    ./mutserve/mutserve annotate --input vcf_consensus-nodups/${rt}.txt --annotation mutserve/rCRS_annotation_2020-08-20.txt --output vcf_consensus-nodups/${rt}_annotated.txt
  fi
done



# get read depth for allele on each strand using bcftools
# mpileup includes all reads in FORMAT/ADF for example, (mapq and baseq filters only apply to genotype calling), so the filtered bam file must be piped to bcftools pileup

mkdir bcftools_out/;

# index reference
samtools faidx ${ref};

for rt in "${rts[@]}"
do

  echo ${rt};

  if [ -f "bcftools_out/${rt}_calls.vcf" ]; then
    echo "${rt} already called";
  else 
    
    echo "Creating mpileup for ${rt}...";
    samtools view bamc/${rt}.bam chrM -h -u | bcftools mpileup - --no-BAQ --max-depth 999999 --fasta-ref ${ref} -q 18 -Q 20 --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR --threads 8 -Ov --output bcftools_out/${rt}_mpileup.vcf;
    echo "Calling point mutations (bcftools call) for ${rt}...";
    #bgzip -i -c bcftools_out/${rt}_mpileup.vcf --threads 8 | bcftools call - --multiallelic-caller --keep-alts --skip-variants indels --regions chrM -Ov --ploidy 1 --output bcftools_out/${rt}_calls.vcf;
  fi
done

module purge;
