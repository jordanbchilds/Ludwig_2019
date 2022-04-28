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


# set dir and group
bamdir="bam_cnodups"
outdir="bam_cnodups_BAQ_nodups"
j="SRP149534"


mkdir vcf_${outdir}/


  ## check if mutserve is installed
if [ -f "software/bin/mutserve" ]; then 
  echo "mutserve already installed";
else
  echo "Mutserve is not installed, installing mutserve...";
  source install_software.sh;
fi


  ## Variant call
# read bulk ATAC-seq from TF1 cells into array
readarray -t rts < data/group_${j}_SRRs.txt;

if test -f "nuc/parent_consensus.fa"; then
  ref="nuc/parent_chrM_consensus.fa"
else
  ref="software/bin/rCRS.fasta"
fi


echo "Reference fasta file: ${ref}"

#for rt in "${rts[@]}"
#do
#
#  echo ${rt};
#
#  if [ -f "vcf_${outdir}/${rt}_annotated.txt" ]; then
#    echo "${rt} already called";
#  else 
#    
#    # default settings: min heteroplasmy level=0.01, mapping quality=20, base quality=20, alignment quality=30
#    echo "Calling variants for ${rt}...";
#    mutserve call ${bamdir}/${rt}.bam --threads 8 --baseQ 30 --mapQ 18 --level 0.001 --reference ${ref} --output vcf_${outdir}/${rt}.vcf.gz ;
#
#    gunzip vcf_${outdir}/${rt}.vcf.gz
#
#    echo "left aligning with bcftools";
#    bcftools norm vcf_${outdir}/${rt}.vcf -f ${ref} --multiallelics -any -o vcf_${outdir}/${rt}_normalised.vcf -Ov
#    
#    echo "Annotating mutserve .txt output";
#    mutserve annotate --input vcf_${outdir}/${rt}.txt --annotation software/bin/rCRS_annotation_2020-08-20.txt --output vcf_${outdir}/${rt}_annotated.txt
#  fi
#done
#


    ###############  Pileups  ##################

# get read pileups with bcftools: includes depth for allele for each strand 
# mpileup includes all reads in FORMAT/ADF for example, (mapq and baseq filters only apply to genotype calling), so the filtered bam file must be piped to bcftools pileup

mkdir mpileups_${outdir}/

# index reference
samtools faidx ${ref};


# Create pileups:
for rt in "${rts[@]}"
do

  echo ${rt};

  if [ -f "mpileups_${outdir}/${rt}_calls.vcf" ]; then
    echo "${rt} already called";
  else 
    
    echo "Creating mpileup for ${rt}...";
    samtools view ${bamdir}/${rt}.bam chrM -h -u -F DUP | bcftools mpileup - --max-depth 999999 --fasta-ref ${ref} -q 18 -Q 30 --annotate FORMAT/SP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR --threads 8 -Ou | bcftools norm -f ${ref} -m -any -Ov --output mpileups_${outdir}/${rt}_mpileup.vcf 
#- | bcftools query -f '%CHROM\t%POS\t%REF\t%DP\t%RO\t%AO\t%ALT\t%SAF\t%SAR\t%SRF\t%SRR\n' --output mpileups_${bamdir}/${rt}_mpileup.vcf 
#FORMAT/AD,FORMAT/ADF,FORMAT/ADR,
    # separate final column with AD, ADF, ADR and SP info into separate columns: replace ";" separator with tabs (\t) 
    #sed "s/;/\t/g" -i mpileups_${bamdir}/${rt}_mpileup.vcf
    grep -v -P '\t<\*>' mpileups_${outdir}/${rt}_mpileup.vcf > mpileups_${outdir}/${rt}_mpileup_nodels.vcf
    # replace "-" in #headers row with tab separated names of columns
    #sed "s/INFO    FORMAT  -/INFO    FORMAT  AD	ADF	ADR	SP/g" -i mpileups_${bamdir}/${rt}_mpileup.vcf
    
    #echo "Calling point mutations (bcftools call) for ${rt}...";
    #bgzip -i -c mpileups_${bamdir}/${rt}_mpileup.vcf --threads 8 | 
    #bcftools call mpileups_${bamdir}/${rt}_mpileup.vcf --multiallelic-caller --keep-alts --skip-variants indels --regions chrM -Ov --ploidy 1 --output bcf_calls_${bamdir}/${rt}_calls.vcf;
  fi
done

module purge;
