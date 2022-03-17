#!/bin/bash
#SBATCH --chdir=/nobackup/proj/clsclmr/Ludwig_2019
#SBATCH -p defq
#SBATCH -A clsclmr
#SBATCH -t 04:00:00
#SBATCH -c 8
#

  ## load modules

module load BCFtools/1.10.2-foss-2019b;
module load SAMtools/1.12-GCC-10.2.0; 


  ## Create consensus sequence to use instead of rCRS reference sequence in hg38.fa ##
# aligning and variant calling from rCRS reference mitochondrial genome calls lots of homozygous mutations that arose and became fixed in these clones after diverging from their (hypothetical) most recent common ancestor with rCRS genome. Calculating and aligning to a consensus sequence of the parent clone allows only heterozygous mutations undergoing clonal expansion to be called.

# bcftools consensus doesn't account for allele frequency, just replaces reference allele with all alternative alleles in vcf when copying the reference.fa. Consequently homozygous variant calls (AF > 0.95, likely already fixed in parent population) are extracted first:

# extract fixed homozygous variant calls (AF>0.95) from bulk parent .vcfs (SRR7245880.vcf.gz, SRR7245881.vcf.gz) into combined vcf.
grep '^\#' vcf/SRR7245880.vcf > nuc/fixed_calls.vcf  # extract header lines in vcf: all lines begin (^) with '#'.
grep $'\t1:' vcf/SRR7245880.vcf >> nuc/fixed_calls.vcf  # extract homozygous calls: genotype column in GT:AF:DP = \t 1:(AF):(DP)
grep $'\t1:' vcf/SRR7245881.vcf >> nuc/fixed_calls.vcf  # extract in bulk sequencing replicate.
# Remove duplicate variant calls
sort -u -t$'\t' -k2n nuc/fixed_calls.vcf > nuc/parent_fixed.vcf
rm nuc/fixed_calls.vcf
# get consensus sequence of bulk parent clones from vcf SRR7245880.vcf.gz  
cat nuc/hg38.fa | bcftools consensus nuc/parent_fixed.vcf > nuc/parent_consensus.fa


