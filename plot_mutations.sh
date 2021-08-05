#!/bin/bash
#
#SBATCH --chdir=/nobackup/proj/clsclmr/Ludwig_2019
#SBATCH -p defq
#SBATCH -A clsclmr
#SBATCH -t 00:30:00
#SBATCH -c 8
#
module load R/3.6.2-fosscuda-2019b

mkdir plots

loc=$PWD
echo $loc
Rscript mutation_plot_script.R $loc

