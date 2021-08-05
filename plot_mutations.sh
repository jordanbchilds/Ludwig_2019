#!/bin/bash
#
#SBATCH --chdir=/nobackup/proj/clsclmr/Ludwig_2019
#SBATCH -p defq
#SBATCH -A clsclmr
#SBATCH -t 00:30:00
#SBATCH -c 8
#
module load R/3.6.0-foss-2019a 

mkdir plots
mkdir .R_local_lib/
 

loc=$PWD
echo $loc
Rscript mutation_plot_script.R $loc

