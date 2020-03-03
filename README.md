# Ludwig_2019
 Some re-analysis of data from [Ludwig et al. (2019)](https://doi.org/10.1016/j.cell.2019.01.022).

**analyse.sh** is a BASH shell script that downloads the reference human genome and sequencing reads from the article, hosted by NCBI.  Splits reference sequence into nuclear and mitochondrial reference sequences, builds indices and aligns reads using [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).  First aligns reads to nuclear reference, saving reads which do not aligned.  These are then aligned to the mitochondrial reference.  In this way, we heed the note of caution raised by [Santibanez-Koref et al. (2019)](https://doi.org/10.1016/j.mito.2018.08.003).  

**parse.py** is a custom python script called by analyse.sh which gets the connection between GSM IDs (e.g. from metadata [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115218) and SRX IDs, which are needed to download the actual raw data.

**split_genome.py** is a custom python script called by analyse.sh which separates individual sequences (e.g. chromosomes & mtDNA sequence) into a nuclear reference (nuc.fna) and a mitochondrial reference sequence (mito.fna).

**ExamplePath.txt** is a text file containing the IDs of samples which make up one route through eight passages of clonal TF1 cultures, highhlighted in yellow asterisks in this adapted version of Figure 1C:

[<img src="reports/LudwigFigs.png">](https://doi.org/10.1016/j.cell.2019.01.022)

**findmutations.py** is a custom python script which takes alignment files and converts these into read coverage along the genome as well as mutation load estimates at each nucleobase.  This script also generates .pdf reports plotting coverage and mutation load for each sample.

<img src="reports/bulk_mutation_load.png">

It also generates plots of mutation load profiles across all eight samples for all 16,569 locations on the genome.

<img src="reports/01492.png">
