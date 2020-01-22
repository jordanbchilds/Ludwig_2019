# Ludwig_2019
 Some re-analysis of data from [Ludwig et al. (2019)](https://doi.org/10.1016/j.cell.2019.01.022).

**analyse.sh** is a BASH shell script that downloads sequencing reads from the article, hosted by NCBI and aligns them to the human mitochondrial genome using [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).

**parse.py** is a custom python script called by analyse.sh which gets the connection between GSM IDs (e.g. from metadata [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115218) and SRX IDs, which are needed to download the actual raw data.

**ExamplePath.txt** is a text file containing the IDs of samples which make up one route through eight passages of clonal TF1 cultures presented in Figure 1C.

**findmutations.py** is a custom python script which takes alignment files and converts these into read coverage along the genome as well as mutation load estimates at each nucleobase.  This script also generates .pdf reports plotting coverage and mutation load for each sample.  It also generates plots of mutation load profiles across all eight samples for all 16,569 locations on the genome.
