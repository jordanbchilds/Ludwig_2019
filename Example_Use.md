### Example Usage for bulk ATAC-seq of TF1 cells

**analyse.sh**
====================================

Maps the sequencing reads from the .fasta files to the reference genome's sequence, to create an aliginment in the .bam format. Once aligned, mutations can be identified, eg.
                   TTGGGGACTCTGG
                   TTGGGGACTC  <- successfully aligned read
                  TTAGGGGAC  <- successfully aligned read with a potential T->C mutation 
             TCGCGTTTGGG
         GGATTCGC
ref seq:   ATTCGCGTTTGGGGACTCT

### Stages in the script:
1. Parse arguments
2. Index reference sequence
  - Indexing makes it easier for the alignment program to search reference sequence more easily, like bookmarking and skipping to a page instead of searching a whole book.
3. Read alignment
  - Involves sorting the reads, identifying read pairs and exact duplicates
4. Extract alignment statistics
  - eg. total no. reads, no. reads aligned, no. duplicates etc.

### Output
.bam files, not human readable unless converted and need to be summarised to get information from them.
alignment_stats/alignment_summary.txt, gives overall alignment % per clone.

### Information needed to run: (`bash analyse.sh -h` for help/options)
- name of the reference sequence fasta file: either "nuc/hg38" (a version of the whole human genome) or "nuc/parent_consensus" (consensus sequence of parent clones)
- group name / previously chosen keyword used to choose and extract sequencing runs from data/SraRunTable.txt, eg. B11 (five clones in the B11 lineage), SRP149534 (sequencing runs of ATAC-seq of TF1 clones)
- choose a name for the outdir, eg. "bam_hg38_B11". Don't include the "/".

### Execute
eg. `bash analyse.sh --reference hg38 --group-name B11 --bam-directory bam_hg38_B11`


**post_alignment.sh**
===============================================
Outputs files with detailed stats of alignment. For each clone: mean coverage, base quality, per genomic position coverage and base quality. Repeated for different thresholds.

This script summarises the alignment but is less important now, because the variant caller gives the same information but only for the exact reads used in the variant calling. 

### Stages in the script
1. Overall, Per position, and after filtering out reads with mapping and base quality filters.

### Output
alignment\_stats/\*

### Information needed to run: 
- group name / previously chosen keyword used to choose and extract sequencing runs from data/SraRunTable.txt, eg. B11 (five clones in the B11 lineage), SRP149534 (sequencing runs of ATAC-seq of TF1 clones)

### Execute
- Change the bamdir to the group name in the script (Doesn't use command line arguments yet): eg. bamdir="B11" 
- Execute script: `bash post_alignment.sh`


**variant\_call.sh**
===============================================

Takes the alignment in the .bam files and calls all mutations and their allele frequencies. Pileup includes all raw mutations including multiallelics (after filtering for base and mapping qualities). It doesn't call mutations/"decide" which ones are valid by applying any other filters, thresholds or models: just all alleles and the no. times they occur at any site. Also outputs information in per site coverage/depth, a calculation of strand bias, no. reads supporting an allele on each strand etc.

### Stages in the script
1. In bam/ directory name, group name, Makes out directory 
2. Reads in names of clones
3. Checks for reference fasta
4. Creates pileup of all alleles at each position

### Execute
Need to know:
- name of input bam directory
- choose name of output directory (name after bam dir is informative)
- group name 
Then execute:
- change input and output directories in the script, and the group name.
- `bash variant_call.sh`
