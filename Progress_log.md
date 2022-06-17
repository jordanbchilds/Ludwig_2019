## Progress Log  
### Week 1  
- Set up Github page
- Added description to Github 
- Re-learnt basic Git commands 
- Learnt Git branching and merging theory and commands
- Re-learnt Git markdown
- Re-learnt how to connect via and transfer files via ssh and sftp connections
- Connected to c0078068@unix.ncl.ac.uk server via ssh, then the rocket.hpc via ssh
- Forked CnrLwlss/Ludwig\_2019 into my Github and cloned into group folder in rocket cluster
- Met with Conor and Dasha
- Read through and made notes on BASH scripting basics: http://linuxsig.org/files/bash_scripting.html
- Wrote and added this log

Aims
- Overall - get Conor's script working 
  - Automatically download SRA tools
  - (Only use python tools already installed on hpc cluster otherwise - biopython?)

- Downloaded and installed sra-tools, exported to $PATH.

Questions

1. How to stop remote cluster from not responding after no use - is using `top` cpu usage check inconsiderate?
2. Aligning tool hisat2 download 

---


### Week 2
- Fix stupid laptop wifi
- Follow bash tutorial
- Use vim during all text editing
- Fix git push port problem
- Create test directory and analyse_no_download script with necessary downloads in directory
- ready download of hisat2 - unnecesary

Aims
- Run script without automatic downloading - specifically alignment
- quality assess mitochondrial assembly (FastQC?)
- Run with nuclear genome assembly indexes built and reads aligned
- fully understand SRA, GSE files and info.
- implement checks for sra tools and fastq files 
- read through hpc rocket cluster tutorial 
- add instructions for cluster eg. how long  
- run line by line and see if it runs as expected 
- write up details of each line, and how to improve
- add load modules to beginning of script
 - get running then estimate time, no cores etc.

Questions
1. importance of `make` compiling
2. ~~python3 not on system???~~
3. ~~how to view SDERR~~
4. how to view progress of alignment
5. why does .gitignore not exclude sra-toolkit
6. add unload modules to bottom of script?
7. how to calculate no. cores to use in a job and no. tasks.
---


### Week 3 + 4

- add instructions to load modules
- read hpc rocket cluster tutorial
- add preliminary SBATCH commands
- add instructions to install python GEOparse module
- Configure sra-tools properly
- Prefetch working! (only with interactive installation of sra-tools)
- fastq-dumkp working - see test_script.sh

Aims
- implement check (and installation) of python module GEOparse <- without use of venv?
- output just .bam format
- ~~configure sra-toolkit, see which files changed most recently and how~~
- ~~make sure fast(er?)q-dump takes from prefetched files~~
- Prefetch all sra files and convert to fasta
- Split scripts into appropriate scripts eg. download everything in one script
- align to nuclear and mitochondrial genome
- change permissions of project folder for Dasha
- adjust for circularity of mitochondrial genome

Questions 
1. Estimate size - 2733 SRA accessions. Only ~57 .sra /.sra.cache prefetched in a bit less than 1 hr
2. semi-colon use after loops


---

### Week 5 
- Prefetch all files
- add module versions
- fastq-dump all files
- hisat2 alignment working
- samtools output of alignment working

Aims
- update README.md 
- ~~fastq-dump all files~~
- align reads to indices
- link SRA (accessions) to SRR (runs) to check if all are downloaded
- adjust for cirularity of mt genome
- check quality of alignments and reads

Questions
- quality control - check using fastqc but impractical. Trim using timmomatic, cutadapt, **fastp**?
- suitable PHRED quality threshold for trimming - 20?

---

### Week 6
- fastqc on all 
- multiqc on for subsets



---

### Week 7
- improve prefetch: validate .sra files
- fastq-dump: gzipped and split forward and reverse reads
- fastq-dump very slow. Parallelized, then downloaded scATAC-seq files (SRP149536)
- QC.sh: check quality of scATAC-seq files
- Focus on ATAC-seq: scRNA-seq also shows mutations due to RNA-editing
- 

Aims
- filter for quality - trim.sh, and recheck change in quality, read length, adaptor contamination
- align subgroups with soft-clipping eg. scATAC in TF1 clones
- variant call 
- read depth and coverage
- filter alleles for MAF and ---
- nj tree and FSTs between lineages


Questions
- read length of split paired reads is 38bp - short? quality scores suggest otherwise EDIT: paper 2x38 scATAC-seq
- barcoding and demultiplexing

---

### Week 9
- align successfully
- align to whole reference genome - no splitting into mt genome and nuc
- correct output bam files
- qualimap - quality of alignment

Aims
- understand qualimap output
- variant call for only mitochondrial region (-r chrM)
- bcftools stats and plot for vcf files - use to filter for quality 

Questions
- quality filtering - 10 for read alignment, variant quality level? Use plots

---

### Week 10 
- check SRP149534 (TF1 bulk ATAC-seq) qualities: mix of read lengths: only some sra files converted to fastq.gz using `fastq-dump --split-files` 
- validate .sra files, delete fastq.gz files and re-convert all in group to __paired__ .fastq.gz files
- re-ran fastqc on each file and multiqc: still seeminly random mix of read length and
- Combined (for each of 69 SRRs): sample lineage info, generation info, read length and no., to metadata.
- Compared to processed read data provided by Ludwig et al.

Aims (look at notes from last meeting)
- minimum read depth for high confidence 
- format log
- ~~save script for RNA mapping and post-quality checks~~
- ~~re-`fastq-dump` all files~~
- check re-dumped files in group for quality: notable features/differences?- __why__
- remove nextera adaptor sequences contamination
- set up ssh key on laptop and make notes
- SRR7245916 missing - dump and include, redownloaded SRR_Acc_list.txt as 5 SRRs missing from list
- mt circularity adustment: GSNAP?

Questions
1. Duplicate removal - if the eventual goal is to estimate mutation rate and replication rate then allele frequency/mutational load is important, but that is affected by duplicate reads: duplicates should be removed?
2. Minimum no. reads acceptable for high confidence variant calling - depends on mitochondrial coverage?
3. SNP missed by ATAC-seq because of allele specific enrichment?
4. Adaptor content plot


---

### Week 11
- multiqc with missing SRR7245916
- align with bowtie2 - no splicing
- streamline samtools piping instead of creating files
- coverage: qualimap and samtools

Aims
- ~~align~~
- variant call
- circle map?
- update reference genome
- 16 pages - what's included
- Latex?
- figures - __ask for help__
- draft plan and figures
- mutation level with time for some variants
- heteroplasmy - polyploid variant caller
- not very sensitive - default
- (array of jobs)
- overall plan for tomorrow
- 10 alignment mapq? look at default
- Visualisation https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6728899/
- add QC for SRR7245916 to table
- concordant and discordant alignment meaning and impact

Questions
1. Demonstration of reproducibility
2. "using BOWTIE2 using the parameter –X2000 allowing fragments of up to 2 kb to align" -X 2000 allows a gap between f anr reads of 2000 bp inclusive? why
3. haplogrep / haplocheck


---

### Week 12
- realign using --sensitive, default quality filter, and no filter on sam flags
- research variant calling, plotting software
- plan figures

## links
- https://www.protocols.io/view/low-frequency-variant-calling-from-high-quality-mt-nfkdbkw guide for mt NGS analysis including circularization adjustment with commands
- (https://www.scielo.br/j/zool/a/FjqvPKvPcspcCydQTgV5Fdx/?lang=en popgen analysis of mt NGS data, familar software)
- https://mseqdr.org/mitobox.php analylsis software

Questions
1. We only want variants and changes in heteroplasmy in relation to the bulk generation 0 cells, rather than all differences from the rCRS. Easier assemble bulk gen 0 to rCRS, then align samples to bulk gen 0 assembly, than filtering variants afterwards?
2. Variant selection - they found 44 variants, BQ cutoff 23.8. allele frequencies > 2.5% in at least one sample.
3. Local realignment around indels
4. strand bias and read placement bias filtering, 
5. BQS recalibration?
6. circular mapper https://circularmapper.readthedocs.io/en/latest/contents/userguide.html
7. 

---

### Week 13
- variant call
- annotate variants
- preliminary individual plots
- coverages with samtools depth

Aims
- alignment stats
- variant calling stats
- interpretations
- stats and assess variant calling. No 
- duplicate read removal, BQS recalibration, and __realignment around indels__
- group plots by lineage - patterns
- summary of variants: 
- heatmap: all variants, all samples?
- ~~coverage~~
- normalise multli allelic sites
- blacklisted regions and coverage
- automate plots
- remove cap on depth
- reproduce Conor's plots
- compare to Ludwig paper variants
- percent reads that 
- limits to low level calling eg bias coverage, other sources of noise
- not interested in homoplasmy
- pos pos+1
- noise -like dynamic
- rnorm 10000 normally distibuted

Questions
1. filter when calculating coverage? since alignment was not filtered
2. max depth calculated is 8000 default - rerun with no limit



---

### Week 14
Priorities:
1. ~~save plots~~ ish --> doesn't work on hpc
2. ~~coverages~~
3. alignment stats: percent aligned. rest mostly covered by samtools depth/coverage?
    1. no. and percent reads aligned per sample: extract from alignment_summary.txt concordant, discordant, how many of either - which is important.
    3. mean coverage: per sample, per position, any of note?
    4. mean mapping quality and base quality: per sample, per position. Any of note? 
4. variant stats: what's useful?
5. filter for only interesting variants and plot DEFINE interesting: 
    1. __PASS filter for blacklisted regions and strand bias 
    2. heteroplasmic/low level (mutserve \"type\" 2: comprehensive filtering from mutserve mtDNA-Server see "Next-generation sequencing data analysis of human mitochondrial DNA in the cloud" (Weissensteiner et al., 2016).  
    3. What range within heteroplasmic variants is most useful? <- filter further 
    4. Visual stochastic drift in general (random only if there is no selection: only applies to neutral variants?) - look at single position mutation load plots. 
    5. Visual difference from low level (? what mean to set) normal (?) distribution expected from sequencing noise__ <- should be adjusted for by complex mutserve filters, comment for discussion?                             
6. ~~all lineages~~
7. DRAFT ALREADY


Aims
- coverages
- Improve plots: better scale on x axis, plot subsets of SRR_data_table eg. only PASS filter, mutations of interest; split into neutral vs missense...
- make combined variant table: filter variants by eg. starting heteroplasmy level <- research suitable noise threshold DON'T mark all variants at that position 
- Find and choose interesting variants
- put plots and R?? on windows partition of laptop for meeting
- (plot threshold lines on mut_load profiles of specific positions)
- how does mutserve define heteroplasmic/low level? sufficient coverage + variant level < 0.990?
- rerun variant calling with lowere base quality filter: depending on coverages mean quality filter and coverage plots
- add commands to analyse.sh to copy stout from slurm-id to alignement_stats.txt

Questions
1. what counts as heteroplasmic (noise threshold from variant lvl 0 or 1)
2. what should be in the git repo ie. any results (eg. coverages, plots) or just scripts to get results?




---


Week 16

- add coverage to background of plots


Aims
- adjust scales and grids on plots, especially log breaks
- save post-alignment plots
- background for individual position plots etc
- TABLES:
- prealignment: dups, reads, 
- post-alignment: reads, coverage, map and base qualities
- variant stats


Questions
1. ~~Why are is the mean depth calculated from depths.txt different from samtools coverage meandepth???~~
2. did not filter for only uniquely mapping reads - false positives? Or would mapping quality be much lower? samtools view -F 260 <- to exclude secondary alignement and unmapped reads. (flags 256 + 4)


---

Week 17
- Filter variants so HOMOPLASMIC variants (variant level is > 0.990 (default filter as homoplasmic or type 1 in mutserve annotations) ) excluded -> INTERESTING
- improve colours, fixed axis scales 
- Normalised vcf variants.

Aims
- Edit lineage tree image to use SRR names
- COMPARE TO LUDWIG SNPs
- Look at HET_LOWLVL_nofilt carefully - all lineage paths - plot mut load plots for potential positions

Questions
1. Since aim is to build a framework for further exploration of the data - is a section in my discussion to discuss the scripts I've adapted/produced appropriate? - eg. R was a bad choice because.. discuss automation of and completeness of script, eg. some scripts don't check that files/programs are present, tracebacks are practically non-existent, storage waste etc. Functionality of pipeline. context of general standards and for Conor's future use and analyses. 
2. README.md very incomplete
3. MAIN PROBLEM: STRAND BIAS. Can we assume that SNPs at the same position in different samples but some pass and some are filtered for strand bias, should not be filtered for strand bias? Could the strand bias be explained by duplicatated reads - RESEARCH AND UNDERSTAND BOTH BETTER
3. Removed variants where reference was "N" (insertion in cell line or rather deletion in rCRS) - either homoplasmic or multiallelic, and mutserve file format complicates plotting as only major variant level is used. Time - simpler to remove (only ~2 positions removed) 
4. Has keeping duplicates skewed strand bias - Conor's plot at 1492 is filtered out due to strand bias in my data (1493).
5. Ludwig report heteroplasmy as square root of allele frequency - why
6. ~~Why our positions vary relative to Ludwig so much sometimes 0, sometimes +1, +2? Ludwig_positions_convertions_and_lineages.csv~~ https://haplogrep.i-med.ac.at/2014/09/08/rcrs-vs-rsrs-vs-hg19/  <- differences between hg19 (Ludwig ref) and rCRS (our ref)
7. Terminology best practice: use "our" data? or "this study's" data
8. How to represent genotypes - total number, and number in each lineage. eg. __ heteroplasmic SNPs total, 35 in G11, 17 in B3, etc
9. terminology: __point mutations__, SNPs, variants?  
10. bulk variant calling inconsistency - what does this say about our pipeline when Ludwig's variant calling was reliable?
11. Statistcal difference between ours and Ludwig's allele freqs.
12. In general there are lots of very low level point mutations in the first generation: worse sequencing earlier on? Are these useful? Low confidence as only present in 1st gen/1st gen and one of the bulks/1st gen and both bulks?
13. Also, first generation jump in allele freq - single cell from bulk popultion cause this in dramatic bottlneck?


---


### Week 18



Aims
- get per strand per allele no. reads
- add coverage to mut_load_plots
- write pre-variant calling draft
- update readme
- reproduce 2x bulk replicate plot and correlation
- correlation with Ludiwig allele freqs
- rerun variant calling with no/lower limit minimum allele frequency: dips in allele frequency changes to zero may be due to anything <1 not output by mutserve and plotted as zero. Oddly why are some allele frequencies <0.01 in the data? Consistent timecourse apart from below 1% in eg. 5007 in G11_longest_upper
- update README.md

Questions

### Week 20 
RESULTS TODO 28th
- Alignment stats: FINISH figure and results. Ask about which plots for coverage representation.
- Variant calling: exactly what figure and table of point mutations in results. Fill out with current data, even though filtering may change. Improve mut_load_plots with coverage and colour, re-label bulks. Panel of example mut_load_plots
- Comparison: correlation per sample and per position: bullet points for disussion eg. why is each position consistently different to Ludwigs? - no. duplicates across position? Difference in variant calling filters eg. baseq.

- STRAND BIAS: can we validate mutations, but not allele frequencies

- reran with baseq 24: threshold of allele frequencies that fall on the Ludwig distribution which are not potentially sequencing errors. allele frequencies
- Mutserve exludes alleles with less than 3 reads per strand!!!
-  


Questions
1. What kind of pattern in the mutation load profiles are we looking for, how to describe? What should the drift look like eg. what is the maximum likely change in allele frequency from one to generation to another (research).
2. Spearman's rank for correlation? - non-parametric, monotonic (differences in our variant calling pipleine means relationshipbetween our AFs and theirs aren't directly linear, see figure. Per position should be directly linear - Pearson's rank instead?)
per sample, per mutation position - what is statistically informative? worth looking at correlation between number of duplicates at the position and each position's (of the 40) deviation from Ludwig?
3. which coverage representation figure
4. The variant calling pipeline of this framework shows lower sensitivity to point mutations(???) 40 out of 44, but removed fewer false negatives due to validation across the lineage. Difference in variant calling parameters. See discussion.
5. ~~Why sqrt allele frequency for "heteroplasmy"?~~ To emphasise AF at low levels
6. Why did the the longer reads have a lower mapping quality???
7. Can I add a read only link to google sheets table as a table for the supplementary material?
8. Ludwig's duplication rates
9. How to refer to our data vs Ludwig's. cite every time? "Our"/ "this framework/pipleine", vs "their", "Ludwig et al.'s", previous analysis/ previously identified.








# Resume log: week 4 of temp research assistant
## Wed 6th April
- [ ] refactor mutation plot script:
  - [x] section titles and spacing
  - [x] use a variable to choose vcf folder and alignment\_stats/ files
  - [ ] use a variable to choose results/ folder
- [x] find cause of odd G2 lable in G3 in plots - easiest to check by lineage
- [ ] then fix: if ^ not solved reconstruct SRR\_lineage\_generation and check lineage\_paths.txt SRRs
- [ ] remove na position in lin\_mut\_load\_change\_lin\_val (?) - lineage mut plots
- [ ] COMPARE effect of dups, and consensus sequence on STRAND BIAS 
- [ ] X why in supposedly validated mutation load lineage plots in some samples the AF doesn't reach threshold
- [x] Improve lineage validation- if mutation occurs more than once in the lineage path, rather than whole lineage eg. B3\_shortest\_upper instead of over all B3. 
- [x] Removed bulks and mix as lineages from lineage\_paths.txt 
- [x] add A9 lineage to all analyses

## Thurs 7th April
- [x] why in supposedly validated mutation load lineage plots in some samples the AF doesn't reach threshold
- [x] fixed VariantLevel selection when plotting mutation load plots!!! Many fewer missing calls
- [x] fix Bulks and A9 in all\_variants\_in\_lineage
- [ ] Modify pipeline to work outside of rocket.hpc and slurm on a typical linux system
  - [x] install software script
    - [x] samtools
    - [x] fastqc
  - [x] get SLURM module versions
  - [x] prefetch_files.sh
  - [ ] QC.sh
  - [x] post_alignment.sh
  - [ ] analyse.sh
- [ ] het in consensus lower maybe 0.90
- [x] Modify SraRunTable.txt for easier selection of runs - add keyword column
- [ ] Put individual pos plots into folder or one pdf
- [ ] extract bcftools calls and depths
- [ ] fix G2 label in G3

NOTES
- Positions 16172,4,5,6,7 close together suspicious
- G3 spikes which ones

## Friday 8th April
 - [ ] Modify to work outside of slurm
   - [ ] Downloads
     - [x] bowtie2
     - [x] multiqc
     - [ ] sratools - ubuntu vs redhat/centos
     - [x] mutserve
     - [x] bcftools
   - [x] QC.sh
   - [x] analyse.sh
   - [ ] variant\_call.sh
   - [ ] plot\_mutations.sh
 - [x] Successfully run test up to post\_alignment.sh on my laptop
 - [ ] Meet and run with Jordan

## Monday 11th and Tuesday 12th April
 - [x] save all individual plots into multi-page PDFs
 - [ ] Variant call again with mapping quality higher
 - [x] Improve output of pileup: delete no alt/indel positions, output F,R,total reads for ref and alt alleles
 - [x] Read pileups into R for small test group
 - [x] Reformate pileups into mutserve-like table
 - [ ] Run Rscript for test group
 - [x] re-run pileups for all samples

## Wednesday 13th April
 - [x] Research building own variant caller
 - [x] Research and start argument parsing in bash scripts 

## Wednesday 20th April
 - [x] Variant analysis usable with massive mpileups - variant stats and lineage validation
 - [x] Run whole script with mpileups
 - [x] Create vectors
 - [ ] Weight AFs?
 - [ ] Check plots are not affected by multiallelic sites
 - [x] implement argument parsing in analyse.sh
 - [ ] Check 54 from cnodups against Ludwig
 - [ ] look at bulk replicates plot for pileups
 - [ ] Test for no correlation between lineages? eg large spikes.

## Thursday 21st April
 - [ ] Compare mapping quality and other metrics of third generations against the rest - make table when making plots. 
   - [ ] Discount alleles where depth exceeds twice the mean depth usually caused by artifacts
   - [ ] test mean of g3 against the rest. Use difference to filter.
 - [ ] apply BAQ when making pileup
 - [ ] apply --ambig-reads drop when making pileups: reads that don't span all of a short tandem repeat region
 - [ ] add name of allele to individual plots eg. actg del
 - [ ] why homoplasmic alleles still present

## Friday 22nd April
 - [x] Restructure vectors
 - [x] Research bcftools options, including BAQ and --ambig-reads
 - [x] Compare mapping quality and other metrics of third generations against the rest - make table when making plots. 
 - [x] apply BAQ when making pileup
 - [x] apply --ambig-reads drop when making pileups: reads that don't span all of a short tandem repeat region
 
## Monday 25th April
 - [x] Correct extraction of allele depths from allele positions caused by multiallelics
 - [ ] Exclude gen 3 and compare
 - [ ] Correct plotting for multiallelics
 - [ ] Compare results of BAQ and --ambig-reads
 - [x] Get full pilup with all positions not just mutations to add all reference depths to vectors for Jordan
 - [ ] Send vectors to Jordan

## Tuesday 26th April
 - [x] **Send vectors**
 - [ ] **Call IT desk**
 - [ ] Duplicate reads excluded?
 - [ ] Check SRR number in lineage - direct not doubly indirect on Ludwig tree

## Wednesday 27th April
 - [x] Remove Gen 3
 - [x] Correct plotting for multiallelics
 - [x] Induction meeting
 - [x] samtools view exclude duplicates for pileups
 - [ ] compare effect excluding duplicates with samtools view on pileups

## Thursday 28th April 
- [ ] Email Conor
- [ ] Research regions of suspected misalignments
- [x] Increase mapping q threshold to see effect on regions
  - introduced more "validate positions
- [x] Reresearch mapping quality thresholds for bowtie2
  - http://biofinysics.blogspot.com/2014/05/how-does-bowtie2-assign-mapq-scores.html
  - MAPQ >= 23 allows up to 3 mismatches > bq20, (higher the base quality, more impact on mapq: more trustworthy. 
- [x] Lower mapq to 4: only "uniquely" mapping reads (what are multi and maxi), up to 5 mismatches with bq >20
- [x] compare effect of excluding duplicates with samtools view on pileups
  - (No Difference - bcftools mpileups excludes dups by default)
- [x] Check SRR number in lineage - direct not doubly indirect on Ludwig tree
  - [x] save lineage tree with Snumbs into images/
- [ ] why homoplasmic alleles still present
- [x] base quality 23.8 like Ludwig - see effect on bulks

## Friday 29th April
- [x] Find further differences between mq4 and mq 18:
 - Too many reads removed: < 3 mismatches on one read causes read removal, including mutations supporting real mutations: eg. 4214 4446 5007 5862
  - Mutations missed entirely by lin val: 6055 6025
- [x] Add ratios column: find positions with =< 2x mean coverage
  - Gen 3 spikes coverages 
- [ ] strand bias calculation
- [x] Difference in BAQ
- [ ] **CALL IT DESK DON'T FORGET**
- [ ] View suspected misalignment regions
  - 3107 
  - 4548: 2 values per gen for all
- [ ] Add gen 3 back for mq4 - any improvement?
- [ ] Check/rerun alignment to consensus for homoplasmic alleles
### Suggestions:
 - 0.01 threshold not so useful: min no. reads in at least 1 sample: eg. 4552
 - Mapping Quality much lower, based on how bowtie2 assigns mapping quality, and improved plots.
 - baseq from Ludwig's gaussian distribution threshold for now. 

## Tuesday 3rd May
- [ ] Mark potential sequencing errors
- [ ] Strand bias score calculation
- [x] Further research on expected sites: detect errors through autocorrelation - compare to errors at expected sites 
  - [ ] Correlation of expected sites and strand bias score
- [ ] Re-run alignment to check it's aligned to parent\_consensus - homoplasmic sites shouldn't be present
- [x] copy chrM of parent\_consensus to laptop and view suspected misalignment regions
  - 3107 
  - 4548: 2 values per gen for all
- [ ] Ludwig - technical replicate or biological?
  - "two technical replicates of the bulk TF1 sample"
- [ ] Ludwig bulk parent population clone?
  - "parental colony" then refer to clones and sub-clones
  - Also referred to as the bulk TF1 sample
- [x] replicate plot on log scale, smaller dots, 
  - [ ] overlay Ludwig and 
  - [x] my selected positions. 
- [ ] angela.pyle@ncl.ac.uk <- ticket IT:wq

## Wednesday 4th May
- [x] research other potential datasets
- [x] strand bias calculation
- [x] email IT desk, open new ticket

## Friday 6th May
- [x] Weighted allele frequency (redundant)
- [x] Strand bias filter: 0.3 < fwd/(rv+fwd) <0.7
- [x] Overlay Ludwig's positions
- [ ] Overlay our AFs for Ludwig positions

## Monday 9th May
- [x] Overlay our AFs for Ludwig positions
- [x] Overlay Ludwig's AFs with ours.
- [ ] Find cause of weird 4548 plot
  -> insertion in one sample caused misalignment and false calls from 4541:4555

## Tuesday 10th May 
- [x] Write email to angela about IT ticket
- [ ] No mismatches per read - 35 bp reads with 4 mismatches worse than 75bp reads with 4. Change mapping quality with read length?
  - [ ] increase mapping quality to 1 less mismatch for 75bp reads: peaks in Gen 2 and 3 (only some in Gen 2 and 3 have longer reads) 
  - [ ] Look at Gen3 read lengths and base q
- [x] Fix correlation with Ludwig AFs for larger data, multiallelics
- [x] Recreate Ludwigs correlation plot from their data
  - NOT IDENTICAL TO THEIR PUBLISHED PLOT?
- [x] Plot just bulks correlation
  - [ ] Repeat for 
- [ ] Plot with strand bias sites removed

## Wednesday 11th May
- [ ] STRAND BIAS
  - [x] Plot mut load plots with strand bias sites removed
  - [ ] Research flexible strand bias removal 
    - [ ] look at the strand bias of the same site in different samples?
    - [ ] measures of strand bias eg. bcftools, GATK
    - occassionally removes sites in generations with genuine mutation eg.4214
    - occassionally doesn't remove sequencing error sites (
    - detect misalignment?
  - [ ] Plot correlation plots with strand bias sites removed
  - Implement unequal strand weight check/unequal strand AF - whichever is not atm
- [x] Find reason for AF difference between plots: STRAND\_BIAS sites removed (PASS) vs nofilt
  - No difference in AF just some seq err peaks removed!!
- [ ] How low an AF before correlation plot blobs?
- [ ] Each strand equally weighted or by no. reads?
  - atm the strand with more coverage influences AF the most. If removing sites with low coverage/unequal strand weight
   See results

## Friday 13th May
- [x] Create vectors without supporting reads of strand bias variants
- [ ] log plots of correlations

## Monday 16th May 
- [ ] 

## Tuesday 17th May
- [ ] Run pipeline wih weighted AFs
  - [ ] compare plots
- [ ] tidy and make frd and rv AFs, types of strand bias, difference columns 

## Wednesday 18th
- [ ] implement min no reads on pileup data
- [ ] degree difference in AF depend on no supporting reads

## Thursday 19th
- [ ] check strand bias filtering for all positions that have 0 on one strand
- [x] write parse\_longest\_paths.txt for Jordan
- [x] email Jordan with lineage tree reading

## Friday 20th May
- [ ] Run pipeline wih weighted AFs
  - [ ] compare plots
- [ ] difference between weighted and unweighted, and AF strand bias 
- [ ] check strand bias filtering for all positions that have 0 on one strand
- [ ] bcftools SB values with strand difference in AF

- strand bias adjustment should incorporate no supporting reads and AF

## Thursday
- [ ] change lineage validation to if >=1 read on both strands, not >0.01 AF in at least onei
  - [ ] then >= 1 read on both strands in at least one sample
- [ ] option to change results dir
- [ ] implement GATK strand bias measure
- [ ] email Conor
  - change alpha of plots
  - change axis of plots 
  - thoughts on working remotely

## Monday 30th May
- [ ] min 1 read per strand

## Tuesday 31st May
- [x] implement GATK strand bias 
- [ ] how to implement filters into stats? 
  - raw reads vs filtered?
  - for any autocorrelated, a test for correlation between lineages would help identify sequencing error eg. 16175 in min_1_read, 309 310
- [ ] CALL IT DESK
- [ ] remove lineage labels from plots which aren't in lineage_paths.txt
- [ ] Proportion test
- [x] View results of all SB methods


## Monday 6th June
- [ ] Proportion test vs others
- [ ] Message Jordan
- [x] actually call the IT desk
- [ ] g score between autocorrelation of all vs autocorrelation of strand bias calls
- [ ] write methods for Jordan

## Tuesday 7th June
- [ ] compare and decide on strand bias measure
  - [ ] strand_bias_ratio_AF_SqrDiff
  - [ ] how each works
Questions
- How to implement strand bias with autocorrelation stats

## Wednesday 8th June
- [ ] get missing values for bcftools
- [ ] error log $>
- [ ] add checks the file is present
