

  ## NOTE: this script should be called via the bash script "plot_mutations.sh" 
   # Or working directory should be set manually.

# set working directory
#args <- commandArgs(trailingOnly = T)
#print(args)
#setwd(args[1])
setwd("/home/thomas/Documents/projects/Research_proj/Ludwig_2019/")
#setwd("/home/thomas/Documents/projects/Research_proj/pipeline_test/Ludwig_2019/")


  ## Load packages ##
# Check if packages are installed, install to .R_local_lib in working directory if needed.
#local_lib_path <- paste0(args[1],"/.R_local_lib/")
#print(local_lib_path)
#.libPaths(c(local_lib_path, .libPaths()))

packages <- c("tidyr","tidyverse","ggplot2","gridExtra","ggrepel","egg","grid","BiocManager", "circlize", "reshape2","cowplot", "data.table", "dplyr")
lapply(packages, FUN = function(i) {
  if (!require(i, character.only = TRUE)) {
    install.packages(i, dependencies = TRUE, lib = local_lib_path, repos="https://www.stats.bris.ac.uk/R/")
    library(i, character.only = TRUE)
    }
  }
)

if (!require("ComplexHeatmap")){
  BiocManager::install("ComplexHeatmap")
  library("ComplexHeatmap")
}
library("ComplexHeatmap")


  #########################  Options  ##########################

# Change boolian to choose
use_pileups <- TRUE
create_vectors <- FALSE
pre_plots <- FALSE
post_plots <- TRUE
exploratory_plots <- FALSE
Ludwig_comparison <- FALSE
position_specific_plots <- TRUE
print(paste0("Make large exploratory plots: ", exploratory_plots))

  ###########################  DATA  ############################
  # Controls for which vcf/ directory and which post-alignment files in alignment_stats/ 
  # are imported. Allows simpler comparison between datasets eg. calls from reads aligned 
  # to rCRS or a consensus sequence of parent clones, duplicates or removed duplicates etc.
# Change string to choose:
# append_string: files made using post_alignment_QC.sh will have the same string attached to the end of the following files (alignment_stats/):
# "depths", "depths_qfilt", "all_coverages", "all_coverages_qfilt", "mean_coverage", "mean_coverage_qfilt"
# the "append_string", then ".txt". This is typically the bam dir name eg. "bam_consensus-nodups".
# Make sure to add a preceding "_". eg. "_bam_c"
append_string <- "_bam_cnodups"
vcfdir <- paste0("vcf", append_string)
bcfdir <- paste0("mpileups", append_string, "_mq4_bq23.8")
group_name <- "SRP149534"


  ####################  Read SRR files  ########################
filenames <- list.files(vcfdir, pattern="*_annotated.txt")
filepath <- paste0("data/group_", group_name, "_SRRs.txt")
SRR_names <- as.list(read.table(filepath, stringsAsFactors = FALSE))[1]
SRR_names <- SRR_names$V1
# Create list of data frame names without the ".txt" part 
#SRR_names <-substr(SRR_names,1,10)
#SRR_names <- c("SRR7245880", "SRR7245881","SRR7245883", "SRR7245897", "SRR7245917")



     #############  BCFtools  pileups  ###################
 # First compare AFs of variants called with current pipeline: ie. extract AF and allele depths (F and R) of all positions in the lineage called by mutserve.
# all_variants in lineage
if (use_pileups == TRUE){
  bcf_mpileups <- list()  # All 

  #######  Vectors of allele read proportions ######
  # for each nucleotide ACGT (solves multiallelic problem of more than one AF per genomic position)
  # (weighted?) supporting reads over total reads
  
# Create function. (Vectors themselves created when making bcf_SRR_table_lists below, so sites with no variants can be removed from large dfs)
AF_vectors <- function(SRR_table, base){
  #mt_vectors <- data.frame()
  #colnames(mt_vectors) <- c("A","C","G","T")
  #for (base in c("A","C","G","T")){
  mt_vector <- rep(0,16569)
  for (pos in SRR_table$Pos[which(SRR_table$Variant == base)]){
    #print(paste0("Pos = ",pos))
    #print(SRR_table$Pos[which(SRR_table$Pos == pos)])
    mt_vector[pos] <- as.numeric(SRR_table[SRR_table$Variant == base & SRR_table$Pos == pos, "Variant_AD"])
  }
  for (pos in SRR_table$Pos[which(SRR_table$Ref == base)]){
    #print(paste0("ref pos = ", pos))
    #print(SRR_table$Pos[which(SRR_table$Pos == pos)])
    mt_vector[pos] <- as.numeric(SRR_table[SRR_table$Ref == base & SRR_table$Pos == pos, "Ref_AD"][1])
  }
  #mt_vectors[[base]] <- mt_vector
  return(mt_vector)
}


# Read in mpileups:
for(i in SRR_names){
  filepath <- file.path(bcfdir,paste0(i,"_mpileup.vcf"))
  bcf_mpileups[[i]] <- read.table(filepath, sep = "\t", header = F, stringsAsFactors = T, comment.char = "#") 
  colnames(bcf_mpileups[[i]]) <- c("CHROM", "Pos", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "V10")
  bcf_mpileups[[i]] <- bcf_mpileups[[i]] %>% separate(col = "V10", into = c("Phred-scaled-GT-likelihoods", "STRAND_BIAS_pval", "ADF", "ADR", "AD"), sep = ":")
  bcf_mpileups[[i]] <- bcf_mpileups[[i]] %>% separate(col = "ADF", into = c("ref_ADF","alt_ADF"), sep = ",")
  bcf_mpileups[[i]] <- bcf_mpileups[[i]] %>% separate(col = "ADR", into = c("ref_ADR","alt_ADR"), sep = ",")
  bcf_mpileups[[i]] <- bcf_mpileups[[i]] %>% separate(col = "AD", into = c("ref_AD","alt_AD"), sep = ",")
  bcf_mpileups[[i]] <- bcf_mpileups[[i]][,c(1:11,16,12,14,17,13,15)]
  bcf_mpileups[[i]]$Depth <- as.numeric(bcf_mpileups[[i]]$ref_AD) + as.numeric(bcf_mpileups[[i]]$alt_AD)
  #bcf_mpileups[[i]] <- bcf_mpileups[[i]] %>% separate(col = "INFO", into = c("raw_DP", "I16", "QS", "VDB", "SGB", "RDB"), sep = ";")
  #print(range(bcf_mpileups[[i]]$STRAND_BIAS_pval)
    
}
  
# Convert mpileups df into same format as Mutserve vcfs (SRR_table_list), and calculate VariantLevels, Coverages etc.
# AND create vectors of reads supporting each base for every position

# vcf format of bcftools mpileup (tab separated):
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SRR7245880.bam
# mutserve format (tab separated)
#ID      Filter  Pos     Ref     Variant VariantLevel    MajorBase       MajorLevel      MinorBase       MinorLevel      Coverage            CoverageFWD     CoverageREV     Type

bcf_SRR_table_list <- list()
vector_list <- list()
for (i in SRR_names){
  bcf_SRR_table_list[[i]] <- data.frame(bcf_mpileups[[i]]$Pos)
  colnames(bcf_SRR_table_list[[i]]) <- "Pos"
  if (bcf_mpileups[[i]]$STRAND_BIAS_pval > (40+bcf_mpileups[[i]]$Depth/2)){
    bcf_SRR_table_list[[i]]$Filter <- "STRAND_BIAS"
  }
  else{
    bcf_SRR_table_list[[i]]$Filter <- "PASS"
  }
  bcf_SRR_table_list[[i]]$Ref <- bcf_mpileups[[i]]$REF
  bcf_SRR_table_list[[i]]$Variant <- bcf_mpileups[[i]]$ALT
  bcf_SRR_table_list[[i]]$Variant_AD <- bcf_mpileups[[i]]$alt_AD
  bcf_SRR_table_list[[i]]$Variant_ADF <- bcf_mpileups[[i]]$alt_ADF
  bcf_SRR_table_list[[i]]$Variant_ADR <- bcf_mpileups[[i]]$alt_ADR
  bcf_SRR_table_list[[i]]$Ref_AD <- bcf_mpileups[[i]]$ref_AD
  bcf_SRR_table_list[[i]]$Ref_ADF <- bcf_mpileups[[i]]$ref_ADF
  bcf_SRR_table_list[[i]]$Ref_ADR <- bcf_mpileups[[i]]$ref_ADR
  bcf_SRR_table_list[[i]]$Coverage <- as.numeric(bcf_mpileups[[i]]$ref_AD) + as.numeric(bcf_mpileups[[i]]$alt_AD)
  
  # Correct multiallelic allele depths (Coverage): sum of multiallelic allele depths plus ref allele depth
  for (pos in bcf_SRR_table_list[[i]]$Pos[duplicated(bcf_SRR_table_list[[i]]$Pos)]){
    bcf_SRR_table_list[[i]]$Coverage[which(bcf_SRR_table_list[[i]]$Pos == pos)] <- sum(as.numeric(c(bcf_SRR_table_list[[i]]$Variant_AD[which(bcf_SRR_table_list[[i]]$Pos == pos)], bcf_SRR_table_list[[i]]$Ref_AD[which(bcf_SRR_table_list[[i]]$Pos == pos)][1])))
  }
  mean_cov <- mean(bcf_SRR_table_list[[i]]$Coverage)
  bcf_SRR_table_list[[SRR_name]]$meanCovRatio <- as.numeric(as.list(bcf_SRR_table_list[[SRR_name]]$Coverage))/mean_cov
  bcf_SRR_table_list[[i]]$VariantLevel <- as.numeric(bcf_mpileups[[i]]$alt_AD) / as.numeric(bcf_SRR_table_list[[i]]$Coverage)
  bcf_SRR_table_list[[i]]$RefLevel <- as.numeric(bcf_mpileups[[i]]$ref_AD) / as.numeric(bcf_SRR_table_list[[i]]$Coverage)
  bcf_SRR_table_list[[i]]$Type <- 2
  #apply(bcf_mpileups[[i]], 1, FUN = bcf_mpileups[[i]]["ref_AD"]+bcf_mpileups[[i]]["alt_AD"])
  #    as.data.frame(as.numeric(bcf_mpileups[[i]]$ref_AD) + as.numeric(bcf_mpileups[[i]]$alt_AD))
  
  # Create vectors, then remove massive number of variants with no supporting reads (all <*> sites)
  
  if (create_vectors == TRUE) {
    vector_list[[i]] <- data.frame("As"=AF_vectors(bcf_SRR_table_list[[i]], "A"))
    vector_list[[i]]$Cs <- AF_vectors(bcf_SRR_table_list[[i]], "C")
    vector_list[[i]]$Gs <- AF_vectors(bcf_SRR_table_list[[i]], "G")
    vector_list[[i]]$Ts <- AF_vectors(bcf_SRR_table_list[[i]], "T")
    #vector_list[[deparse(substitute(SRR))]] <- data.frame(As=vectorsA, Cs=vectorsC, Gs=vectorsG, Ts=vectorsT)
    filestring <- paste0("results/", i,"_reads_vector.csv")
    write.csv(vector_list[[i]], file = filestring, quote = F)
  }
  bcf_SRR_table_list[[i]] <- bcf_SRR_table_list[[i]][which(bcf_SRR_table_list[[i]]$Variant_AD >= 1 ),]
  
  bcf_mpileups[[i]] <- NULL

}

}  # end if(use_pileups == T)


    ###############  Read additional files  #####################

  ## Read coverage files ##
file_string <- paste0("alignment_stats/depths_qfilt", append_string, ".txt")  # temporarily filtered depths file as well - depths_qfilt.
depths <- read.table(file_string, sep = "\t", header = F, stringsAsFactors = T)
file_string <- paste0("alignment_stats/depths_qfilt", append_string, ".txt")
depths_qfilt <- read.table(file_string, sep = "\t", header = F, stringsAsFactors = T)
colnames(depths_qfilt) <- c("chr", "Pos", SRR_names)
colnames(depths) <- c("chr", "Pos", SRR_names)

  ## Read lineage path/s ##
# path through a lineage specified and read from 'lineage_paths.txt'
paths <- list()
paths <- as.list(strsplit(readLines("lineage_paths.txt"), " "))

# Name, lineage and generation for each sample (for plotting colours and annotations)
lineages <- c("bulk", "bulk", "A9","B11","C7","D3","F4","G11","B3","B5","B9","C4","C10","D2","C9","D6","G10","B11","B11","B5","B5","F4","F4","A9","A9","B3","B3","D2","D2","G11","G11","B3","B3","D2","D2","G11","G11","B11","B11","B5","B5","F4","F4","B3","B3","B3","B3","D2","D2","G11","G11","G11","G11","B3","B3","B3","B3","G11","G11","G11","G11","G11","G11","G11","G11","G11","G11","mix","mix")
generation <- c("0","0","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","2","2","2","2","2","2","2","2","2","2","2","2","2","2","3","3","3","3","3","3","3","3","3","3","3","3","4","4","4","4","4","4","4","4","4","4","5","5","5","5","5","5","5","5","6","6","7","7","8","8","","")
generation_axis_labs <- c("bulk","bulk","G1","G1","G1","G1","G1","G1","G1","G1","G1","G1","G1","G1","G1","G1","G1","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G2","G3","G3","G3","G3","G3","G3","G3","G3","G3","G3","G3","G3","G4","G4","G4","G4","G4","G4","G4","G4","G4","G4","G5","G5","G5","G5","G5","G5","G5","G5","G6","G6","G7","G7","G8","G8","mix","mix")
SRR_lineage_generation <- data.frame(SRR_names,lineages,generation,generation_axis_labs)

lineage_cols <- c("bulk"="royalblue4", "G11"="magenta3", "B3"="orange", "D2"="yellow", "F4"="grey", "B5"="burlywood4", "B11"="palevioletred2", "A9"="lightseagreen", "D3"="palegreen1", "C7"="lightgoldenrod", "C4"="pink2","C10"="cyan", "B9"="plum1","G10"="steelblue2", "D6"="springgreen4","C9"="red", "mix"="darkgreen")

# table of all sequencing runs in Ludwig paper - only TF1 bulk ATAc-seq needed
raw_sample_info <- read.csv("data/SraRunTable_1.csv", header = T)

# raw sequencing information
pre_multiqc <- read.table("multiQC/group_SRP149534_multiQC_report_data/multiqc_general_stats.txt", header = T) 
colnames(pre_multiqc) <- c("SRR_sample", "percent_dup", "percent_gc", "sequence_lengths", "percent_fails", "num_seqs")


  ######################   Pre-alignment   #########################
  #####################  Plots and tables  #########################

if (pre_plots == TRUE) {

median.default(pre_multiqc$num_seqs)
min(pre_multiqc$num_seqs)
max(pre_multiqc$num_seqs)
median(pre_multiqc$percent_dup)
min(pre_multiqc$percent_dup)
max(pre_multiqc$percent_dup)
mean(pre_multiqc$percent_gc)
sd(pre_multiqc$percent_gc)

# see post-alignment for meanbaseq
# histogram of estimated number of duplicated sequences in raw data
pre_dup_hist <- ggplot(data = pre_multiqc, aes(percent_dup)) +
  geom_histogram(fill = "dodgerblue3", colour = "black", binwidth = 2) +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  labs(x ="% duplicate reads in clone") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=13)) +
  geom_vline(xintercept = median(pre_multiqc$percent_dup), linetype = "dotted", size = 1.5, colour = "lightblue")
file_string <- "results/pre_alignment_dup_seqs_hist.png"
ggsave(file=file_string, plot=pre_dup_hist, width = 8, height = 4, units = "in")

# histogram of raw reads
pre_num_seqs_hist <- ggplot(data = pre_multiqc, aes(num_seqs)) +
  geom_histogram(fill = "dodgerblue3", colour = "black", bins = 30) +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  scale_x_continuous(breaks = c(5000000,10000000,15000000,20000000,25000000), labels = c("5M","10M","15M","20M","25M")) +
  labs(x ="Number of reads in clone") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=13)) #+
  #annotate(geom = "text", x = 2000000, y = 16, label = "A")
file_string <- "results/pre_alignment_num_seqs_hist.png"
ggsave(file=file_string, plot=pre_num_seqs_hist, width = 8, height = 4, units = "in")

# combine for figure 2
pre_plots <- list(pre_num_seqs_hist,pre_dup_hist)
pre_plots <- ggarrange(plots = pre_plots, nrow = 2, align = "hv")
ggsave(file="results/pre_plots.png", plot = pre_plots, width = 8, height = 8, units = "in")

}  # end if (pre_plots <- TRUE)

  #########################  POST-ALIGNMENT  ############################
  #####################  Coverage plots and stats  ######################

if (post_plots == TRUE) {
## x axis as sample and as pos.
depths_max <- lapply(depths[,3:ncol(depths)], max)
depths_qfilt_max <- lapply(depths_qfilt[,3:ncol(depths_qfilt)], max)
third_y_lim_maxcoverage <- max(as.data.frame(lapply(depths[,3:ncol(depths)], max)))
third_y_lim_maxcoverage_qfilt <- max(as.data.frame(lapply(depths_qfilt[,3:ncol(depths_qfilt)], max)))


coverage_plots <- list()
for (i in SRR_names){
  #depths_qfilt_log2[[i]] <- log2(depths_qfilt[[i]])

  coverage_plots[[i]] <- ggplot() +
    geom_line(data = depths_qfilt, aes(Pos, depths_qfilt[[i]])) +
    coord_trans(y="log2") +
    scale_y_continuous(trans='log2')
}


# overall alignment rate for each clone
percent_alignment <- read.table("alignment_stats/alignment_and_duplicate_summary.txt", sep = " ", header = T, skip = 1, nrows = 69)

mean(percent_alignment$Overall_alignment_rate)
min(percent_alignment$Overall_alignment_rate)

# summary statistics for coverages (before and after filtering for mapping and base quality):
# information included (column names): "X.rname", "startpos","endpos","numreads","covbases","coverage","meandepth","meanbaseq","meanmapq","SRRs","read_lengths","calculated_mean_depth","calculated_sd_depth" and min and max depths for qfilt
file_string <- paste0("alignment_stats/all_coverages_qfilt", append_string, ".txt")
all_coverages_qfilt <- read.csv(file_string, header = T, sep = "\t")
all_coverages_qfilt$SRRs <- SRR_names

file_string <- paste0("alignment_stats/all_coverages_qfilt", append_string, ".txt")  # temporarily also filtered version
all_coverages <- read.csv(file_string, header = T, sep = "\t")
all_coverages$SRRs <- SRR_names
df_SRR_names <- data.frame(SRR_names)
all_coverages$read_lengths <- merge(data.frame(SRR_names), raw_sample_info, all.x = T, by.x = "SRR_names", by.y = "Run")[,3]
all_coverages_qfilt$read_lengths <- merge(data.frame(SRR_names), raw_sample_info, all.x = T, by.x = "SRR_names", by.y = "Run")[,3]

all_coverages$calculated_mean_depth <- as.numeric(lapply(depths[,3:ncol(depths)], mean))
all_coverages$calculated_sd_depth <- as.numeric(lapply(depths[,3:ncol(depths)], sd))
all_coverages_qfilt$calculated_mean_depth <- as.numeric(lapply(depths_qfilt[,3:ncol(depths_qfilt)], mean))
all_coverages_qfilt$calculated_sd_depth <- as.numeric(lapply(depths_qfilt[,3:ncol(depths_qfilt)], sd))
all_coverages_qfilt$min_depth <- as.numeric(lapply(depths_qfilt[,3:ncol(depths_qfilt)], min))
all_coverages_qfilt$max_depth <- as.numeric(lapply(depths_qfilt[,3:ncol(depths_qfilt)], max))

# create dataframe for boxplot of per sample coverages
all_depths_qfilt <- list()
SRR_rep_names <- list()
for (SRR in SRR_names) {
  tmp <- as.list(depths_qfilt[[SRR]])
  all_depths_qfilt <- c(all_depths_qfilt,tmp)
  SRR_rep_names <- c(SRR_rep_names, rep(SRR, nrow(depths_qfilt)))
}
depths_qfilt_bplot_data <- data.frame(matrix(nrow = 1143261, ncol = 2))
depths_qfilt_bplot_data$X1 <- SRR_rep_names
depths_qfilt_bplot_data$X2 <- as.numeric(all_depths_qfilt)

depths_qfilt_bplot_data <- reshape2::melt(depths_qfilt[,-1], id.vars = c("Pos"))
depths_qfilt_bplot_data["read_depths" <=1] <- 1
colnames(depths_qfilt_bplot_data) <- c("Pos", "SRRs", "read_depths")


fold_coverage <- (mean(pre_multiqc$sequence_lengths)/16569)*mean(all_coverages_qfilt$numreads)
fold_coverage
mean(all_coverages$meandepth)

  ## pre map and base quality filtering ##
mean(all_coverages$meandepth)
mean(all_coverages_qfilt$meandepth)
mean(all_coverages$meanbaseq)
mean(all_coverages_qfilt$meanbaseq)
mean(all_coverages_qfilt$meanmapq)
mean(all_coverages$meanmapq)

# coverage plot for all bases, all samples
yax <- c(1,2,4,8,16)
yax <- 2^yax

pos_coverage_all_plot <- ggplot(depths_qfilt_bplot_data) +
  geom_line(aes(Pos,log2(read_depths), colour = SRRs),size = 0.5)  +
  scale_y_continuous(trans='log2',breaks = c(1,2,4,8,16), labels = yax) +
  scale_x_continuous(breaks = c(0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,11000,12000,13000,14000,15000,16000))+
  #coord_trans(y="log2") +
  geom_hline(yintercept = log2(60), colour = "red", size = 1,linetype="dotted")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=13), legend.position = "none")+
  labs(x= "Mitochondrial genome position", y = "log2(Read depth)")

ggsave(file="results/genome_coverage_plot.png", plot=pos_coverage_all_plot, width = 8, height = 4, units = "in")

# filtered read histograms: coverage, no.reads, base quality, mapping quality
mean_coverage_hist_qfilt <- ggplot(data = all_coverages_qfilt, aes(meandepth)) +
  geom_histogram(fill = "dodgerblue3", colour = "black") +
  theme_bw()
mean_reads_hist_qfilt <- ggplot(data = all_coverages_qfilt, aes(numreads)) +
  geom_histogram(fill = "dodgerblue3", colour = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  scale_x_continuous(labels = scales::comma) +
  theme_bw()
mean_baseq_hist_qfilt <-  ggplot(data = all_coverages_qfilt, aes(meanbaseq)) +
  geom_histogram(fill = "dodgerblue3", colour = "black", bins = 7) +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  scale_x_continuous(labels = scales::comma, breaks = c(34.7,34.8,34.9,35)) +
  theme_bw()
mean_mapq_hist_qfilt <- ggplot(data = all_coverages_qfilt, aes(meanmapq)) +
  geom_histogram(fill = "dodgerblue3", colour = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  scale_x_continuous(labels = scales::comma, breaks = c(30,31,32,33,34,35,36,37,38,39,40)) +
  theme_bw()

# unfiltered read histograms: coverage, no.reads, base quality, mapping quality
mean_coverage_hist <- ggplot(data = all_coverages, aes(meandepth)) +
  geom_histogram(fill = "orange", colour = "black") +
  theme_bw()
mean_reads_hist <- ggplot(data = all_coverages, aes(numreads)) +
  geom_histogram(fill = "orange", colour = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  scale_x_continuous(labels = scales::comma) +
  theme_bw()
mean_baseq_hist <-  ggplot(data = all_coverages, aes(meanbaseq)) +
  geom_histogram(fill = "orange", colour = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  scale_x_continuous(labels = scales::comma) +
  theme_bw()
mean_mapq_hist <- ggplot(data = all_coverages, aes(meanmapq)) +
  geom_histogram(fill = "orange", colour = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  scale_x_continuous(labels = scales::comma, breaks = c(16,18,20,22,24,26,28,30,32)) +
  theme_bw()


 # filtered read plots: coverage (x axis: SRR), no.reads, base quality, mapping quality
sample_coverage_plot_qfilt <- ggplot(data = all_coverages_qfilt, aes(SRRs, calculated_mean_depth, calculated_sd)) +
  geom_col(colour = "black", fill = "dodgerblue3") +
  geom_errorbar(aes(ymin=calculated_mean_depth-calculated_sd_depth, ymax=calculated_mean_depth+calculated_sd_depth), width=0) +
  scale_y_continuous(trans='log2', expand = expansion(mult = c(0, .1))) +
  theme(axis.text.x = element_text(angle = 45, vjust=1.05, hjust = 1.0),
        axis.ticks.x = element_line(),
        panel.background = element_rect(fill = "white"))
sample_reads_plot_qfilt <- ggplot(data = all_coverages_qfilt, aes(SRRs, numreads)) +
  geom_col(colour = "black", fill = "dodgerblue3") +
  scale_y_continuous(expand = expansion(mult = c(0, .1)), labels = scales::comma) +
  theme(axis.text.x = element_text(angle = 45, vjust=1.05, hjust = 1.0),
        axis.ticks.x = element_line(),
        panel.background = element_rect(fill = "white")) 
sample_baseq_plot_qfilt <-  ggplot(data = all_coverages_qfilt, aes(SRRs, meanbaseq)) +
  geom_col(colour = "black", fill = "dodgerblue3") +
  scale_y_continuous(expand = expansion(mult = c(0, .1)), labels = scales::comma) +
  theme(axis.text.x = element_text(angle = 45, vjust=1.05, hjust = 1.0),
        axis.ticks.x = element_line(),
        panel.background = element_rect(fill = "white")) 
sample_mapq_plot_qfilt <- ggplot(data = all_coverages_qfilt, aes(SRRs, meanmapq)) +
  geom_col(colour = "black", fill = "dodgerblue3") +
  scale_y_continuous(expand = expansion(mult = c(0, .1)), labels = scales::comma) +
  theme(axis.text.x = element_text(angle = 90, vjust=1.05, hjust = 1.0),
        axis.ticks.x = element_line(),
        panel.background = element_rect(fill = "white")) 

# unfiltered read plots: coverage (x axis: SRR), no.reads, base quality, mapping quality
sample_coverage_plot <- ggplot(data = all_coverages, aes(SRRs, calculated_mean_depth, calculated_sd_depth)) +
  geom_col(colour = "black", fill = "orange") +
  geom_errorbar(aes(ymin=calculated_mean_depth-calculated_sd_depth, ymax=calculated_mean_depth+calculated_sd_depth), width=0) +
  scale_y_continuous(trans='log2', expand = expansion(mult = c(0, .1))) +
  theme(axis.text.x = element_text(angle = 45, vjust=1.05, hjust = 1.0, size = 10),
        axis.ticks.x = element_line(),
        panel.background = element_rect(fill = "white"))
sample_reads_plot <- ggplot(data = all_coverages, aes(SRRs, numreads)) +
  geom_col(colour = "black", fill = "orange") +
  scale_y_continuous(expand = expansion(mult = c(0, .1)), labels = scales::comma) +
  theme(axis.text.x = element_text(angle = 45, vjust=1.05, hjust = 1.0),
        axis.ticks.x = element_line(),
        panel.background = element_rect(fill = "white")) 
sample_baseq_plot <-  ggplot(data = all_coverages, aes(SRRs, meanbaseq)) +
  geom_col(colour = "black", fill = "orange") +
  scale_y_continuous(expand = expansion(mult = c(0, .1)), labels = scales::comma) +
  theme(axis.text.x = element_text(angle = 45, vjust=1.05, hjust = 1.0),
        axis.ticks.x = element_line(),
        panel.background = element_rect(fill = "white")) 
sample_mapq_plot <- ggplot(data = all_coverages, aes(SRRs, meanmapq)) +
  geom_col(aes(fill = factor(read_lengths))) +
  scale_fill_manual(values = c("150"="skyblue3","76"="plum3"), name = "Read length", labels = c("2x75bp", "2x35bp")) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)), labels = scales::comma) +
  theme(axis.text.x = element_text(angle = 90, vjust=1.05, hjust = 1.0, size = 10),
        text=element_text(size=13), axis.ticks.x = element_line(),
        panel.background = element_rect(fill = "white"), legend.position=c(.8, .2)) +
  labs(x = "Clone", y = "Mean read mapping quality (unfiltered)")
ggsave(file="results/sample_mapq_plot.png", plot=sample_mapq_plot, width = 12, height = 5, units = "in")


sample_coverage_boxplot_qfilt <- ggplot(data = depths_qfilt_bplot_data, aes(SRRs, read_depths)) +
  geom_boxplot() +
  scale_y_continuous(trans='log2', expand = expansion(mult = c(0, .1))) +
  geom_hline(yintercept=60, size = 1, colour = "red", linetype = "dotted") +
  theme(axis.text.x = element_text(angle = 90, vjust=1.05, hjust = 1.0, size = 10),
        axis.ticks.x = element_line(), text=element_text(size=13),
        panel.background = element_rect(fill = "white")) +
  expand_limits(y = 0) +
  labs(x = "Clone", y="Read depth")


# Unfiltered:
mean_coverage_hist
mean_baseq_hist
mean_mapq_hist
mean_reads_hist

sample_coverage_plot
sample_baseq_plot
sample_mapq_plot
sample_reads_plot

# Base and mapping quality filtered
mean_coverage_hist_qfilt
mean_baseq_hist_qfilt
mean_mapq_hist_qfilt
mean_reads_hist_qfilt

sample_coverage_plot_qfilt
sample_baseq_plot_qfilt
sample_mapq_plot_qfilt
sample_reads_plot_qfilt

fig2_plots <- list(pos_coverage_all_plot,sample_coverage_boxplot_qfilt)
fig2_plots <- ggarrange(plots = fig2_plots, nrow = 2, align = "v")
ggsave(file="results/fig2_plots.png", plot=fig2_plots, width = 12, height = 8, units = "in")

}  # end if (post_plots == TRUE)


  #########################  VARIANT DATAFRAMES  ############################

SRR_table_list <- list()  # All 
SRR_table_list_PASS <- list()  # Filtered
SRR_table_list_HET_OR_LOWLVL <- list()  # Filtered and only heteroplasmic/low level
SRR_table_list_HET_OR_LOWLVL_nofilt <- list()
SRR_table_list_HET_OR_LOWLVL_validated <- list()
SRR_table_list_HET_OR_LOWLVL_potautocor_validated <- list()

#threshold <- 0.05
# Load files into list of data.frames
for(i in SRR_names){
  if (use_pileups==TRUE){
    SRR_table_list[[i]] <- bcf_SRR_table_list[[i]]
  }
  else {
    filepath <- file.path(vcfdir,paste0(i,"_annotated.txt"))
    SRR_table_list[[i]] <- read.table(filepath, sep = "\t", header = T, stringsAsFactors = T)
    }
  # remove positions (deletion eg. at 3107)
  SRR_table_list[[i]] <- SRR_table_list[[i]][!(SRR_table_list[[i]]$Pos==3107),]
  # subset for variants which passed filter
  SRR_table_list_PASS[[i]] <- subset(SRR_table_list[[i]], Filter == "PASS")
  # subset for "HET_OR_LOWLVL" variants (Heteroplasmic or low-level variant)
  SRR_table_list_HET_OR_LOWLVL[[i]] <- subset(SRR_table_list_PASS[[i]], Type == 2)
  # no strand bias filter to see if variants are filtered differently in different samples across lineage paths
  SRR_table_list_HET_OR_LOWLVL_nofilt[[i]] <- subset(SRR_table_list[[i]], Type ==2)
  # Lineage_validated created below 
}


  #################### Variant calling stats ########################

bulk_variant_pos80 <- data.frame(SRR_table_list_HET_OR_LOWLVL$SRR7245880$Pos)#[SRR_table_list$SRR7245880$Type==2])
bulk_variant_pos81 <- data.frame(SRR_table_list_HET_OR_LOWLVL$SRR7245881$Pos)#[SRR_table_list$SRR7245881$Type==2])
bulk_variant_pos <- merge(bulk_variant_pos80, bulk_variant_pos81, by=1, all=T)
colnames(bulk_variant_pos) <- "Bulk_Variants"
write.csv(bulk_variant_pos, file = "results/bulk_variant_positions.csv", quote = F)

### make data frame of positions of all our variants ###

all_variants <-  data.frame(matrix(ncol = 1))
colnames(all_variants) <- "Pos"
all_variants <- list()
for (SRR in SRR_names) {
  SRR_pos_level <- data.frame(SRR_table_list_HET_OR_LOWLVL_nofilt[[SRR]]$Pos, SRR_table_list_HET_OR_LOWLVL_nofilt[[SRR]]$VariantLevel)
  colnames(SRR_pos_level) <- c("Pos", paste0(SRR,"_variant_lvl"))
  all_variants <- c(all_variants, SRR_pos_level$Pos)
  all_variants <- unique(all_variants)
  #all_variants <- rbind(all_variants, SRR_pos_level$Pos)
  #all_variants <- all_variants[!duplicated(all_variants$Pos), ]
  #all_variants <- merge(all_variants$Pos, SRR_pos_level$Pos, by = "Pos", all = TRUE)
}
all_variants <- all_variants[!duplicated(all_variants)]


   ####  all_variants_in_path  ####
# merge to list all variants in lineage path
# replace position with variant level.
# if min variant level is < eg. 0.95 change filter to FIXED_IN_BULK.

# Produce table of the variant levels in each sample of all variants found in the lineage path
for (p in paths){
  if (p[[1]] == "#"){
    print("skipping comment line...")
    next
  }
  all_variants_in_lineage <- data.frame(matrix(ncol = 1))
  colnames(all_variants_in_lineage) <- "Pos"
  n=0

  for (SRR in p){
    # Skip Lineage path name (1st in character vector of paths[[p]] )
    n=n+1
    print(SRR)
    #print(typeof(SRR))
    if (n==1){
      next
    }
    # stop where the variants of interest positions are listed on the line (for 
    # mutation load plots below)
    if (SRR == "VARIANTS_OF_INTEREST") {
      print("Reached VARIANTS_OF_INTEREST for this lineage")
      break
    }
    SRR_pos_level <- data.frame(SRR_table_list[[SRR]]$Pos, SRR_table_list[[SRR]]$VariantLevel)
    print(colnames(SRR_pos_level))
    colnames(SRR_pos_level) <- c("Pos", paste0(SRR,"_variant_lvl"))
    all_variants_in_lineage <- merge(all_variants_in_lineage,SRR_pos_level, by = "Pos", all = TRUE)
  }
file_string <- paste0("results/",p[[1]],"_all_variants.csv")
write.csv(all_variants_in_lineage,file = file_string, quote = F)
print("table of bulk variants in lineage path saved in 'results/'")
}


     ##################  LINEAGE VALIDATION  ###########################

# Repeat to produce table for only unfiltered heteroplasmic or low level variants (ie. incl strand bias), and another table for lineage validated.
## CHECK ## validation_groups.txt:
  # - does not include LUDWIG
  # - does not include "mix" samples
  # - (for lineage mutation load plots:) all lineages start with the two parent "bulks" SRR7245880 and 81.
validation_paths <- list()
#validation_paths <- as.list(strsplit(readLines("validation_groups.txt"), " "))
# Use lineage paths: better validation through each path, rather than whole lineage group.
validation_paths <- as.list(strsplit(readLines("lineage_paths.txt"), " "))

all_variants_in_lineages <- list()
validated_per_lineage <- list()
all_potential_autocor_poses <- c()
all_validated_poses <- c()
all_validated_pos_base <- data.frame(matrix(ncol = 2))
colnames(all_validated_pos_base) <- c("Pos", "Base")
#all_lineages_validated <- data.frame(matrix(ncol = 1)) 
#colnames(all_lineages_validated) <- "Pos"
#all_potential_autocor_poses <- data.frame(matrix(ncol = 1))
#colnames(all_potential_autocor_poses) <- "Pos"
for (p in validation_paths){
  if (p[[1]] == "#"){
    print("skipping comment line...")
    next
  }
  if (str_detect(p[[1]], 'LUDWIG')){
    print(paste0("skipping LUDWIG line: ", p[[1]]))
  }
  
  n_alleles_in_SRR <- 0
  for (base in c("A","C","G","T")){

  n=0
  print("stage 1")

  all_variants_in_lineage_HET_OR_LOWLVL_nofilt <- data.frame(matrix(ncol = 1))
  colnames(all_variants_in_lineage_HET_OR_LOWLVL_nofilt) <- "Pos"
  
  for (SRR in p){
    # Skip Lineage path name (1st in character vector of paths[[p]] )
    n=n+1
    print(SRR)
    if (n==1){
      print(paste0("Lineage name: ", SRR))
      next
    }
    # stop where the variants of interest positions are listed on the line (for 
    # mutation load plots below)
    if (SRR == "VARIANTS_OF_INTEREST") {
      print("Reached VARIANTS_OF_INTEREST for this lineage")
      break
    }
    print("stage 2")
    print(base)
    if (any(SRR_table_list_HET_OR_LOWLVL_nofilt[[SRR]]$Variant == base)){
      SRR_table_1base <- SRR_table_list_HET_OR_LOWLVL_nofilt[[SRR]] %>% filter(Variant == base)
      SRR_pos_level <- data.frame(SRR_table_1base$Pos, SRR_table_1base$VariantLevel)
      #SRR_pos_level <- data.frame(Pos=unique(c(SRR_table_list_HET_OR_LOWLVL_nofilt[[SRR]]$Pos)), VariantLevel=rep(0,16570))
      #for (pos in SRR_table_1base$Pos) {
      #  SRR_pos_level$VariantLevel[SRR_pos_level$Pos == pos] <- SRR_table_1base$VariantLevel[SRR_table_1base$Pos == pos]
      #  #[SRR_table_list_HET_OR_LOWLVL_nofilt[[SRR]]$Variant == base & SRR_table_list_HET_OR_LOWLVL_nofilt[[SRR]]$Pos == pos]
      #}
      print(colnames(SRR_pos_level))
      colnames(SRR_pos_level) <- c("Pos", paste0(SRR,"_variant_lvl"))
      all_variants_in_lineage_HET_OR_LOWLVL_nofilt <- merge(all_variants_in_lineage_HET_OR_LOWLVL_nofilt, SRR_pos_level, by = "Pos", all = T)
    }
  }

  ###  LINEAGE VALIDATION  ###
  
  # only keep variants which have at least one allele with an allele frequency > 0.01, and are present at least twice in the lineage
  lineage_validated_1base <- all_variants_in_lineage_HET_OR_LOWLVL_nofilt %>% filter_at(-1, any_vars(.>0.01))
  lineage_validated_1base <- lineage_validated_1base[rowSums(!is.na(lineage_validated_1base[,-1]))>=2,]
  
  # Merge so the lineage validated of all four bases are in one multiallelic table
  if (dim(lineage_validated_1base)[1] != 0) {
    lineage_validated_1base$Base <- base
    n_alleles_in_SRR <- n_alleles_in_SRR + 1
    print(paste0("No. Duplicate rows: ", duplicated(lineage_validated_1base)))
    if (n_alleles_in_SRR <= 1){  # difficult to create empty df for each lineage at the beginning without knowing length of lineage
      lineage_validated <- lineage_validated_1base
    }
    else{
      lineage_validated <- rbind(lineage_validated, lineage_validated_1base)
    }
  } 
  } # end ACGT loop
  
  #should_be_no_dups_nrow <- nrow()
  print("stage 4")
  file_string <- paste0("results/",p[[1]],"_HET_nofilt_variants.csv")
  write.csv(all_variants_in_lineage_HET_OR_LOWLVL_nofilt,file = file_string, quote = F)
  print("table of variants in lineage path saved in 'results/'")
  all_variants_in_lineages[[ p[[1]] ]] <- all_variants_in_lineage_HET_OR_LOWLVL_nofilt
  
  print("stage 6")
  
  file_string <- paste0("results/",p[[1]],"_lineage_validated_mutations.csv")
  write.csv(lineage_validated,file = file_string, quote = F)
  print("table of bulk variants in lineage path saved in 'results/'")
  
  # Get mutations which could be validated though autocorrelation (occur >=2 x in lineage, not necessary to also have one below 0.01 heteroplasmy)
  potential_autocorrelation <- all_variants_in_lineage_HET_OR_LOWLVL_nofilt[rowSums(!is.na(all_variants_in_lineage_HET_OR_LOWLVL_nofilt[,-1]))>=2,]
  file_string <- paste0("results/",p[[1]],"_pot_valid_by_autocorrelation.csv")
  write.csv(potential_autocorrelation, file = file_string, quote = F)
  print("stage 7")
  
  # Add new lineage validated mutations to table with all lineage mutations, same for potentially validated with autocorrelation
  #all_lineages_validated <- merge(all_lineages_validated, lineage_validated, by.x = "Pos", by.y = "Pos", all = T) 
  all_validated_poses <- unique(c(all_validated_poses, c(lineage_validated$Pos)))
  all_validated_pos_base <- unique(rbind(all_validated_pos_base, lineage_validated[,c("Pos", "Base")]))
  #file_string <- paste0("results/",p[[1]],"_lineage_validated.csv")
  #write.csv(all_, file = file_string, quote = F)
  print("stage 7.5")
  all_potential_autocor_poses <- unique(c(all_potential_autocor_poses, c(potential_autocorrelation$Pos)))
  #all_potential_autocor_poses <- merge(all_potential_autocor_poses, potential_autocorrelation, by = "Pos", all = T)
  #all_potential_autocor_poses <- all_potential_autocor_poses[,1:3]
  print("stage 8")
}

all_validated_pos_base <- all_validated_pos_base %>% drop_na()

# write table of all lineage validated variants from any lineage
#all_lineages_validated <- all_lineages_validated[!duplicated(all_lineages_validated$Pos), ]  # potential removal of variant levels on repeated SRRs?
file_string <- paste0("results/all_variants_lineage_validated.csv")
write.csv(all_validated_poses,file = file_string, quote = F) 
#all_lineages_validated <- subset(all_variants, Pos %in% all_lineages_validated)
#all_lineages_validated <- all_lineages_validated[!duplicated(all_lineages_validated$Pos), ]

# write table of variants which could potentially be validated with autocorrelation (present in more than one sample, not necessarily any with an AF>0.01)
#all_potential_autocor_poses <- subset(all_variants, Pos %in% potential_autocor_poses)
#all_potential_autocor_poses <- all_potential_autocor_poses[!duplicated(all_potential_autocor_poses$Pos), ]
file_string <- paste0("results/all_variants_pot_valid_by_autocorrelation.csv")
write.csv(all_potential_autocor_poses,file = file_string, quote = F)

# Create SRR_table equivalents for validated and potentially autocorrelated.
all_lineages_validated_poses_df <- data.frame(all_validated_poses)
colnames(all_lineages_validated_poses_df) <- "Pos"
all_potential_autocor_poses_df <- data.frame(all_potential_autocor_poses)
colnames(all_potential_autocor_poses_df) <- "Pos"
# lineage validated, 
for(i in SRR_names) {
  SRR_table_list_HET_OR_LOWLVL_validated[[i]] <- merge(all_lineages_validated_poses_df, SRR_table_list_HET_OR_LOWLVL_nofilt[[i]], by.x = "Pos", by.y = "Pos")
  print(nrow(SRR_table_list_HET_OR_LOWLVL_validated[[i]]))
  # Don't remove duplicates in $Pos - due to multiallelic pileups
  #SRR_table_list_HET_OR_LOWLVL_validated[[i]] <- SRR_table_list_HET_OR_LOWLVL_validated[[i]][!duplicated(SRR_table_list_HET_OR_LOWLVL_validated[[i]]$Pos), ]
  print(nrow(SRR_table_list_HET_OR_LOWLVL_validated[[i]]))
  #SRR_table_list_HET_OR_LOWLVL_validated[[i]] <- SRR_table_list_HET_OR_LOWLVL_validated[[i]][,1:18] %>% filter(drop_na())
  #SRR_table_list_HET_OR_LOWLVL_validated[[i]] %>% subset(SRR_table_list_HET_OR_LOWLVL_validated[[i]], Filter !="NA")
  SRR_table_list_HET_OR_LOWLVL_potautocor_validated[[i]] <- merge(all_potential_autocor_poses_df, SRR_table_list_HET_OR_LOWLVL_nofilt[[i]], by.x = "Pos", by.y = "Pos")
  #print(nrow(SRR_table_list_HET_OR_LOWLVL_potautocor_validated[[i]]))
}

#SRR_table_list_HET_OR_LOWLVL_validated[[i]][!is.na(SRR_table_list_HET_OR_LOWLVL_validated[[i]]$ID),]


  # Get total no validated mutations
n=0
for (i in SRR_table_list_HET_OR_LOWLVL_validated){
  n=n+nrow(i)
}
print(n)

  # Get total no nofilt mutations
n=0
for (i in SRR_table_list_HET_OR_LOWLVL_nofilt){
  n=n+nrow(i)
}
print(n)

  # no het or lowlvl pass (not strand bias)
n=0
for (i in SRR_table_list_HET_OR_LOWLVL){
  n=n+nrow(i)
}
print(n)

  # No. potentially validated by autocorrelation
n=0
for (i in SRR_table_list_HET_OR_LOWLVL_potautocor_validated){
  n=n+nrow(i)
}
print(n)

  ####################  Known sequencing errors  ####################
# for NextSeq 500:

## Context:
# a base of the same type as the one preceding the error, eg. CG -> CC
# EXCEPTIONS for NextSeq 500: TA -> TT (for previous base, A is overrepresented instead), AC -> AA, AT -> AA,
## post-homopolymer error
# eg. AAAT
#        ^ if eg. C-> T: homopolymer of 3 preceding the error  

## Prev 3 bases:
# GGT
# CGT
# AGT
# CGA
# CCA
# GCT
# 

# string detect
## get index of base
# if index in linval$pos, add to list 





   ###################### Variant stats ########################

#variant_stats <- data.frame(matrix(nrow = length(SRR_names), ncol = 10))
#colnames(variant_stats) <- c("SRR", "No.Variants", "No.het", "No.hom", "No.lineage.validated","No.transition", "No.transversion", "Ts/Tv", "No.missense", "Strand bias")
#variant_stats$SRR <- SRR_names
myfunc <- function(x){
  string = "transversion"
  ret <- filter(x, x$Substition == string) %>% nrow
  return(ret)
}

variant_summaries <- list()
variant_summaries_dfnames <- c("HET_OR_LOWLVL_nofilt", "HET_OR_LOWLVL", "HET_OR_LOWLVL_validated")
#n=0
#table_list <- list(SRR_table_list_HET_OR_LOWLVL_nofilt)#, SRR_table_list_HET_OR_LOWLVL, SRR_table_list_HET_OR_LOWLVL_validated)
#for (table in table_list){
  #n=n+1
  variant_stats_by_SRR <- data.frame(SRR_names)
  #No.transversions = sapply(table, myfunc)
  #print(No.transversions)
  #No.Strand.Bias.Per.SRR <- 
  #Mean.Coverage <- 
  No.Variants.Per.SRR <- sapply(SRR_table_list_HET_OR_LOWLVL_nofilt, nrow)
  No.Variants <- sum(sapply(SRR_table_list_HET_OR_LOWLVL_nofilt, nrow))
  No.Positions <- length(all_variants)
  No.Validated.Variants.Per.SRR <- sapply(SRR_table_list_HET_OR_LOWLVL_validated, nrow)
  No.Validated.Variants <- sum(sapply(SRR_table_list_HET_OR_LOWLVL_validated, nrow))
  No.Validated.Positions <- length(all_validated_poses)
  No.Potential.AutoCor.Valid.Variants.per.SRR <- sapply(SRR_table_list_HET_OR_LOWLVL_potautocor_validated, nrow)
  No.Potential.AutoCor.Valid.Variants <- sum(sapply(SRR_table_list_HET_OR_LOWLVL_potautocor_validated, nrow))
  No.Potential.AutoCor.Valid.Postions <- length(all_potential_autocor_poses)
  No.Variants.Per.SRR.PASS.FILTERS <- sapply(SRR_table_list_HET_OR_LOWLVL, nrow)
  No.Variants.PASS <- sum(sapply(SRR_table_list_HET_OR_LOWLVL, nrow))
  No.Variants.Per.SRR.nofilt <- sapply(SRR_table_list_HET_OR_LOWLVL_nofilt, nrow)
  No.Variants.nofilt <- sum(sapply(SRR_table_list_HET_OR_LOWLVL_nofilt, nrow))
  No.Variants.Per.SRR.STRAND.BIAS <- No.Variants.Per.SRR.nofilt - No.Variants.Per.SRR.PASS.FILTERS  # implement strand bias check and filter
  No.Variants.STRAND.BIAS <- No.Variants.nofilt - No.Variants.PASS
  
  variant_summaries <- data.frame(No.Variants, No.Positions, No.Validated.Variants, No.Validated.Positions, No.Potential.AutoCor.Valid.Variants, No.Potential.AutoCor.Valid.Postions, No.Variants.STRAND.BIAS)
  print(variant_summaries)
  file_string <- paste0("results/variant_summaries.csv")
  write.csv(variant_summaries, file = file_string, quote = F)
  #table[sapply(nrow)]
#  for (SRR_name in SRR_names){
    #No.Variants <- nrow(table[[SRR_name]])
    #print(paste("No. Variants: ", No.Variants))
    #print(paste("No. mutations (Sum No.Variants): ", sum(No.Variants)))
    #No.transversions <- c(No.transversions, nrow(table[[SRR_name]][which(Substitution == "transversion"), ]))
    #No.transversions <- c(No.transversions, nrow(table[[SRR_name]][which(Substitution =="transition"), ]), No.transversions)
    #No.Strand.Bias <- c(No.Strand.Bias, nrow(table[[SRR_name]][which(table[[SRR_name]]$Filter == "STRAND_BIAS"), ]))
    #variant_stats$`Ts/Tv`[[i]] <- variant_stats$No.transition[[i]]/variant_stats$No.transversion[[i]]
  #}
  #variant_stats_by_SRR <- data.frame(SRR_names)
  #variant_stats_by_SRR$No.Variants <- No.Variants #st(as.data.frame(No.Variants))
  #variant_stats_by_SRR$No.transitions <- No.transitions
  #variant_stats_by_SRR$No.transversions <- t(as.data.frame(No.transversions))
  #variant_stats_by_SRR$No.Strand.Bias <- t(as.data.frame(No.Strand.Bias))
  #print(variant_stats_by_SRR)
#  variant_summaries <- c(variant_summaries, variant_stats_by_SRR)
  #variant_stats[[table]]
#}
#colnames(variant_stats) <- c("statistic", variant_summaries_dfnames)
#variant_stats
#variant_stats$no.mutations <- nrow


#### Lineage stats #####


   ###############  Comparison with Ludwig's variants  #################
if (Ludwig_comparison == TRUE){
# Ludwig's research aligned to the Hg19 mitochondrial genome, which has indels compared to rCRS. Positions off by 0,-2,-1,-2 in different parts of the chromosome. See converted positions below.
# read in Ludwigs variants, and variant level in each sample
Ludwig_variants <- read.csv("data/LUDWIG_TF1_clones_ATAC_alleleFrequencies.csv", header = T)
colnames(Ludwig_variants) <- c("Ludwig_variant_positions", SRR_names)
Ludwig_variants$tobecombined_Pos <- Ludwig_variants$Ludwig_variant_positions

# Combine with Ludwigs data (all unconverted variant positions: Hg19-rCRS)
all_variants_df <- data_frame(Pos=all_variants)
all_variants_and_Ludwigs <- merge(all_variants_df, Ludwig_variants, by.x = "X1", by.y = "tobecombined_Pos", all = T)
all_pos_and_Ludwigs <- data.frame(all_variants_and_Ludwigs$Pos, all_variants_and_Ludwigs$OurPos, all_variants_and_Ludwigs$Ludwig_variant_positions)

# Repeat for differently filtered data
all_variants_HET_OR_LOWLVL <-  data.frame(matrix(ncol = 1))
colnames(all_variants_HET_OR_LOWLVL) <- "Pos"

for (SRR in SRR_names) {
  SRR_pos_level <- data.frame(SRR_table_list_HET_OR_LOWLVL[[SRR]]$Pos, SRR_table_list_HET_OR_LOWLVL[[SRR]]$VariantLevel)
  colnames(SRR_pos_level) <- c("Pos", paste0(SRR,"_variant_lvl"))
  all_variants_HET_OR_LOWLVL <- merge(all_variants_HET_OR_LOWLVL, SRR_pos_level, by = "Pos", all = T)
}
all_variants_HET_OR_LOWLVL <- all_variants_HET_OR_LOWLVL[!duplicated(all_variants_HET_OR_LOWLVL$Pos), ]
all_variants_HET_OR_LOWLVL$OurPos <- all_variants_HET_OR_LOWLVL$Pos
all_variants_HET_OR_LOWLVL_and_Ludwigs <- merge(all_variants_HET_OR_LOWLVL, Ludwig_variants, by.x = "Pos", by.y = "tobecombined_Pos", all = T)
all_pos_HET_OR_LOWLVL_and_Ludwigs <- data.frame(all_variants_HET_OR_LOWLVL_and_Ludwigs$Pos, all_variants_HET_OR_LOWLVL_and_Ludwigs$OurPos, all_variants_HET_OR_LOWLVL_and_Ludwigs$Ludwig_variant_positions)


all_variants_HET_OR_LOWLVL_nofilt <-  data.frame(matrix(ncol = 1))
colnames(all_variants_HET_OR_LOWLVL_nofilt) <- "Pos"

for (SRR in SRR_names) {
  SRR_pos_level <- data.frame(SRR_table_list_HET_OR_LOWLVL_nofilt[[SRR]]$Pos, SRR_table_list_HET_OR_LOWLVL_nofilt[[SRR]]$VariantLevel)
  colnames(SRR_pos_level) <- c("Pos", paste0(SRR,"_variant_lvl"))
  all_variants_HET_OR_LOWLVL_nofilt <- merge(all_variants_HET_OR_LOWLVL_nofilt, SRR_pos_level, by = "Pos", all = T)
}
all_variants_HET_OR_LOWLVL_nofilt$OurPos <- all_variants_HET_OR_LOWLVL_nofilt$Pos
all_variants_HET_OR_LOWLVL_nofilt_and_Ludwigs <- merge(all_variants_HET_OR_LOWLVL_nofilt, Ludwig_variants, by.x = "Pos", by.y = "tobecombined_Pos", all = T)
all_pos_HET_OR_LOWLVL_nofilt_and_Ludwigs <- data.frame(all_variants_HET_OR_LOWLVL_nofilt_and_Ludwigs$Pos, all_variants_HET_OR_LOWLVL_nofilt_and_Ludwigs$OurPos, all_variants_HET_OR_LOWLVL_nofilt_and_Ludwigs$Ludwig_variant_positions)

# List of Ludwig converted variant positions
Ludwig_pos <-  c("182","309","824","849","1412","1495","1797","1972","2110","2818","3174","3911","4038","4114","4215","4447","4513","5008","5564","5863","6076","6963","7075","7790","8003","8207","8922","10372","11185","11404","11712","12062","12254","12790","12839","13289","13413","13709","14437","15089","15489","15641","15798","16252")
rCRS_Ludwig_pos <- c("182","309","822","847","1410","1493","1795","1970","2108","2816","3173","3910","4037","4413","4214","4446","4512","5007","5563","5862","6075","6962","7074","7789","8002","8206","8921","10371","11184","11403","11711","12061","12253","12789","12838","13288","13412","13708","14436","15088","15488","15640","15797","16250")
Ludwig_positions <- data.frame(Ludwig_pos,rCRS_Ludwig_pos)

# add column of converted positions to Ludwig_variants
Ludwig_variants <- merge(Ludwig_positions, Ludwig_variants, by.x = "Ludwig_pos", by.y = "Ludwig_variant_positions")

# Extract rCRS Ludwig variants in our data
HET_OR_LOWLVL_nofilt_Ludwig_variants <-  data.frame(rCRS_Ludwig_pos)
HET_OR_LOWLVL_nofilt_Ludwig_variants$rCRS_Ludwig_pos <- rCRS_Ludwig_pos
our_Ludwig_variants_nofilt <-  data.frame(rCRS_Ludwig_pos)
our_Ludwig_variants_nofilt$rCRS_Ludwig_pos <- rCRS_Ludwig_pos
for (SRR in SRR_names){
  SRR_pos_level <- data.frame(SRR_table_list_HET_OR_LOWLVL_nofilt[[SRR]]$Pos,SRR_table_list_HET_OR_LOWLVL_nofilt[[SRR]]$VariantLevel)
  colnames(SRR_pos_level) <- c("Pos", SRR)
  HET_OR_LOWLVL_nofilt_Ludwig_variants <- merge(HET_OR_LOWLVL_nofilt_Ludwig_variants, SRR_pos_level, by.x="rCRS_Ludwig_pos", by.y = "Pos", all.x=T)
  
  SRR_pos_level <- data.frame(SRR_table_list[[SRR]]$Pos,SRR_table_list[[SRR]]$VariantLevel)
  colnames(SRR_pos_level) <- c("Pos", SRR)
  our_Ludwig_variants_nofilt <- merge(our_Ludwig_variants_nofilt, SRR_pos_level, by.x="rCRS_Ludwig_pos", by.y = "Pos", all.x=T)
}


  ## Overall correlation ##
# prepare data
melted_our_Ludwig_variants <- reshape2::melt(HET_OR_LOWLVL_nofilt_Ludwig_variants)
Ludwig_variants <- Ludwig_variants[,-72]
melted_Ludwig_variants <- reshape2::melt(Ludwig_variants[,-1])
colnames(melted_our_Ludwig_variants) <- c("rCRS_Ludwig_pos", "SRR_names", "our_variant_level")
colnames(melted_Ludwig_variants) <- c("rCRS_Ludwig_pos", "SRR_names", "Ludwigs_variant_level")
melted_our_Ludwig_variants[is.na(melted_our_Ludwig_variants)] = 0
melted_correlation_data <- merge(melted_Ludwig_variants, melted_our_Ludwig_variants, by = c("rCRS_Ludwig_pos", "SRR_names"))
melted_correlation_data$Ludwigs_variant_level <- sqrt(melted_correlation_data$Ludwigs_variant_level)
melted_correlation_data$our_variant_level <- sqrt(melted_correlation_data$our_variant_level)

corr_plot_all <- ggplot(melted_correlation_data, aes(Ludwigs_variant_level, our_variant_level)) +
  geom_point(aes(colour=factor(rCRS_Ludwig_pos))) +
  scale_x_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) +
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) +
  #geom_smooth(method = "lm", se = TRUE, color = 'black', aes(alpha = 0.5)) +
  theme(text = element_text(size = 13), legend.key.size = unit(0.2, "cm"), legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10), legend.position = "none") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=13), legend.position = "none") +
  labs(x = "sqrt(Allele Frequency): Ludwig et al, 2019", y = "sqrt(Allele Frequency)", colour = "Variants")+
  expand_limits(x=1,y=1) +
  geom_abline(intercept = 0, slope = 1, alpha = 0.4)
ggsave(file="results/correlation_with_Luwig_by_pos.png", plot=corr_plot_all, width = 8, height = 4, units = "in")

#shapiro.test(melted_correlation_data$Ludwigs_variant_level)

cor.test(melted_correlation_data$Ludwigs_variant_level, melted_correlation_data$our_variant_level, method = 'pearson')
cor.test(melted_correlation_data$Ludwigs_variant_level, melted_correlation_data$our_variant_level, method = 'spearman')


  ## Spearman's rank correlation per position ##

Ludwig_spearman_by_pos <- data.frame(rCRS_Ludwig_pos)
spear_res_tmp <- data.frame(rCRS_Ludwig_pos)
estimates <- list()
pvalues <- list()
for (i in rCRS_Ludwig_pos) {
  print(i)
  n=n+1
  spear_data <- melted_correlation_data[melted_correlation_data$rCRS_Ludwig_pos==i,]
  spear_res_tmp<- cor.test(spear_data$Ludwigs_variant_level, spear_data$our_variant_level, method = 'spearman')
  estimates <- c(estimates, spear_res_tmp$estimate)
  pvalues <- c(pvalues,spear_res_tmp$p.value)
}
Ludwig_spearman_by_pos$rho <- estimates
Ludwig_spearman_by_pos$p.values <- pvalues

Ludwig_spearman_by_pos <- apply(Ludwig_spearman_by_pos,2,as.character)
file_string <- paste0("results/Correlation_spearman_perPos.csv")
write.csv(Ludwig_spearman_by_pos,file = file_string, quote = F)


  ## Pearson's correlation per position ##

Ludwig_pearson_by_pos <- data.frame(rCRS_Ludwig_pos)
pearson_res_tmp <- data.frame(rCRS_Ludwig_pos)
estimates <- list()
pvalues <- list()
for (i in rCRS_Ludwig_pos) {
  print(i)
  n=n+1
  pearson_data <- melted_correlation_data[melted_correlation_data$rCRS_Ludwig_pos==i,]
  pearson_res_tmp<- cor.test(pearson_data$Ludwigs_variant_level, pearson_data$our_variant_level, method = 'pearson')
  estimates <- c(estimates, pearson_res_tmp$estimate)
  pvalues <- c(pvalues,pearson_res_tmp$p.value)
}
Ludwig_pearson_by_pos$rho <- estimates
Ludwig_pearson_by_pos$p.values <- pvalues

Ludwig_pearson_by_pos <- apply(Ludwig_pearson_by_pos,2,as.character)
file_string <- paste0("results/Correlation_pearson_perPos.csv")
write.csv(Ludwig_pearson_by_pos,file = file_string, quote = F)


 ############## Heatmap of Ludwig's variant positions ###################

#HET_OR_LOWLVL_nofilt_Ludwig_variants[is.na(HET_OR_LOWLVL_nofilt_Ludwig_variants)] <- as.numeric(0)
#htmp_HET_OR_LOWLVL_nofilt_Ludwig_variants <- HET_OR_LOWLVL_nofilt_Ludwig_variants[,-1]
#rownames(htmp_HET_OR_LOWLVL_nofilt_Ludwig_variants) <- HET_OR_LOWLVL_nofilt_Ludwig_variants$rCRS_Ludwig_pos
#colnames(htmp_HET_OR_LOWLVL_nofilt_Ludwig_variants) <- SRR_names
#htmp_HET_OR_LOWLVL_nofilt_Ludwig_variants <- sqrt(htmp_HET_OR_LOWLVL_nofilt_Ludwig_variants)

#het_lvl_cols <- colorRamp2(breaks = c(0.05,0.4), 
#                     colors = c("white", "red"))
#lineage_cols <- list(Lineage = c("bulk"="royalblue4", "G11"="magenta3", "B3"="orange", "D2"="yellow", "F4"="grey", "B5"="burlywood4", "B11"="palevioletred2", "A9"="lightseagreen", "D3"="palegreen1", "C7"="lightgoldenrod", "C4"="pink2","C10"="cyan", "B9"="plum1","G10"="steelblue2", "D6"="springgreen4","C9"="red", "mix"="darkgreen"))

#ha <- HeatmapAnnotation(Lineage = SRR_lineage_generation$lineages, col = lineage_cols)
#heatmap_ludwig_variants <- Heatmap(htmp_HET_OR_LOWLVL_nofilt_Ludwig_variants, name = "sqrt(allele frequency)", col = het_lvl_cols, na_col = "white", top_annotation = ha, row_names_gp = gpar(fontsize = 9),column_names_gp = gpar(fontsize = 9))
#png(filename="results/heatmap_Ludwig_variants.png", width = 1024*1.5, height = 1024, units = "px")
#heatmap_ludwig_variants
#dev.off()

}  # for Ludwig_comparison == T

  
##########################  PCA  ##########################################

#all_variants[is.na(all_variants)] <- 0
#all_variants <- all_variants[-71]
#t_all_variants <- data.frame(transpose(all_variants[-1]))
#rownames(t_all_variants) <- colnames(all_variants[-1])
#pca_all <- prcomp(t_all_variants)
#pca_df <- data.frame(pca_all$x)
#pca_plot <- ggplot(data = pca_df, aes(PC1, PC2)) +
#  geom_point()


############ Function to get SRR numbers from S00-- numbers in metadata ##########

# get SRRs for lineage_paths.txt from S00__ number labels on lineage tree (S1d_lineage_tree.png)

get_SRRs_from_Snumbs <- function(Snumb_path){  # See S1d_lineage_tree.png (labelled with S#### sample names). Get list of SRRs to place in lineage_paths.txt (Don't forget to choose and add a name in front of the path list).
  Snumbs_all <- c("bulk","bulk","S0003","S0004","S0005","S0006","S0007","S0008","S0009","S0010","S0011","S0012","S0013","S0014","S0015","S0016","S0017","S0018","S0019","S0020","S0021","S0022","S0023","S0024","S0025","S0026","S0027","S0028","S0029","S0030","S0031","S0032","S0033","S0034","S0035","S0036","S0037","S0038","S0039","S0040","S0041","S0042","S0043","S0044","S0045","S0046","S0047","S0048","S0049","S0050","S0051","S0052","S0053","S0054","S0055","S0056","S0057","S0058","S0059","S0060","S0061","S0062","S0063","S0064","S0065","S0066","S0067","S0068","S0069")
  SRR_path <- list()
  i <- match(Snumb_path, Snumbs_all)
  print(typeof(i))
  for (x in i){
    print(x)
    SRR_path <- paste(SRR_path, SRR_names[x])
  }
  print(SRR_path)
  return(SRR_path)
}

# add Snumbs here. S MUST BE CAPITALIZED. S000 and S0001 not recognised - use "bulk" instead.
Snumb_path <- list("bulk", "bulk", "S0008", "S0030", "S0036", "S0050", "S0058", "S0063", "S0065", "S0066")
SRR_path <- get_SRRs_from_Snumbs(Snumb_path)


  ###################  Modify tables for plots  ######################
## IMPORTANT: add empty rows to SRR_table_list of variant information, so there is one row for every position (for x axis of mutation plots)
barplot_lims <- data.frame(0:16569, rep(1,16570))
colnames(barplot_lims) <- c("Position", "ylimit")
for (i in SRR_names){
  print(i)
  SRR_table_list[[i]]  <- merge(SRR_table_list[[i]], barplot_lims, by.x = "Pos", by.y = "Position", all = T)
  SRR_table_list_HET_OR_LOWLVL_validated[[i]]  <- merge(SRR_table_list_HET_OR_LOWLVL_validated[[i]], barplot_lims, by.x = "Pos", by.y = "Position", all = T)
  SRR_table_list_HET_OR_LOWLVL_nofilt[[i]]  <- merge(x = SRR_table_list_HET_OR_LOWLVL_nofilt[[i]], y = barplot_lims, by.x = "Pos", by.y = "Position", all = T)
  SRR_table_list_HET_OR_LOWLVL[[i]]  <- merge(x = SRR_table_list_HET_OR_LOWLVL[[i]], y = barplot_lims, by.x = "Pos", by.y = "Position", all = T)
}


######################   Exploratory Mutation Plots   ##########################

# only run whole section if exploratory_plots is TRUE:
if (exploratory_plots) {
  
# function to return monotonic values for second y axis transformation (coverage)
f <- function(y){
  log_max <- log2(third_y_lim_maxcoverage_qfilt)
  if (y<=(1/third_y_lim_maxcoverage_qfilt)){
    mono_y <- 2^((y*log_max))
  } else {
  mono_y <- 2^(y*log_max)
  }
  return(mono_y)
}
#f(0.8)

  ## Combine figures by lineage ##

# For each path specified (per line in lineage_paths.txt): 
#   for each SRR in lineage path: 
#     create plot,
#   combine plots on top of each other and save to file

for (p in paths){
    if (p[[1]] == "#"){
    print("skipping comment line...")
    next
    }
  plots_in_lineage <- list()
  n=0
  
  for (SRR in p){
# Skip Lineage path name (1st in character vector of paths[[p]] )
    n=n+1
    print(SRR)
    print(typeof(SRR))
    if (n==1){
      next
    }
# stop where the variants of interest positions are listed on the line (for 
# mutation load plots below)
    if (SRR == "VARIANTS_OF_INTEREST") {
      print("Reached VARIANTS_OF_INTEREST for this lineage")
      break
    }
# Colour according to filter: PASS, STRAND_BIAS, or BLACKLISTED, as not all 
# plots contain blacklisted variants.
    print(SRR)
# Create individual plot:    
    plots_in_lineage[[SRR]] <- ggplot(data = SRR_table_list_HET_OR_LOWLVL_validated[[SRR]], aes(Pos, VariantLevel)) + 
      geom_col(width = 0.9, aes(colour = factor(Filter))) + 
      scale_color_manual(values = c("PASS" = "dodgerblue3",
                                    "STRAND_BIAS"="red",
                                    "BLACKLISTED"="black")) +
      #geom_point(aes(colour = factor(Filter)), size = 0.8) +
      theme_minimal() + 
      ylab(SRR) +
      theme(axis.text.x = element_text(),
            axis.ticks.x = element_line(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 8),
            legend.position = "none",
            plot.margin = margin(t=0.1, r=0.1, b=0.1, l=0.1, "cm"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) +
      geom_text_repel(aes(label = Pos), size = 1.5, nudge_y = 0.05, label.padding = 0.03, box.padding = 0.03, max.overlaps = 13) +
      scale_x_continuous(breaks = seq(0, 16569, by = 1000)) +
      scale_y_continuous(breaks = seq(0, 1.1, by = 0.2), sec.axis = sec_axis(~f(.), name = "log2 coverage", breaks = waiver(), labels = scales::comma)) +
      geom_hline(yintercept=0, size = 0.2) +
      # coverage plot overlay
      geom_line(data = depths_qfilt, aes(Pos, (log2(depths_qfilt[[i]]))/(log2(third_y_lim_maxcoverage_qfilt))), alpha=0.7, size = 0.15) +
      coord_cartesian(ylim = c(0, 1))
    
  }

# combine all the plots in the lineage path: stack on top of each other.
  lineage_plot <- ggarrange(plots = plots_in_lineage, nrow = length(plots_in_lineage), align = "hv")

# Add x and y labels to grid of stacked plots
  y.grob <- textGrob("Variant Level", 
                       gp=gpar(col="black", fontsize=12), rot=90)
  x.grob <- textGrob("Position in mitochondrial genome", 
                     gp=gpar(col="black", fontsize=12))
  lab_lineage_grob <- arrangeGrob(lineage_plot, left = y.grob, bottom = x.grob)
  
# save plot
  file_string <- paste0("results/",p[[1]],"validated_nofilt.png")
  px_height <- 1.2*length(plots_in_lineage)+0.8
  ggsave(file=file_string, plot=lab_lineage_grob, width = 8, height = px_height, units = "in")
}





   ###  Repeat for SRR_table_list_HET_OR_LOWLVL_nofilt  ###

# More easy to select specific position to plot mutation load profiles using only heteroplasmic/low-level variants.
# Unfiltered as some variants seem to switch filter status between generations

# add empty rows to SRR_table_list of variant information, so there is one row for every position (for x axis of mutation plots)
# ^ ALREADY DONE - JUST BEFORE LARGE EXPLORATORY PLOTS
#for (i in SRR_names){
#  print(i)
#  SRR_table_list_HET_OR_LOWLVL_nofilt[[i]]  <- merge(x = SRR_table_list_HET_OR_LOWLVL_nofilt[[i]], y = barplot_lims, by.x = "Pos", by.y = "Position", all = T)
#}

#print("length of SRR_table_list_HET_OR_LOWLVL_nofilt:")
#print(length(SRR_table_list_HET_OR_LOWLVL_nofilt))
#print("structure of SRR 80 in SRR_table_list_HET_OR_LOWLVL_nofilt[[SRR 80]]")
#print(str(SRR_table_list_HET_OR_LOWLVL_nofilt[["SRR7245880"]]))
#print(nrow(SRR_table_list_HET_OR_LOWLVL_nofilt[["SRR7245880"]]$Pos))
#print(levels(SRR_table_list_HET_OR_LOWLVL_nofilt[["SRR7245880"]]$Filter))




## Combine figures by lineage ##

# For each path specified (per line in lineage_paths.txt): 
#   for each SRR in lineage path: 
#     create plot,
#   combine plots on top of each other and save to file

for (p in paths){
  if (p[[1]] == "#"){
    print("skipping comment line...")
    next
  }
  plots_in_lineage <- list()
  n=0
  
  for (SRR in p){
    # Skip Lineage path name (1st in character vector of paths[[p]] )
    n=n+1
    print(SRR)
    print(typeof(SRR))
    if (n==1){
      next
    }
    # stop where the variants of interest positions are listed on the line (for 
    # mutation load plots below)
    if (SRR == "VARIANTS_OF_INTEREST") {
      print("Reached VARIANTS_OF_INTEREST for this lineage")
      break
    }
    # Colour according to filter: PASS, STRAND_BIAS, or BLACKLISTED, as not all 
    # plots contain blacklisted variants.
    print(SRR)
    # Create individual plot:    
    plots_in_lineage[[SRR]] <- ggplot(data = SRR_table_list_HET_OR_LOWLVL_validated[[SRR]], aes(Pos, VariantLevel)) + 
      geom_col(width = 0.9, aes(colour = factor(Filter))) + 
      scale_color_manual(values = c("PASS" = "dodgerblue3",
                                    "STRAND_BIAS"="red",
                                "BLACKLISTED"="black")) + 
      #geom_point(aes(colour = factor(Filter)), size = 0.8) +
      theme_minimal() + 
      ylab(SRR) +
      theme(axis.text.x = element_text(),
            axis.ticks.x = element_line(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 8),
            legend.position = "none",
            plot.margin = margin(t=0.3, r=0.1, b=0.1, l=0.1, "cm"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) +
      geom_text_repel(aes(label = Pos), size = 1.5, nudge_y = 0.05, label.padding = 0.03, box.padding = 0.03, max.overlaps = 13) +
      scale_x_continuous(breaks = seq(0, 16569, by = 1000)) +
      scale_y_continuous(breaks = seq(0, 1.1, by = 0.2)) +
      geom_line(data = depths_qfilt, aes(Pos, (log2(depths_qfilt[[i]]))/(log2(third_y_lim_maxcoverage_qfilt))), alpha=0.4, size = 0.15) + # coverage track
      geom_hline(yintercept=0, size = 0.2) +
      coord_cartesian(ylim = c(0, 1))
  }
  
  # combine all the plots in the lineage path: stack on top of each other.
  lineage_plot <- ggarrange(plots = plots_in_lineage, nrow = length(plots_in_lineage), align = "hv")
  
  # Add x and y labels to grid of stacked plots
  y.grob <- textGrob("Variant Level", 
                     gp=gpar(col="black", fontsize=12), rot=90)
  x.grob <- textGrob("Position in mitochondrial genome", 
                     gp=gpar(col="black", fontsize=12))
  lab_lineage_grob <- arrangeGrob(lineage_plot, left = y.grob, bottom = x.grob)
  
  # save plot
  file_string <- paste0("results/",p[[1]],"_HET_OR_LOWLVL_nofilt.png")
  px_height <- 500*length(plots_in_lineage)+370
  ggsave(file=file_string, plot=lab_lineage_grob, width = 3600, height = px_height, units = "px")
}


#remove(SRR_table_list_HET_OR_LOWLVL_nofilt)




###  Repeat for SRR_table_list_HET_OR_LOWLVL  ###
# More easy to select specific position to plot mutation load profiles using only heteroplasmic/low-level variants.
# Filtered

# add empty rows to SRR_table_list of variant information, so there is one row for every position (for x axis of mutation plots)
# ^ ALREADY DONE JUST BEFORE LARGE EXPLORATORY PLOTS
#for (i in SRR_names){
#  print(i)
#  SRR_table_list_HET_OR_LOWLVL[[i]]  <- merge(x = SRR_table_list_HET_OR_LOWLVL[[i]], y = barplot_lims, by.x = "Pos", by.y = "Position", all = T)
#  #SRR_table_list[[i]]  <- merge(x = SRR_table_list[[i]], y = barplot_lims, by.x = "Pos", by.y = "Position", all = T)
#}

## Combine figures by lineage ##

# For each path specified (per line in lineage_paths.txt): 
#   for each SRR in lineage path: 
#     create plot,
#   combine plots on top of each other and save to file

for (p in paths){
  if (p[[1]] == "#"){
    print("skipping comment line...")
    next
  }
  plots_in_lineage <- list()
  n=0
  
  for (SRR in p){
    # Skip Lineage path name (1st in character vector of paths[[p]] )
    n=n+1
    print(SRR)
    print(typeof(SRR))
    if (n==1){
      next
    }
    # stop where the variants of interest positions are listed on the line (for 
    # mutation load plots below)
    if (SRR == "VARIANTS_OF_INTEREST") {
      print("Reached VARIANTS_OF_INTEREST for this lineage")
      break
    }
    # Colour according to filter: PASS, STRAND_BIAS, or BLACKLISTED, as not all 
    # plots contain blacklisted variants.
    print(SRR)
    # Create individual plot:    
    plots_in_lineage[[SRR]] <- ggplot(data = SRR_table_list_HET_OR_LOWLVL[[SRR]], aes(Pos, VariantLevel)) + 
      geom_col(width = 0.9, aes(colour = factor(Filter))) + 
      scale_color_manual(values = c("PASS" = "dodgerblue3",
                                    "STRAND_BIAS"="red",
                                    "BLACKLISTED"="black")) + 
      #geom_point(aes(colour = factor(Filter)), size = 0.8) +
      theme_minimal() + 
      ylab(SRR) +
      theme(axis.text.x = element_text(),
            axis.ticks.x = element_line(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 8),
            legend.position = "none",
            plot.margin = margin(t=0.3, r=0.1, b=0.1, l=0.1, "cm"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) +
      geom_text_repel(aes(label = Pos), size = 1.5, nudge_y = 0.05, label.padding = 0.03, box.padding = 0.03, max.overlaps = 13) +
      scale_x_continuous(breaks = seq(0, 16569, by = 1000)) +
      scale_y_continuous(breaks = seq(0, 1.1, by = 0.2)) +
      geom_line(data = depths_qfilt, aes(Pos, (log2(depths_qfilt[[i]]))/(log2(third_y_lim_maxcoverage_qfilt))), alpha=0.7, size = 0.15) + # coverage across genome
      geom_hline(yintercept=0, size = 0.2) +
      coord_cartesian(ylim = c(0, 1))
    
  }
  
  # combine all the plots in the lineage path: stack on top of each other.
  lineage_plot <- ggarrange(plots = plots_in_lineage, nrow = length(plots_in_lineage), align = "hv")
  
  # Add x and y labels to grid of stacked plots
  y.grob <- textGrob("Variant Level", 
                     gp=gpar(col="black", fontsize=12), rot=90)
  x.grob <- textGrob("Position in mitochondrial genome", 
                     gp=gpar(col="black", fontsize=12))
  lab_lineage_grob <- arrangeGrob(lineage_plot, left = y.grob, bottom = x.grob)
  
  # save plot
  file_string <- paste0("results/",p[[1]],"_HET_OR_LOWLVL.png")
  px_height <- 500*length(plots_in_lineage)+370
  ggsave(file=file_string, plot=lab_lineage_grob, width = 3600, height = px_height, units = "px")
}

}  #(if (exploratory_plots is TRUE){}) bracket


  ###############  Position-specific mutation load plots  #####################
  #####################  for variants of interest  ############################

if (position_specific_plots == TRUE) {
#long_lineage_cols <- unique(df[, c("zone", "color.names")])
#pos_of_interest <- 4769

# For each specified path, move over and record SRRs until VARIANTS_OF_INTEREST 
# are listed. Then for each VARIANT_OF_INTEREST position, plot line graph of 
# how variant's mutation load changes across the lineage path.

#generation_axis_labs <- c("bulk","bulk","G1","G2","G3","G4","G5","G6","G7","G8")
lin_mut_load_change <- data.frame(matrix(nrow = 0, ncol = 8))
colnames(lin_mut_load_change) <- c("SRR", "Generation", "Generation_labs", "Pos", "Lineage", "Lineage_group", "VariantLevel", "Coverage")
for (p in paths){
  if (p[[1]] == "#"){
    print("skipping comment line...")
    next
  } 
  if (str_detect(p[[1]], 'LUDWIG')){
    print(paste0("skipping LUDWIG line: ", p[[1]]))
    next
  }
  SRRs_in_path <- list()
  index=0
  at.positions = F
  # open dev for PDF for lineage:
  file_string <- paste0("results/",p[[1]],"_selected_pos_of_interest.pdf")
  pdf(file_string, onefile = TRUE)

  for (string in p){
    print(string)
    index=index+1
    #print(index)
    #print(at.positions)
# skip lineage name
    if (index==1){
      next
    }
# From lineage_paths.txt: add SRR to list, until the positions of variants of interest are listed
# instead. This is indicated by "VARIANTS_OF_INTEREST" instead of SRR name, 
# and followed by the positions of variants of interest.
# eg. 
#    LINEAGE_PATH_NAME SRR1 SRR2 SRR3 VARIANTS_OF_INTEREST 1495 12788

    if (string == "VARIANTS_OF_INTEREST") {
      print("Reached VARIANTS_OF_INTEREST for this lineage")
      at.positions <- T
      #print(paste("n =",n))
      n=0
      next
    }
    paste(string, at.positions)
    if (at.positions == F){
      SRRs_in_path <- c(SRRs_in_path, string)
    }
    if (at.positions == T){
      print(paste("n=",n))
      pos_of_interest <- as.numeric(string)
      #print(paste0("Plotting position: ", pos_of_interest, ", for lineage path: ", p[[1]]))
      
      # make new data frame for new variant position
      mut_load_change <- data.frame(matrix(nrow = length(SRRs_in_path), ncol = 8))
      colnames(mut_load_change) <- c("SRR", "Generation", "Generation_labs", "Pos", "Lineage", "Lineage_group", "VariantLevel", "Coverage")
      #print(paste("SRRs_in_path: ", SRRs_in_path))
      for (SRR_name in SRRs_in_path){
        #for (base in c("A","C","G","T")) {
          #if (SRR_table_list[[SRR_name]]$Variant[which(SRR_table_list[[SRR_name]]$Pos == pos_of_interest)]) == base)
          n=n+1
          print(SRR_name)
          mut_load_change$SRR[n] <- SRR_name
          mut_load_change$Generation[n] <- n-1
          mut_load_change$Generation_labs[n] <- paste(SRR_lineage_generation$generation_axis_labs[SRR_lineage_generation$SRR_names==SRR_name])
          mut_load_change$VariantLevel[n] <- SRR_table_list[[SRR_name]]$VariantLevel[which(SRR_table_list[[SRR_name]]$Pos == pos_of_interest)] #SRR_table_list[[SRR_name]]$VariantLevel[pos_of_interest+1]
          mut_load_change$Coverage[n] <- depths_qfilt[[SRR_name]][pos_of_interest]
          mut_load_change$Pos[n] <- pos_of_interest
          mut_load_change$Lineage <- p[[1]]
          mut_load_change$Lineage_group <- paste(SRR_lineage_generation$lineages[SRR_lineage_generation$SRR_names==SRR_name])
      }
      last_SRR <- SRRs_in_path[length(SRRs_in_path)]
      lin <- as.character(SRR_lineage_generation$lineages[SRR_lineage_generation$SRR_names==last_SRR])
      lin_col <- lineage_cols[[lin]]
      
      # Plot for new variant position
      mut_load_change[is.na(mut_load_change)] <- 0
      plot_title <- paste0(p[[1]],": ", pos_of_interest)
      mut_plot <- ggplot(data = mut_load_change, aes(x=Generation,y=VariantLevel)) +
        geom_line(colour = lin_col) +
        geom_text_repel(aes(label = Coverage, size = 13)) +
        geom_point(colour = lin_col) +
        theme_minimal() +
        theme(plot.background = element_rect(fill = "white",
                                colour = "white"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               axis.line = element_line(colour = "black"), legend.position = "none", text = element_text(size=20)) +
        ggtitle(plot_title) +
        scale_x_continuous(breaks = mut_load_change$Generation, labels = mut_load_change$Generation_labs) +
        expand_limits(y = c(0,0.01)) +
        geom_hline(yintercept=0.01, size = 1,linetype="dotted", colour = "red")+
        labs(y = "Allele Frequency")
      
      # save plot
      print(mut_plot)
      file_string <- paste0("results/",p[[1]],"_pos_",pos_of_interest, ".png")
      #ggsave(file=file_string, plot=mut_plot)
      n=0
      #print("n reset")
      #print(lin_mut_load_change)
      # Append mut_load_change (for one position) to lin_mut_load_change (all positions_of_interest, all lineages)
      lin_mut_load_change <- rbind(lin_mut_load_change, mut_load_change)
      }
    else {  # if VARIANTS_OF_INTEREST hasn't been reached (at.positions=F)
      next
    }
  }
  
  # Plot ALL interesting variant positions on one graph per lineage
  mut_load_change[is.na(mut_load_change)] <- 0
  plot_title <- paste0(p[[1]], ": ALL positions of interest")
  mut_plot <- ggplot(data = lin_mut_load_change[lin_mut_load_change$Lineage == p[[1]], ], aes(x=Generation,y=VariantLevel,group=Pos,colour=Lineage_group)) +
    geom_line() +
    scale_colour_manual(values=lineage_cols, breaks=colnames(lineage_cols)) +
    geom_text_repel(aes(label = Coverage, size = 13)) +
    geom_point(aes(colour=lin_col)) +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white",
                                         colour = "white"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"), legend.position = "none", text = element_text(size=20)) +
    ggtitle(plot_title) +
    scale_x_continuous(breaks = lin_mut_load_change$Generation, labels = lin_mut_load_change$Generation_labs) +
    expand_limits(y = c(0,0.01)) +
    geom_hline(yintercept=0.01, size = 1,linetype="dotted", colour = "red") +
    labs(y = "Allele Frequency")
  
  # save plot
  print(mut_plot)
  file_string <- paste0("results/",p[[1]],"_ALL_pos_of_interest.png")
  #ggsave(file=file_string, plot=mut_plot)
  n=0
  #print("n reset")
  dev.off()
}


###################  Mutation load plots:  #########################
#############  each single position of interest  ###################
############  across multiple lineages on one plot  ################

# With any positions of interest chosen in lineage_paths.txt that occur in more than one lineage
cross_lineage_positions_interest <- unique(lin_mut_load_change[ ,c('Pos', 'Lineage')]) %>% 
  filter(!str_detect(Lineage, 'LUDWIG')) %>% count(Pos) #%>% filter(n>=2)
cross_lineage_positions_interest <- merge((lin_mut_load_change %>% filter(!str_detect(Lineage, 'LUDWIG'))), 
                                          cross_lineage_positions_interest, by = 'Pos') %>% filter(n>=2)
nrow(cross_lineage_positions_interest)
nrow(lin_mut_load_change)

pdf("results/selected_pos_of_intrest_across_lineages.pdf", onefile = TRUE)
for (pos_of_interest in unique(cross_lineage_positions_interest$Pos)){
  # Plot for new variant position
  #mut_load_change[is.na(mut_load_change)] <- 0
  cross_lineages <- unique(cross_lineage_positions_interest[ ,c('Pos', 'Lineage')]) %>% filter(Pos == pos_of_interest)
  #print(c(pos_of_interest, cross_lineages$Lineage))
  # Plot ALL interesting variant positions on one graph per lineage
  mut_load_change[is.na(mut_load_change)] <- 0
  plot_title <- paste0("Position: ", pos_of_interest, " across lineages")
  mut_plot <- ggplot(data = cross_lineage_positions_interest[cross_lineage_positions_interest$Pos == pos_of_interest, ], 
                     aes(x=Generation,y=VariantLevel,group=Lineage,color=Lineage_group)) +
    geom_line() +
    scale_colour_manual(values=lineage_cols, name="Lineage") +
    geom_text_repel(aes(label = Coverage), size = 3.5, show.legend = FALSE) +
    geom_point() +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white", colour = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"), text = element_text(size=15), legend.text = element_text(size = 9), legend.title = element_text(size = 10)) +
    ggtitle(plot_title) +
    scale_x_continuous(breaks = lin_mut_load_change$Generation, labels = lin_mut_load_change$Generation_labs) +
    expand_limits(y = c(0,0.01)) +
    geom_hline(yintercept=0.01, size = 1,linetype="dotted", colour = "red") +
    labs(y = "Allele Frequency")
    
  # save plot
  #file_string <- paste0("results/pos_", pos_of_interest, "_across_lineages.png")
  #ggsave(file=file_string, plot=mut_plot)
  n=0
  print(mut_plot)
}
dev.off()


# png("myplot.png")  # opens file to write to
# plot somthing
# dev.off()

# multiple plots into one file
# pdf()


#################  Mutation load plots across lineages ######################
#####################  for ALL LINEAGE VALIDATED ############################

  ## Plot all validated variants that occur in more than one lineage on one graph ##
## With all validated variant calls that occur in more than one lineage
# First make table like lin_mut_load_change but with all validated positions, 
#if they're validated in the lineage, then filter for positions which are valid in more than one lineage afterwards.
# "SRR", "Generation", "Generation_labs", "Pos", "Lineage", "Lineage_group", "VariantLevel", "Coverage"

# Create table:
lin_mut_load_change_lin_val <- data.frame(matrix(nrow = 0, ncol = 8))
next_row <- data.frame(matrix(nrow = 1, ncol = 8))
colnames(lin_mut_load_change_lin_val) <- c("SRR", "Generation", "Generation_labs", "Pos", "Lineage", "Lineage_group", "VariantLevel", "Coverage")
colnames(next_row) <- c("SRR", "Generation", "Generation_labs", "Pos", "Lineage", "Lineage_group", "VariantLevel", "Coverage")
n=0
for (p in paths){
  if (p[[1]] == "#"){
    print("skipping comment line...")
    next
  }
  if (str_detect(p[[1]], 'LUDWIG')){
    print(paste0("skipping LUDWIG line: ", p[[1]]))
    next
  }
  SRRs_in_path <- list()
  index=0
  at.positions = F
  
  for (string in p){
    print(string)
    index=index+1
    print(index)
    print(at.positions)
    # skip lineage name
    if (index==1){
      next
    }
    # From lineage_paths.txt: add SRR to list, until the positions of variants of interest are listed
    # instead. This is indicated by "VARIANTS_OF_INTEREST" instead of SRR name, 
    # and followed by the positions of variants of interest. eg.:
    # LINEAGE_PATH_NAME SRR1 SRR2 SRR3 VARIANTS_OF_INTEREST 1495 12788
    if (string == "VARIANTS_OF_INTEREST") {
      print("Reached VARIANTS_OF_INTEREST for this lineage")
      at.positions <- T
      print(paste("n =",n))
      next
    }
    paste(string, at.positions)
    if (at.positions == F){
      SRRs_in_path <- c(SRRs_in_path, string)
    }
    if (at.positions == T){
      break
    }
  }
  last_SRR <- SRRs_in_path[length(SRRs_in_path)]
  print(paste("n=",n))
  for (row_num in 1:nrow(all_validated_pos_base)){  # do row so multiple poses can
    base <- all_validated_pos_base[row_num, "Base"]
    pos <- all_validated_pos_base[row_num, "Pos"]
    n=0
    print(paste("Position:", pos, "Lineage:", p[[1]]))
    for (SRR_name in SRRs_in_path){
      #if (pos %in% (SRR_table_list_HET_OR_LOWLVL_validated[[SRR_name]] %>% drop_na(VariantLevel) %>% select(Pos))){
        n=n+1
        next_row$SRR <- SRR_name
        next_row$Pos <- pos
        next_row$Generation <- n-1  #as.numeric(SRR_lineage_generation$generation[SRR_names == SRR_name])
        next_row$Generation_labs <- SRR_lineage_generation$generation_axis_labs[SRR_names == SRR_name]
        print("Stage 2")
        #next_row$VariantLevel <- SRR_table_list[[SRR_name]] %>% filter(Pos == pos & Variant == base) %>% select(VariantLevel)
        if (base %in% SRR_table_list[[SRR_name]]$Variant[which(SRR_table_list[[SRR_name]]$Pos == pos)]) { 
          next_row$VariantLevel <- SRR_table_list[[SRR_name]]$VariantLevel[which(SRR_table_list[[SRR_name]]$Pos == pos & SRR_table_list[[SRR_name]]$Variant == base)]
          next_row$Coverage <- SRR_table_list[[SRR_name]]$Variant_AD[which(SRR_table_list[[SRR_name]]$Pos == pos & SRR_table_list[[SRR_name]]$Variant == base)] #SRR_table_list[[SRR_name]]$Coverage[which(SRR_table_list[[SRR_name]]$Pos == pos & SRR_table_list[[SRR_name]]$Variant == base)] 
        }
        else {
          next_row$VariantLevel <- NA
          next_row$Coverage <- NA
        }
        print("Stage 3")
        next_row$Lineage <- p[[1]]
        print("Stage 4")
        next_row$Lineage_group <- SRR_lineage_generation$lineages[SRR_lineage_generation$SRR_names==last_SRR]
        lin_mut_load_change_lin_val <- rbind(lin_mut_load_change_lin_val, next_row)
      #}
    }
  }
}

# filter lin_mut_load_change_lin_val for positions of any validated mutations that occur in more than one lineage.
cross_lineage_positions_lin_val <- unique(lin_mut_load_change_lin_val %>% 
  filter(!str_detect(Lineage, 'LUDWIG')) %>% drop_na(VariantLevel) %>% select(Pos, Lineage)) %>% count(Pos) #%>% filter(n>=2)
cross_lineage_positions_lin_val <- merge(
  (lin_mut_load_change_lin_val %>% filter(!str_detect(Lineage, 'LUDWIG'))), 
  cross_lineage_positions_lin_val, by = 'Pos') %>% filter(n>=2)

# Convert allele frequency (VariantLevel) NAs to 0.00s for plotting
cross_lineage_positions_lin_val[is.na(cross_lineage_positions_lin_val)] <- 0

# Calculate differences between generation 3 and the rest - spike in number and AF of alleles for some positions
# general tests (mapping qualities, base qualities), position specific tests (depths, SP)
gen3_poerpos_stats <- data.frame(Pos=lin_mut_load_change_lin_val$Pos)
not_gen3_SRRs_stats <- data.frame(Pos=lin_mut_load_change_lin_val$Pos)
gen3SRRs <- c()
notgen3SRRs <- c()
for (SRR_name in SRR_names){
  if (SRR_lineage_generation$generation[which(SRR_lineage_generation$SRR_names == SRR_name)] == 3) {
    gen3SRRs <- c(gen3SRRs, SRR_name)
  }
  else {
    notgen3SRRs <- c(notgen3SRRs, SRR_name)
  }
}
gen3_stats <- data.frame(SRR_names=gen3SRRs)
mapqs <- c()
baseqs <- c()
meandepths <- c()
for (SRR_name in gen3_stats$SRR_names){
  # mapq
  mapqs <- c(mapqs, all_coverages_qfilt$meanmapq[which(all_coverages_qfilt$SRRfile == SRR_name)])
  baseqs <- c(baseqs, all_coverages_qfilt$meanbaseq[which(all_coverages_qfilt$SRRfile == SRR_name)])
  meandepths <- c(meandepths, all_coverages_qfilt$meandepth[which(all_coverages_qfilt$SRRfile == SRR_name)])
}
gen3_stats[["mapq"]] <- mapqs
gen3_stats[["baseq"]] <- baseqs
gen3_stats[["meandepth"]] <- meandepths
gen3_stats$isgen3 <- "Gen 3"

notgen3_stats <- data.frame(SRR_names=notgen3SRRs)
mapqs <- c()
baseqs <- c()
meandepths <- c()
for (SRR_name in notgen3_stats$SRR_names){
  # mapq
  mapqs <- c(mapqs, all_coverages_qfilt$meanmapq[which(all_coverages_qfilt$SRRfile == SRR_name)])
  baseqs <- c(baseqs, all_coverages_qfilt$meanbaseq[which(all_coverages_qfilt$SRRfile == SRR_name)])
  meandepths <- c(meandepths, all_coverages_qfilt$meandepth[which(all_coverages_qfilt$SRRfile == SRR_name)])
}
notgen3_stats[["mapq"]] <- mapqs
notgen3_stats[["baseq"]] <- baseqs
notgen3_stats[["meandepth"]] <- meandepths
notgen3_stats$isgen3 <- "not Gen 3"

gen3_stats <- rbind(gen3_stats, notgen3_stats)

wilcox.test(mapq ~ isgen3, data = gen3_stats)
par(mfrow = c(2, 2))
ggplot(gen3_stats, aes(isgen3, mapq)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width=0.1) +
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2, colour = "red")
ggplot(gen3_stats, aes(isgen3, baseq)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width=0.1)+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2, colour = "red")
ggplot(gen3_stats, aes(isgen3, meandepth)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width=0.1)+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2, colour = "red")
ggplot(gen3_stats, aes(isgen3, meandepth)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width=0.1)+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2, colour = "red")

for (pos in lin_mut_load_change$Pos){
  # mapq
  
  # baseq
  # depth
  # soft-clipped
  # Read-position bias
  # mapping quality bias
}


# Plot:
pdf("results/validated_pos_all_noGen3_noSupportReads.pdf", onefile = TRUE)
plot_list <- list()
# for each lineage validated position that occurs more than once: 
n=0
for (pos in sort(unique(lin_mut_load_change_lin_val$Pos))){
  # Plot for new variant position
  n=n+1
  #cross_lineages <- unique(lin_mut_load_change_lin_val[ ,c('Pos', 'Lineage')]) %>% filter(Pos == pos)
  #print(c(pos, cross_lineages$Lineage))
  # Plot ALL interesting variant positions on one graph per lineage
  lin_mut_load_change_lin_val[is.na(lin_mut_load_change_lin_val)] <- 0
  plot_title <- paste0("Position: ", pos, " across lineages")
  mut_plot <- ggplot(data = lin_mut_load_change_lin_val[lin_mut_load_change_lin_val$Pos == pos, ], 
                     aes(x=Generation,y=VariantLevel,group=Lineage,color=Lineage_group)) +
    geom_line() +
    scale_colour_manual(values=lineage_cols, name="Lineage") + #, breaks=colnames(lineage_cols)) +
    geom_text_repel(aes(label = Coverage), size = 3, show.legend = FALSE, segment.color = "black", segment.alpha = 0.5, segment.size = 0.25) +
    #geom_point() +  
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white", colour = "white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               axis.line = element_line(colour = "black"), text = element_text(size=15), legend.text = element_text(size = 9), legend.title = element_text(size = 10)) +
    ggtitle(plot_title) +
    scale_x_continuous(breaks = lin_mut_load_change_lin_val$Generation, labels = lin_mut_load_change_lin_val$Generation_labs) +
    expand_limits(y = c(0,0.01)) +
    geom_hline(yintercept=0.01, size = 1,linetype="dotted", colour = "red") +
    labs(y = "Allele Frequency")
  print(mut_plot)
  #plot_list[[((n-1)%%6 + 1)]] <- mut_plot
  #if (n %% 6 == 0){
  #  plotem <- ggarrange(plots = plot_list, nrow = 3, ncol = 2)
  #  print(plotem)
  #  plot_list <- list()
  #}
  #else if (n == length(unique(cross_lineage_positions_lin_val$Pos))){
  #  plotem <- ggarrange(plots = plot_list, nrow = 3, ncol = 2)
  #  print(plotem)
  #}
  
  # save plot
  #file_string <- paste0("results/validated_pos_", pos, "_across_lins.png")
  #ggsave(file=file_string, plot=mut_plot)
  
  #print(mut_plot)
  #n=0
  #print("n reset")
}
dev.off()

#length(plot_list)
#arranged_plot_list <- marrangeGrob(plot_list, nrow = 3, ncol = 2)
#file_string <- "results/validated_pos_across_lins.pdf"
#ggsave(file=file_string, plot=arranged_plot_list)

}  # end if (position_specific_plots == TRUE)


####################### Replicate correlation plot #######################

# combine AF of replicates
bulk_replicates_all_nofilt <- merge(SRR_table_list_HET_OR_LOWLVL_nofilt$SRR7245880[, c("Pos", "VariantLevel")], SRR_table_list_HET_OR_LOWLVL_nofilt$SRR7245881[, c("Pos", "VariantLevel")], by = "Pos", all = T)
colnames(bulk_replicates_all_nofilt) <- c("Pos", "SRR7245880", "SRR7245881")
bulk_replicates_all_nofilt[is.na(bulk_replicates_all_nofilt)] <- 0


#bulk_replicates_Ludwigs_nofilt <- data.frame(our_Ludwig_variants_nofilt$SRR7245880, our_Ludwig_variants_nofilt$SRR7245881)
#colnames(bulk_replicates_Ludwigs_nofilt) <- c("Bulk_SRR7245880", "Bulk_SRR7245881")
#bulk_replicates_Ludwigs_nofilt[is.na(bulk_replicates_Ludwigs_nofilt)] <- 0

# plot sqrt AF
bulk_rep_corr_plot_all_nofilt <- ggplot(bulk_replicates_all_nofilt, aes(sqrt(SRR7245880),sqrt(SRR7245881))) +
  geom_point() +
  #geom_point(data=bulk_replicates_Ludwigs_nofilt, aes(sqrt(Bulk_SRR7245880),sqrt(Bulk_SRR7245881), colour = "red")) +
  labs(x ="sqrt(AF) of Bulk replicate SRR7245880", y ="sqrt(AF) of Bulk replicate SRR7245881") +
  scale_x_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) +
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), text=element_text(size=13), legend.position = "none") +
  geom_abline(intercept = 0, slope = 1, alpha = 0.4)
ggsave(file="results/bulk_replicate_scatter.png", plot=bulk_rep_corr_plot_all_nofilt, width = 4, height = 4, units = "in")


#Fig 3?
comparison_plots <- list(corr_plot_all, bulk_rep_corr_plot_all_nofilt)
comparison_plots <- ggarrange(plots = comparison_plots, ncol = 2, align = "hv")
ggsave(file="results/comparison_plots.png", plot=comparison_plots, width = 8, height = 4, units = "in")
#comparison_plots_and_heatmap <- plot_grid(
#  heatmap_ludwig_variants, comparison_plots,
#  labels = "AUTO", ncol = 1
#)


  #####################  Vectors of allele read proportions ########################
 ##  (MOVED: created with bcf_SRR_table_lists)

# for each nucleotide ACGT (solves multiallelic problem of more than one AF per genomic position)
# (weighted?) supporting reads over total reads

  ## Test AF_vector ##
##> SRR_table_list$SRR7245880[SRR_table_list$SRR7245880$Pos == 50,]
## Pos Filter Ref Variant Variant_AD Variant_ADF Variant_ADR Coverage VariantLevel Ref_AD Ref_ADF Ref_ADR  RefLevel Type ylimit
##56  50   PASS   T       G          2           2           0     2112 0.0009469697   2110    2110       0 0.9990530    2      1
##57  50   PASS   T       A          1           1           0     2111 0.0004737091   2110    2110       0 0.9995263    2      1
#
#SRR_880_Gs <- AF_vectorise(SRR_table_list$SRR7245880, "G")
#SRR_880_As <- AF_vectorise(SRR_table_list$SRR7245880, "A")
#SRR_880_Gs[50]
##[1] 0.0009469697
#SRR_880_As[50]
##0.0004737091


## return vectors for each lineage
vectors_per_lineage <- list()
for (p in paths){
  if (p[[1]] == "#"){
    print("skipping comment line...")
    next
  }
  if (str_detect(p[[1]], 'LUDWIG')){
    print(paste0("skipping LUDWIG line: ", p[[1]]))
  }
  SRRs_in_path <- list()
  index=0
  at.positions = F
  
  for (string in p){
    print(string)
    index=index+1
    print(index)
    print(at.positions)
    # skip lineage name
    if (index==1){
      next
    }
    # From lineage_paths.txt: add SRR to list, until the positions of variants of interest are listed
    # instead. This is indicated by "VARIANTS_OF_INTEREST" after the last SRR name, 
    # and followed by the positions of variants of interest. eg.:
    # LINEAGE_PATH_NAME SRR1 SRR2 SRR3 VARIANTS_OF_INTEREST 1495 12788
    if (string == "VARIANTS_OF_INTEREST") {
      print("Reached VARIANTS_OF_INTEREST for this lineage")
      at.positions <- T
      print(paste("n =",n))
      next
    }
    paste(string, at.positions)
    if (at.positions == F){
      SRRs_in_path <- c(SRRs_in_path, string)
    }
    if (at.positions == T){
      break
    }
  }
  last_SRR <- SRRs_in_path[length(SRRs_in_path)]
  
  vectors_per_lineage[[ p[[1]] ]] <- list()
  for (n in seq.int(1,length(SRRs_in_path),1)){
    vectors_per_lineage[[ p[[1]] ]][n] <- vector_list[ SRRs_in_path[n] ]
  }
    
}

vector_list

