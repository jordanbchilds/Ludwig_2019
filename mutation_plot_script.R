

  ## NOTE: this script should be called via the bash script "plot_mutations.sh" ##

# set working directory
#args <- commandArgs(trailingOnly = T)
#print(args)
#setwd(args[1])
setwd("/home/thomas/Documents/Research_proj/Ludwig_2019/")
  





  ## Load packages ##
# Check if packages are installed, install to .R_local_lib in working directory if needed.
#local_lib_path <- paste0(args[1],"/.R_local_lib/")
#print(local_lib_path)
#.libPaths(c(local_lib_path, .libPaths()))

packages <- c("tidyr","ggplot2","gridExtra","ggrepel","egg","grid")
lapply(packages, FUN = function(i) {
  if (!require(i, character.only = TRUE)) {
    install.packages(i, dependencies = TRUE, lib = local_lib_path, repos="https://www.stats.bris.ac.uk/R/")
    library(i, character.only = TRUE)
    }
  }
)



  ## Read SRR files ##
filenames <- list.files("vcf/", pattern="*_annotated.txt")
# Create list of data frame names without the ".txt" part 
SRR_names <-substr(filenames,1,10)

  ## Read coverage files ##
depths <- read.table("coverages/depths.txt", sep = "\t", header = F, stringsAsFactors = T)
depths_qfilt <- read.table("coverages/depths_qfilt.txt", sep = "\t", header = F, stringsAsFactors = T)
colnames(depths_qfilt) <- c("chr", "Pos", SRR_names)
colnames(depths) <- c("chr", "Pos", SRR_names)

  ## Read lineage path/s ##
# path through a lineage specified and read from 'lineage_paths.txt'
paths <- list()
paths <- as.list(strsplit(readLines("lineage_paths.txt"), " "))


  ####################  Pre-alignment plots and tables #########################

raw_sample_info <- read.csv("SraRunTable_1.csv", header = T)
pre_multiqc <- read.table("multiQC/group_SRP149534_multiQC_report_data/multiqc_general_stats.txt", header = T) 
colnames(pre_multiqc) <- c("SRR_sample", "percent_dup", "percent_gc", "sequence_lengths", "percent_fails", "num_seqs")

pre_dup_hist <- ggplot(data = pre_multiqc, aes(percent_dup)) +
  geom_histogram(fill = "dodgerblue3", colour = "black", binwidth = 2) +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  labs(x ="% duplicate reads") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

pre_num_seqs_hist <- ggplot(data = pre_multiqc, aes(num_seqs)) +
  geom_histogram(fill = "dodgerblue3", colour = "black", bins = 30) +
  scale_y_continuous(expand = expansion(mult = c(0, .05))) +
  scale_x_continuous(breaks = c(5000000,10000000,15000000,20000000,25000000), labels = scales::comma) +
  labs(x ="Number of reads") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
file_string <- "results/pre_alignment_num_seqs_hist.png"
ggsave(file=file_string, plot=pre_num_seqs_hist)




  #####################  Coverage plots and stats  ###################### (TODO)
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


    ###################### Post alignment quality ########################

all_coverages_qfilt <- read.csv("coverages/all_coverages_qfilt.txt", header = T, sep = "\t")
all_coverages_qfilt$SRRs <- SRR_names
all_coverages <- read.csv("coverages/all_coverages.txt", header = T, sep = "\t")
all_coverages$SRRs <- SRR_names

all_coverages$calculated_mean <- as.numeric(lapply(depths[,3:ncol(depths)], mean))
all_coverages$calculated_sd <- as.numeric(lapply(depths[,3:ncol(depths)], sd))
all_coverages_qfilt$calculated_mean <- as.numeric(lapply(depths_qfilt[,3:ncol(depths_qfilt)], mean))
all_coverages_qfilt$calculated_sd <- as.numeric(lapply(depths_qfilt[,3:ncol(depths_qfilt)], sd))


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
sample_coverage_plot_qfilt <- ggplot(data = all_coverages_qfilt, aes(SRRs, calculated_mean, calculated_sd)) +
  geom_col(colour = "black", fill = "dodgerblue3") +
  geom_errorbar(aes(ymin=calculated_mean-calculated_sd, ymax=calculated_mean+calculated_sd), width=0) +
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
  theme(axis.text.x = element_text(angle = 45, vjust=1.05, hjust = 1.0),
        axis.ticks.x = element_line(),
        panel.background = element_rect(fill = "white")) 

# unfiltered read plots: coverage (x axis: SRR), no.reads, base quality, mapping quality
sample_coverage_plot <- ggplot(data = all_coverages, aes(SRRs, calculated_mean, calculated_sd)) +
  geom_col(colour = "black", fill = "orange") +
  geom_errorbar(aes(ymin=calculated_mean-calculated_sd, ymax=calculated_mean+calculated_sd), width=0) +
  scale_y_continuous(trans='log2', expand = expansion(mult = c(0, .1))) +
  theme(axis.text.x = element_text(angle = 45, vjust=1.05, hjust = 1.0),
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
  geom_col(colour = "black", fill = "orange") +
  scale_y_continuous(expand = expansion(mult = c(0, .1)), labels = scales::comma) +
  theme(axis.text.x = element_text(angle = 45, vjust=1.05, hjust = 1.0),
        axis.ticks.x = element_line(),
        panel.background = element_rect(fill = "white"))




  #######################  VARIANT DATAFRAMES  ####################

SRR_table_list <- list()  # All 
SRR_table_list_PASS <- list()  # Filtered
SRR_table_list_HET_OR_LOWLVL <- list()  # Filtered and only heteroplasmic/low level
SRR_table_list_HET_OR_LOWLVL_nofilt <- list()

#threshold <- 0.05
# Load files into list of data.frames
for(i in SRR_names){
  filepath <- file.path("vcf",paste(i,"_annotated.txt",sep=""))
  SRR_table_list[[i]] <- read.table(filepath, sep = "\t", header = T, stringsAsFactors = T)
  # remove variants where the reference = "N". multiallelic, but only major level plotted - stacked bars
  SRR_table_list[[i]] <- SRR_table_list[[i]][!(SRR_table_list[[i]]$Ref=="N"),]
  # subset for variants which passed filter
  SRR_table_list_PASS[[i]] <- subset(SRR_table_list[[i]], Filter == "PASS")
  # subset for "HET_OR_LOWLVL" variants (Heteroplasmic or low-level variant)
  SRR_table_list_HET_OR_LOWLVL[[i]] <- subset(SRR_table_list_PASS[[i]], Type == 2)
  # no filter to see if variants are filtered differently in different samples across lineage paths
  SRR_table_list_HET_OR_LOWLVL_nofilt[[i]] <- subset(SRR_table_list[[i]], Type ==2)
  }
print("Before merging structure of SRR_table_list")
print(str(SRR_table_list[["SRR7245880"]]))


  #################### Variant calling stats ########################

variant_stats <- data.frame(matrix(nrow = length(SRR_table_list), ncol = 8))
colnames(variant_stats) <- c("SRR","No.Variants", "No.Unfiltered_Variants", "No.het", "No.hom","No.transition", "No.transversion", "No.missense")
variant_stats$SRR <- SRR_names
#variant_stats$

#for (i in SRR_names){
#  variant_stats$No.Variants[[i]] <- nrow(SRR_table_list[[i]])
#}


bulk_variant_pos80 <- data.frame(SRR_table_list_HET_OR_LOWLVL$SRR7245880$Pos)#[SRR_table_list$SRR7245880$Type==2])
bulk_variant_pos81 <- data.frame(SRR_table_list_HET_OR_LOWLVL$SRR7245881$Pos)#[SRR_table_list$SRR7245881$Type==2])
bulk_variant_pos <- merge(bulk_variant_pos80, bulk_variant_pos81, by=1, all=T)
colnames(bulk_variant_pos) <- "Bulk_Variants"
write.csv(bulk_variant_pos, file = "results/bulk_variant_positions.csv", quote = F)

   ####  all_variants_in_path
# merge to list all variants in lineage path
# replace postition with variant level.
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
    
    SRR_pos_level <- data.frame(SRR_table_list[[SRR]]$Pos, SRR_table_list[[SRR]]$VariantLevel)
    print(colnames(SRR_pos_level))
    colnames(SRR_pos_level) <- c("Pos", paste0(SRR,"_variant_lvl"))
    all_variants_in_lineage <- merge(all_variants_in_lineage,SRR_pos_level, by = "Pos", all = T)
    
  }
file_string <- paste0("results/",p[[1]],"_all_variants.csv")
write.csv(all_variants_in_lineage,file = file_string, quote = F)
print("table of bulk variants in lineage path saved in 'results/'")
}



   ### Compare positions of our variants with Ludwig's positions ###
# read in Ludwigs variants, and variant level in each sample
Ludwig_variants <- read.csv("LUDWIG_TF1_clones_ATAC_alleleFrequencies.csv", header = T)
colnames(Ludwig_variants)[1] <- "Ludwig_variant_positions"
Ludwig_variants$tobecombined_Pos <- Ludwig_variants$Ludwig_variant_positions

# make data fram of positions of all our variants
all_variants <-  data.frame(matrix(ncol = 1))
colnames(all_variants) <- "Pos"


for (SRR in SRR_names) {
  SRR_pos_level <- data.frame(SRR_table_list_HET_OR_LOWLVL_nofilt[[SRR]]$Pos, SRR_table_list_HET_OR_LOWLVL_nofilt[[SRR]]$VariantLevel)
  colnames(SRR_pos_level) <- c("Pos", paste0(SRR,"_variant_lvl"))
  all_variants <- merge(all_variants, SRR_pos_level, by = "Pos", all = T)
}
all_variants$OurPos <- all_variants$Pos
all_variants_and_Ludwigs <- merge(all_variants, Ludwig_variants, by.x = "Pos", by.y = "tobecombined_Pos", all = T)
all_pos_and_Ludwigs <- data.frame(all_variants_and_Ludwigs$Pos, all_variants_and_Ludwigs$OurPos, all_variants_and_Ludwigs$Ludwig_variant_positions)


all_variants_HET_OR_LOWLVL <-  data.frame(matrix(ncol = 1))
colnames(all_variants_HET_OR_LOWLVL) <- "Pos"

for (SRR in SRR_names) {
  SRR_pos_level <- data.frame(SRR_table_list_HET_OR_LOWLVL[[SRR]]$Pos, SRR_table_list_HET_OR_LOWLVL[[SRR]]$VariantLevel)
  colnames(SRR_pos_level) <- c("Pos", paste0(SRR,"_variant_lvl"))
  all_variants_HET_OR_LOWLVL <- merge(all_variants_HET_OR_LOWLVL, SRR_pos_level, by = "Pos", all = T)
}
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




  ########################   Mutation Plots   ############################


barplot_lims <- data.frame(0:16569, rep(1,16570))
colnames(barplot_lims) <- c("Position", "ylimit")

  # get SRRs for lineage_paths.txt from lineage tree (S1d_lineage_tree.png)
Snumb_path <- list("bulk", "S0014", "S0028", "S0034", "S0049")  # add Snumbs here. S MUST BE CAPITALIZED. S000 and S0001 not recognised - use "bulk" instead.

get_SRRs_from_Snumbs <- function(Snumb_path){  # See S1d_lineage_tree.png (labelled with S#### sample names). Get list of SRRs to place in lineage_paths.txt (Don't forget to choose and add a name in front of the path list).
  Snumbs_all <- c("bulk","bulk","S0003","S0004","S0005","S0006","S0007","S0008","S0009","S0010","S0011","S0012","S0013","S0014","S0015","S0016","S0017","S0018","S0019","S0020","S0021","S0022","S0023","S0024","S0025","S0026","S0027","S0028","S0029","S0030","S0031","S0032","S0033","S0034","S0035","S0036","S0037","S0038","S0039","S0040","S0041","S0042","S0043","S0044","S0045","S0046","S0047","S0048","S0049","S0050","S0051","S0052","S0053","S0054","S0055","S0056","S0057","S0058","S0059","S0060","S0061","S0062","S0063","S0064","S0065","S0066","S0067","S0068","S0069")
  SRR_path <- list()
    i <- match(Snumb_path, Snumbs_all)
    print(typeof(i))
    for (x in i){
    print(x)
      SRR_path <- paste(SRR_path, SRR_names[x])
  }
  
  return(SRR_path)
}

SRR_path <- get_SRRs_from_Snumbs(Snumb_path)
print(SRR_path)

# add empty rows to SRR_table_list of variant information, so there is one row for every position (for x axis of mutation plots)
for (i in SRR_names){
  print(i)
  SRR_table_list[[i]]  <- merge(SRR_table_list[[i]], barplot_lims, by.x = "Pos", by.y = "Position", all = T)
}


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
    plots_in_lineage[[SRR]] <- ggplot(data = SRR_table_list[[SRR]], aes(Pos, VariantLevel)) + 
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
  file_string <- paste0("results/",p[[1]],"_nofilt.png")
  px_height <- 1.2*length(plots_in_lineage)+0.8
  ggsave(file=file_string, plot=lab_lineage_grob, width = 8, height = px_height, units = "in")
}





   ###  Repeat for SRR_table_list_HET_OR_LOWLVL_nofilt  ###

# More easy to select specific position to plot mutation load profiles using only heteroplasmic/low-level variants.
# Unfiltered as some variants seem to switch filter status between generations

# add empty rows to SRR_table_list of variant information, so there is one row for every position (for x axis of mutation plots)
for (i in SRR_names){
  print(i)
  SRR_table_list_HET_OR_LOWLVL_nofilt[[i]]  <- merge(x = SRR_table_list_HET_OR_LOWLVL_nofilt[[i]], y = barplot_lims, by.x = "Pos", by.y = "Position", all = T)
}

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
    plots_in_lineage[[SRR]] <- ggplot(data = SRR_table_list_HET_OR_LOWLVL_nofilt[[SRR]], aes(Pos, VariantLevel)) + 
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
      geom_line(data = depths_qfilt, aes(Pos, (log2(depths_qfilt[[i]]))/(log2(third_y_lim_maxcoverage_qfilt))), alpha=0.7, size = 0.15) + # coverage track
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


remove(SRR_table_list_HET_OR_LOWLVL_nofilt)








###  Repeat for SRR_table_list_HET_OR_LOWLVL  ###
# More easy to select specific position to plot mutation load profiles using only heteroplasmic/low-level variants.
# Filtered

# add empty rows to SRR_table_list of variant information, so there is one row for every position (for x axis of mutation plots)
for (i in SRR_names){
  print(i)
  SRR_table_list_HET_OR_LOWLVL[[i]]  <- merge(x = SRR_table_list_HET_OR_LOWLVL[[i]], y = barplot_lims, by.x = "Pos", by.y = "Position", all = T)
}

#print("length of SRR_table_list_HET_OR_LOWLVL:")
#print(length(SRR_table_list_HET_OR_LOWLVL))
#print("structure of SRR 80 in SRR_table_list_HET_OR_LOWLVL[[SRR 80]]")
#print(str(SRR_table_list_HET_OR_LOWLVL[["SRR7245880"]]))
#print(nrow(SRR_table_list_HET_OR_LOWLVL[["SRR7245880"]]$Pos))
#print(levels(SRR_table_list_HET_OR_LOWLVL[["SRR7245880"]]$Filter))




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








  ###############  Position-specific mutation load plots  #####################
  #####################  for variants of interest  ############################



#pos_of_interest <- 4769

# For each specified path, move over and record SRRs until VARIANTS_OF_INTEREST 
# are listed. Then for each VARIANT_OF_INTEREST position, plot line graph of 
# how variant's mutation load changes across the lineage path.

for (p in paths){
  if (p[[1]] == "#"){
    print("skipping comment line...")
    next
  } 
  if (substr(p[[1]], 0, 3)=="B11"){
    lin_col <- "palevioletred1"
  } 
  if (substr(p[[1]], 0, 2)=="B5"){
    lin_col <- "burlywood4"
  } 
  if (substr(p[[1]], 0, 2)=="F4"){
    lin_col <- "darkgrey"
  } 
  if (substr(p[[1]], 0, 2)=="D2"){
    lin_col <- "yellow1"
  } 
  if (substr(p[[1]], 0, 2)=="B3"){
    lin_col <- "orange"
  } 
  if (substr(p[[1]], 0, 3)=="G11"){
    lin_col <- "mediumorchid1"
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
# and followed by the positions of variants of interest.
# eg. 
#    LINEAGE_PATH_NAME SRR1 SRR2 SRR3 VARIANTS_OF_INTEREST 1495 12788

    if (string == "VARIANTS_OF_INTEREST") {
      print("Reached VARIANTS_OF_INTEREST for this lineage")
      at.positions <- T
      print(paste("n =",n))
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
      print(paste0("Plotting position: ", pos_of_interest, ", for lineage path: ", p[[1]]))

# make new data frame for new variant position
      mut_load_change <- data.frame(matrix(nrow = length(SRRs_in_path), ncol = 3))
      colnames(mut_load_change) <- c("SRR", "Generation", "VariantLevel")
      print(paste("SRRs_in_path: ", SRRs_in_path))
      for (SRR_name in SRRs_in_path){
        n=n+1
        print(SRR_name)
        mut_load_change$SRR[[n]] <- SRR_name
        mut_load_change$Generation[n] <- n-1
        mut_load_change$VariantLevel[n] <- SRR_table_list[[SRR_name]]$VariantLevel[pos_of_interest+1]

      }
# make new plot for new variant position
      mut_load_change[is.na(mut_load_change)] <- 0
      plot_title <- paste0(p[[1]],": ", pos_of_interest)
      mut_plot <- ggplot(data = mut_load_change, aes(x=Generation,y=VariantLevel)) +
        geom_line() +
        geom_point(aes(colour = lin_col), size = 3) +
        theme_minimal() +
        theme(plot.background = element_rect(fill = "white",
                                colour = "white"),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               axis.line = element_line(colour = "black"), legend.position = "none") +
        ggtitle(plot_title)
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
# save plot
      file_string <- paste0("results/",p[[1]],"_pos_",pos_of_interest, ".png")
      ggsave(file=file_string, plot=mut_plot)
      n=0
      print("n reset")
      } 
    else {  # if VARIANTS_OF_INTEREST hasn't been reached (at.positions=F)
      next
    }

  }
}
