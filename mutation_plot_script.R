

  ## NOTE: this script should be called via the bash script "plot_mutations.sh" ##

# set working directory
args <- commandArgs(trailingOnly = T)
print(args)
setwd(args[1])
#setwd("/home/thomas/Documents/Research_proj/Ludwig_2019/")
  





  ## Load packages ##
# Check if packages are installed, install if needed.
local_lib_path <- paste0(args[1],"/.R_local_lib/")
print(local_lib_path)
.libPaths(c(local_lib_path, .libPaths()))

packages <- c("tidyr","ggplot2","gridExtra","ggrepel","egg","grid")
lapply(packages, FUN = function(i) {
  if (!require(i, character.only = TRUE)) {
    install.packages(i, dependencies = TRUE, lib = local_lib_path, repos="https://www.stats.bris.ac.uk/R/")
    library(i, character.only = TRUE)
  }
  }
)






  ## Read SRR files ##
filenames <- list.files("vcf/", pattern="*.txt")
# Create list of data frame names without the ".txt" part 
SRR_names <-substr(filenames,1,10)

  ## Read coverage files ##
depths <- read.table("depths.txt", sep = "\t", header = F, stringsAsFactors = T)
depths_qfilt <- read.table("depths_qfilt.txt", sep = "\t", header = F, stringsAsFactors = T)
colnames(depths_qfilt) <- c("chr", "Pos", SRR_names)
colnames(depths) <- c("chr", "Pos", SRR_names)

  ## Read lineage path/s ##
# path through a lineage specified and read from 'lineage_paths.txt'
paths <- list()
paths <- as.list(strsplit(readLines("lineage_paths.txt"), " "))








  #####################  Coverage plots and stats  ###################### (TODO)
## x axis as sample and as pos.
coverage_plots <- list()
for (i in SRR_names){
#  print(mean(depths[[i]]))
  coverage_plots[[i]] <- ggplot(data = depths_qfilt, aes(Pos, depths_qfilt[[i]])) + 
    geom_line()
}

#coverage_plots$SRR7245881

#mean_coverage_plot<- ggplot(data = depths, aes(Pos, mean(depths)) + 
#geom_line()



  # MAIN LISTS OF DATAFRAMES #

SRR_table_list <- list()  # All 
SRR_table_list_PASS <- list()  # Filtered
SRR_table_list_INTERESTING <- list()  # Filtered and only heteroplasmic/low level
SRR_table_list_INTERESTING_nofilt <- list()

#threshold <- 0.05
# Load files into list of data.frames
for(i in SRR_names){
  filepath <- file.path("vcf",paste(i,"_annotated.txt",sep=""))
  SRR_table_list[[i]] <- read.table(filepath, sep = "\t", header = T, stringsAsFactors = T)
  
  # subset for variants which passed filter
  SRR_table_list_PASS[[i]] <- subset(SRR_table_list[[i]], Filter == "PASS")
  # subset for "interesting" variants
  SRR_table_list_INTERESTING[[i]] <- subset(SRR_table_list_PASS[[i]], Type == 2)
  SRR_table_list_INTERESTING_nofilt[[i]] <- subset(SRR_table_list[[i]], Type ==2)
  }
print("Before merging structure of SRR_table_list")
print(str(SRR_table_list[["SRR7245880"]]))


  ## Variant calling stats ##

variant_stats <- data.frame(matrix(nrow = length(SRR_table_list), ncol = 8))
colnames(variant_stats) <- c("SRR","No.Variants", "No.Unfiltered_Variants", "No.het", "No.hom","No.transition", "No.transversion", "No.missense")
variant_stats$SRR <- SRR_names
variant_stats$

#for (i in SRR_names){
#  variant_stats$No.Variants[[i]] <- nrow(SRR_table_list[[i]])
#}

   ####  all_variants_in_path
all_variants_in_path
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
    
    #merge(SRR_table_list$SRR)

}
}


  ########################   Mutation Plots   ############################

barplot_lims <- data.frame(0:16569, rep(1,16570))
colnames(barplot_lims) <- c("Position", "ylimit")

# add empty rows to SRR_table_list of variant information, so there is one row for every position (for x axis of mutation plots)
for (i in SRR_names){
  print(i)
  SRR_table_list[[i]]  <- merge(SRR_table_list[[i]], barplot_lims, by.x = "Pos", by.y = "Position", all = T)
}

#print("length of SRR_table_list:")
#print(length(SRR_table_list))
#print("structure of SRR 80 in SRR_table_list[[SRR 80]]")
#print(str(SRR_table_list[["SRR7245880"]]))
#print(nrow(SRR_table_list[["SRR7245880"]]$Pos))
#print(levels(SRR_table_list[["SRR7245880"]]$Filter))




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
      geom_col(width = 1, aes(colour = factor(Filter))) + 
      scale_color_manual(values = c("PASS" = "green4",
                                    "STRAND_BIAS"="red",
                                    "BLACKLISTED"="black")) +
      geom_point(aes(colour = factor(Filter)), size = 0.8) +
      theme_minimal() + 
      ylab(SRR) +
      theme(axis.text.x = element_text(),
            axis.ticks.x = element_line(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 8),
            legend.position = "none",
            plot.margin = margin(t=0.3, r=0.1, b=0.1, l=0.1, "cm")) +
      geom_text_repel(aes(label = Pos), size = 2, nudge_y = 0.05, label.padding = 0.03, box.padding = 0.03, max.overlaps = 13)

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
  file_string <- paste0("plots/",p[[1]],".png")
  px_height <- 1.2*length(plots_in_lineage)+0.8
  ggsave(file=file_string, plot=lab_lineage_grob, width = 8, height = px_height, units = "in")
}







   ###  Repeat for SRR_table_list_INTERESTING_nofilt  ###
# More easy to select specific position to plot mutation load profiles using only heteroplasmic/low-level variants.
# Unfiltered as some variants seem to switch filter status between generations

# add empty rows to SRR_table_list of variant information, so there is one row for every position (for x axis of mutation plots)
for (i in SRR_names){
  print(i)
  SRR_table_list_INTERESTING_nofilt[[i]]  <- merge(x = SRR_table_list_INTERESTING_nofilt[[i]], y = barplot_lims, by.x = "Pos", by.y = "Position", all = T)
}

#print("length of SRR_table_list_INTERESTING_nofilt:")
#print(length(SRR_table_list_INTERESTING_nofilt))
#print("structure of SRR 80 in SRR_table_list_INTERESTING_nofilt[[SRR 80]]")
#print(str(SRR_table_list_INTERESTING_nofilt[["SRR7245880"]]))
#print(nrow(SRR_table_list_INTERESTING_nofilt[["SRR7245880"]]$Pos))
#print(levels(SRR_table_list_INTERESTING_nofilt[["SRR7245880"]]$Filter))




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
    print(nlevels(SRR_table_list_INTERESTING_nofilt[[SRR]]$Filter))
    if (nlevels(SRR_table_list_INTERESTING_nofilt[[SRR]]$Filter)==3){
      colours <- c("black", "green4", "red")
    }
    if (nlevels(SRR_table_list_INTERESTING_nofilt[[SRR]]$Filter)==2){
      colours <- c("green4", "red")
    }
    if (nlevels(SRR_table_list_INTERESTING_nofilt[[SRR]]$Filter)==1){
      colours <- c("green4")
    }
    # Create individual plot:    
    plots_in_lineage[[SRR]] <- ggplot(data = SRR_table_list_INTERESTING_nofilt[[SRR]], aes(Pos, VariantLevel)) + 
      geom_col(width = 1, aes(colour = factor(Filter))) + 
      scale_color_manual(values = c("PASS" = "green4",
                                    "STRAND_BIAS"="red",
                                "BLACKLISTED"="black")) + 
      geom_point(aes(colour = factor(Filter)), size = 0.8) +
      theme_minimal() + 
      ylab(SRR) +
      theme(axis.text.x = element_text(),
            axis.ticks.x = element_line(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 8),
            legend.position = "none",
            plot.margin = margin(t=0.3, r=0.1, b=0.1, l=0.1, "cm")) +
      geom_text_repel(aes(label = Pos), size = 2, nudge_y = 0.05, label.padding = 0.03, box.padding = 0.03, max.overlaps = 13)

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
  file_string <- paste0("plots/",p[[1]],"_INTERESTING_nofilt.png")
  px_height <- 500*length(plots_in_lineage)+370
  ggsave(file=file_string, plot=lab_lineage_grob, width = 3600, height = px_height, units = "px")
}


remove(SRR_table_list_INTERESTING_nofilt)








###  Repeat for SRR_table_list_INTERESTING  ###
# More easy to select specific position to plot mutation load profiles using only heteroplasmic/low-level variants.
# Filtered

# add empty rows to SRR_table_list of variant information, so there is one row for every position (for x axis of mutation plots)
for (i in SRR_names){
  print(i)
  SRR_table_list_INTERESTING[[i]]  <- merge(x = SRR_table_list_INTERESTING[[i]], y = barplot_lims, by.x = "Pos", by.y = "Position", all = T)
}

#print("length of SRR_table_list_INTERESTING:")
#print(length(SRR_table_list_INTERESTING))
#print("structure of SRR 80 in SRR_table_list_INTERESTING[[SRR 80]]")
#print(str(SRR_table_list_INTERESTING[["SRR7245880"]]))
#print(nrow(SRR_table_list_INTERESTING[["SRR7245880"]]$Pos))
#print(levels(SRR_table_list_INTERESTING[["SRR7245880"]]$Filter))




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
    print(nlevels(SRR_table_list_INTERESTING[[SRR]]$Filter))
    if (nlevels(SRR_table_list_INTERESTING[[SRR]]$Filter)==3){
      colours <- c("black", "green4", "red")
    }
    if (nlevels(SRR_table_list_INTERESTING[[SRR]]$Filter)==2){
      colours <- c("green4", "red")
    }
    if (nlevels(SRR_table_list_INTERESTING[[SRR]]$Filter)==1){
      colours <- c("green4")
    }
    # Create individual plot:    
    plots_in_lineage[[SRR]] <- ggplot(data = SRR_table_list_INTERESTING[[SRR]], aes(Pos, VariantLevel)) + 
      geom_col(width = 1, aes(colour = factor(Filter))) + 
      scale_color_manual(values = c("PASS" = "green4",
                                    "STRAND_BIAS"="red",
                                    "BLACKLISTED"="black")) + 
      geom_point(aes(colour = factor(Filter)), size = 0.8) +
      theme_minimal() + 
      ylab(SRR) +
      theme(axis.text.x = element_text(),
            axis.ticks.x = element_line(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 8),
            legend.position = "none",
            plot.margin = margin(t=0.3, r=0.1, b=0.1, l=0.1, "cm")) +
      geom_text_repel(aes(label = Pos), size = 2, nudge_y = 0.05, label.padding = 0.03, box.padding = 0.03, max.overlaps = 13)
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
  file_string <- paste0("plots/",p[[1]],"_INTERESTING.png")
  px_height <- 500*length(plots_in_lineage)+370
  ggsave(file=file_string, plot=lab_lineage_grob, width = 3600, height = px_height, units = "px")
}








  ###############  Position-specific mutation load plots  #####################
  #####################  for variants of interest  ############################



#pos_of_interest <- 4769

# For each specified path, move over and record SRRs until VARIANTS_OF_INTEREST 
# are listed. Then for each VARIANT_OF_INTEREST position, plot line graph of 
# how variant's mutation load changes inacross each SRR in the lineage path.

for (p in paths){
  if (p[[1]] == "#"){
    print("skipping comment line...")
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
# Add SRR to list, until the positions of variants of interest are listed
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
      mut_plot <- ggplot(data = mut_load_change, aes(x=Generation,y=VariantLevel)) +
        geom_line() +
        theme_minimal() +
        theme(panel.background = element_rect(fill = "white",
                                colour = "white"))
# save plot
      file_string <- paste0("plots/",p[[1]],"_pos_",pos_of_interest, ".png")
      ggsave(file=file_string, plot=mut_plot)
      n=0
      print("n reset")
      } else {  # for if VARIANTS_OF_INTEREST hasn't been reached.
      next
    }

  }
}
