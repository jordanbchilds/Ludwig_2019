#args <- commandArgs(trailingOnly = T)
#print(args)
#setwd(args[1])
setwd("/home/thomas/Documents/Research_proj/Ludwig_2019/")
  

  ## Load packages ##
# add instructions to install packages if not already done.

library(tidyr)
library(ggplot2) 
library(gridExtra)
library(ggrepel)
library(egg)
library(grid)


  ## Read coverage files ##
depths <- read.table("depths.txt", sep = "\t", header = F)
depths_qfilt <- read.table("depths_mapq_20_baseq_20.txt", sep = "\t", header = F)


  ## Read SRR files ##
filenames <- list.files("vcf/", pattern="*.csv")
# Create list of data frame names without the ".csv" part 
SRR_names <-substr(filenames,1,10)
SRR_table_list <- list()


  #####################  Coverage plots and stats  ###################### (TODO)

colnames(depths) <- c("Chromosome", "Pos", SRR_names)
#av sd, min max coverage
coverage_plots <- list()
for (i in SRR_names){
  print(mean(depths[[i]]))
  coverage_plots[[i]] <- ggplot(data = depths, aes(Pos, depths[[i]])) + 
    geom_line()
}

coverage_plots$SRR7245881

#mean_coverage_plot<- ggplot(data = depths, aes(Pos, mean(depths)) + 
#geom_line()


# Load files into list of data.frames
for(i in SRR_names){
  filepath <- file.path("./vcf/",paste(i,"_annotated.csv",sep=""))
  #assign(i, read.table(filepath, sep = "\t", header = T))
  SRR_table_list[[i]] <- read.table(filepath, sep = "\t", header = T)
}


  ## Variant calling stats ##
variant_stats <- data.frame(matrix(nrow = length(SRR_table_list), ncol = 8))
colnames(variant_stats) <- c("SRR","No.Variants", "No.Unfiltered_Variants", "No.het", "No.hom","No.transition", "No.transversion", "No.missense")
variant_stats$SRR <- SRR_names

#for (i in SRR_names){
#  variant_stats$No.Variants[[i]] <- nrow(SRR_table_list[[i]])
#}


  ########################   Mutation Plots   ############################

barplot_lims <- data.frame(0:16569, rep(1,16570))
colnames(barplot_lims) <- c("Pos", "ylimit")

# add empty rows to SRR_table_list of variant information, so there is one row for every position (for x axis of mutation plots)
for (i in SRR_names){
  print(i)
  SRR_table_list[[i]]  <- merge(SRR_table_list[[i]], barplot_lims, by.x = "Pos", by.y = "Pos", all = T)
}


  ## Combine figures by lineage ##

#  read path/s of a lineage from lineage_paths.txt
paths <- list()
paths <- as.list(strsplit(readLines("lineage_paths.txt"), " "))

# For each path specified (per line in lineage_paths.txt): 
#   for each SRR in lineage path: 
#     create plot,
#   combine plots on top of each other and save to file

for (p in paths){
  print(p)
  print(str(p))
    if (p[[1]] == "#"){
    print("skipping comment line...")
    next
    }
  plots_in_lineage <- list()
  n=0
  
  for (SRR in p){
# Skip Lineage path name (1st in character vector of paths[[p]] )
    n=n+1
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
    print(nlevels(SRR_table_list[[SRR]]$Filter))
    if (nlevels(SRR_table_list[[SRR]]$Filter)==3){
      colours <- c("black", "green4", "red")
    }
    if (nlevels(SRR_table_list[[SRR]]$Filter)==2){
      colours <- c("green4", "red")
    }
    if (nlevels(SRR_table_list[[SRR]]$Filter)==1){
    print(SRR, "level of factor filter is 1??")
      }
# Create individual plot:    
    plots_in_lineage[[SRR]] <- ggplot(data = SRR_table_list[[SRR]], aes(Pos, VariantLevel)) + 
      geom_col(width = 1, aes(colour = factor(Filter))) + 
      scale_colour_manual(values = colours) + 
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
# Bottom plot (first in lineage path, second in character vector after lineage 
# name) to have x axis label.
#    if (n==2){
#      plots_in_lineage[[SRR]] <- plots_in_lineage[[SRR]] +
#        theme(axis.title.y = element_text(size = 8),
#              axis.title.x = element_text(size = 8),
#              axis.text.x = element_text(),
#              axis.ticks.x = element_line(),
#              legend.position = "none",
#              plot.margin = margin(t=0.3, r=0.1, b=0.1, l=0.1, "cm"))
#    } else {
#      plots_in_lineage[[SRR]] <- plots_in_lineage[[SRR]] +
#        theme(axis.text.x = element_text(),
#              axis.ticks.x = element_line(),
#              axis.title.x = element_blank(),
#              axis.title.y = element_text(size = 8),
#              legend.position = "none",
#              plot.margin = margin(t=0.3, r=0.1, b=0.1, l=0.1, "cm"))
#    }
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
  file_string <- paste0(p[[1]],".png")
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
  n=0
  at.positions = F
  
  for (string in p){
    n=n+1
# skip lineage name
    if (n==1){
      next
    }
# Add SRR to list, until the positions of variants of interest are listed
# instead. This is indicated by "VARIANTS_OF_INTEREST" instead of SRR name, 
# and followed by the positions of variants of interest.
# eg. 
#    LINEAGE_PATH_NAME SRR1 SRR2 SRR3 VARIANTS_OF_INTEREST 1495 12788

    SRRs_in_path <- c(SRRs_in_path, string)
    if (string == "VARIANTS_OF_INTEREST") {
      print("Reached VARIANTS_OF_INTEREST for this lineage")
      at.positions <- T
      n=0  # reset n
      next
    }
    if (at.positions == T){
      n<-n+1
      pos_of_interest <- string
      print("Plotting position:", pos_of_interest, ", for lineage path:", p[[1]])

# make new data frame for new variant position
      mut_load_change <- data.frame(matrix(nrow = length(SRRs_in_lineage), ncol = 3))
      colnames(mut_load_change) <- c("SRR", "Generation", "VariantLevel")
      
      for (SRR_name in SRRs_in_path){
        print(SRR_name, SRR_table_list[[SRR_name]]$VariantLevel[pos_of_interest+1]) # row 1 = pos 0
        mut_load_change$SRR[n] <- SRR_name
        mut_load_change$Generation[n] <- n-1
        mut_load_change$VariantLevel[n] <- SRR_table_list[[SRR_name]]$VariantLevel[pos_of_interest+1]
        n=n+1
      }
# make new plot for new variant position
      mut_load_change[is.na(mut_load_change)] <- 0
      mut_plot <- ggplot(data = mut_load_change, aes(x=Generation,y=VariantLevel)) +
        geom_line() +
        theme_minimal()
# save plot
      file_string <- paste0("plot_mutload_",p[[1]], "_",pos_of_interest,".png")
      ggsave(file=file_string, plot=mut_plot)
    
      } else {  # for if VARIANTS_OF_INTEREST hasn't been reached.
      next
    }

  }
}









  n<-n+1
  print(i)
  print(SRR_table_list[[i]]$VariantLevel[pos_of_interest+1])
  mut_load_change$SRR[n] <- i
  mut_load_change$Generation[n] <- n-1
  mut_load_change$VariantLevel[n] <- SRR_table_list[[i]]$VariantLevel[pos_of_interest+1]

mut_load_change[is.na(mut_load_change)] <- 0
plot <- ggplot(data = mut_load_change, aes(x=Generation,y=VariantLevel)) +
  geom_line() +
  theme_minimal()
  

plot



