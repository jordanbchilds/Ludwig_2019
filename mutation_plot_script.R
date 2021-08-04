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
#library(cowplot)
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


  ## Coverage plots and stats ##
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

str(SRR_table_list$SRR7245880)
nrow(SRR_table_list$SRR7245880)

  ## Variant calling stats ##
variant_stats <- data.frame(matrix(nrow = length(SRR_table_list), ncol = 8))
colnames(variant_stats) <- c("SRR","No.Variants", "No.Unfiltered_Variants", "No.het", "No.hom","No.transition", "No.transversion", "No.missense")
variant_stats$SRR <- SRR_names

#for (i in SRR_names){
#  variant_stats$No.Variants[[i]] <- nrow(SRR_table_list[[i]])
#}



  ## Mutation Plots ##

barplot_lims <- data.frame(0:16569, rep(1,16570))
colnames(barplot_lims) <- c("Pos", "ylimit")


plots <- list()

for (i in SRR_names){
  #print(colnames(SRR_table_list[[i]]))
  print(i)
  SRR_table_list[[i]]  <- merge(SRR_table_list[[i]], barplot_lims, by.x = "Pos", by.y = "Pos", all = T)
}


#Lineage B3
#bulk
#SRR7245888
#SRR7245905

?#SRRs_in_lineage <- c(
  #"SRR7245881",
  #"SRR7245887",
  #"SRR7245909",
  #"SRR7245915",
  #"SRR7245929",
  #"SRR7245937",
  #"SRR7245942",
  #"SRR7245944",
  #"SRR7245945")
  


  ## Combine figures by lineage ##
paths <- list()
paths <- as.list(strsplit(readLines("lineage_paths.txt"), " "))
  
#plots_in_lineage <- list()
for (p in paths){
  print(p)
  print(str(p))
    if (p[[1]] == "#"){
    print("skipping comment line...")
    next
    }
  plots_in_lineage <- list()
  n=0
  for (i in p){
    # Skip Lineage path name (1st in character vector of paths[[p]])
    n=n+1
    if (n==1){
      next
    }
    # colour according to filter: PASS, STRAND_BIAS, or BLACKLISTED
    if (nlevels(SRR_table_list[[i]]$Filter)==3){
      colours <- c("black", "green4", "red")
    }
    if (nlevels(SRR_table_list[[i]]$Filter)==2){
      colours <- c("green4", "red")
    }
    plots_in_lineage[[i]] <- ggplot(data = SRR_table_list[[i]], aes(Pos, VariantLevel)) + 
      geom_col(aes(colour = factor(Filter)), width = 1) + 
      scale_colour_manual(values = colours) + 
      geom_point(aes(colour = factor(Filter)), size = 0.8) +
      theme_minimal() + 
      ylab(i) +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 8),
            legend.position = "none",
            plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")) +
      geom_text_repel(aes(label = Pos), size = 2.5, nudge_y = 0.05, label.padding = 0.03, box.padding = 0.03, max.overlaps = 13)
    }
  lineage_plot <- ggarrange(plots = plots_in_lineage, nrow = length(plots_in_lineage), align = "hv")

   # Add x and y labels to grid of plots
  y.grob <- textGrob("Variant Level", 
                       gp=gpar(col="black", fontsize=12), rot=90)
  x.grob <- textGrob("Position in mitochondrial genome", 
                     gp=gpar(col="black", fontsize=12))

  lab_lineage_grob <- arrangeGrob(lineage_plot, left = y.grob, bottom = x.grob)
  file_string <- paste0(p[[1]],".png")
  ggsave(file=file_string, lab_lineage_grob)
}

  











for (i in SRR_names){
  #print(colnames(SRR_table_list[[i]]))
  print(i)
  if (nlevels(SRR_table_list[[i]]$Filter)==3){
    colours <- c("black", "green4", "red")
  }
  if (nlevels(SRR_table_list[[i]]$Filter)==2){
    colours <- c("green4", "red")
  }
  plots[[i]] <- ggplot(data = SRR_table_list[[i]], aes(Pos, VariantLevel)) + 
    geom_col(aes(colour = factor(Filter)), width = 1) + 
    scale_colour_manual(values = colours) + 
    geom_point(aes(colour = factor(Filter)), size = 0.8) +
    theme_minimal() + 
    ylab(i) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 8),
          legend.position = "none",
          plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")) +
    geom_text_repel(aes(label = Pos), size = 2.5, nudge_y = 0.05, label.padding = 0.03, box.padding = 0.03, max.overlaps = 13)
}











  
  ## position-specific mutation load plots
mut_load_change <- data.frame(matrix(nrow = length(SRRs_in_lineage), ncol = 3))
colnames(mut_load_change) <- c("SRR", "Generation", "VariantLevel")

pos_of_interest <- 4769
n=0
for (i in SRRs_in_lineage){
  n<-n+1
  print(i)
  print(SRR_table_list[[i]]$VariantLevel[pos_of_interest+1])
  mut_load_change$SRR[n] <- i
  mut_load_change$Generation[n] <- n-1
  mut_load_change$VariantLevel[n] <- SRR_table_list[[i]]$VariantLevel[pos_of_interest+1]
}
mut_load_change[is.na(mut_load_change)] <- 0
plot <- ggplot(data = mut_load_change, aes(x=Generation,y=VariantLevel)) +
  geom_line() +
  theme_minimal()
  

plot



