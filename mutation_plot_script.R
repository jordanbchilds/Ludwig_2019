  ## Load packages ##
# add instructions to install packages if not already done.

library(tidyr)
library(ggplot2) 
library(gridExtra)
library(ggrepel)

  ## Read coverage files ##
depths <- read.table("/home/thomas/Documents/Research_proj/depths.txt", sep = "\t", header = F)
depths_qfilt <- read.table("/home/thomas/Documents/Research_proj/depths_mapq_20_baseq_20.txt", sep = "\t", header = F)


  ## Read SRR files ##
#setwd("/home/thomas/Documents/Research_proj/vcf/")
filenames <- list.files(path = "/home/thomas/Documents/Research_proj/vcf/", pattern="*.csv")
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
  filepath <- file.path("/home/thomas/Documents/Research_proj/vcf/",paste(i,"_annotated.csv",sep=""))
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
  #SRR_table_list[[i]]$Text_Pos <- SRR_table_list[[i]]$VariantLevel + 0.05
  # number levels for filter
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
    geom_text_repel(aes(label = Pos), size = 2.5, nudge_y = 0.03, label.padding = 0.03, box.padding = 0.03)#, angle = "45") #min.segment.length = 0,
}


plots$SRR7245880
#
#
#  # Check levels of nlevels of Filter #
#lvls <- list()
#for (i in SRR_names){
#  lvls <- c(lvls, nlevels(SRR_table_list[[i]]$Filter))
#}
#


#Lineage B3
#bulk
#SRR7245888
#SRR7245905

  ## Combine figures by lineage ##

SRRs_in_lineage <- c(
"SRR7245880",
"SRR7245887",
"SRR7245909",
"SRR7245915",
"SRR7245929",
"SRR7245937",
"SRR7245942",
"SRR7245944",
"SRR7245945")

n_samples <- length(SRRs_in_lineage)
fig_pos <- list()
n=0
for (i in SRRs_in_lineage){
  ystart <- n*1/n_samples
  yend <- n*2/n_samples
  fig_param <- c(0,1,ystart,yend)
  par(fig=fig_param, new=T)
  plots[[i]]
  n <- n+1
}
  
par(mfrow = c(1,length(SRRs_in_lineage)))
for (i in SRRs_in_lineage){
  plots[[i]]
}



par(fig=c(0,1,0,1), new=T)
plots$SRR7245880


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
  geom_line()
  

plot



