# Plotting gene counts for singular genes

library(DESeq2)
library(ggplot2)

dds <- readRDS("scripts/R/objects/dds.rds")

# Function to plot the normalized gene counts given a gene identifier
plot_counts <- function(candidate_gene = "AT3G52680") {
    
  # Give nitrate levels
  nitratelvl <- c("Low","High")
  
  # Retreive the count data
  geneCounts <- plotCounts(dds, gene = candidate_gene, 
                             intgroup = c("nitrate","plasticity","genotype"),
                             returnData = TRUE) %>% 
    mutate(nitrate = ifelse(nitrate == "LN", "Low", "High"))
  
  
  # Plot the counts
  ggplot(geneCounts, aes(x = nitrate, y = count, colour = plasticity, shape = genotype)) +
    geom_point(position = position_dodge(0.1), size = 2) +
    geom_pointrange(stat = "summary", fun.data = "mean_cl_boot", aes(group = plasticity),
                    position = position_dodge(0.1)) +
    geom_line(stat = "summary", fun.y = "mean", aes(group = plasticity),
                    size = 1, position = position_dodge(0.1)) +
    scale_x_discrete(limits = nitratelvl) +
    #Aesthetics...
    labs(title = candidate_gene, x = "Nitrate", y = "Counts") +
    scale_colour_manual(values=c("deepskyblue4", "red"), 
                        name  ="Plasticity", 
                        labels=c("Non-Plastic", "Plastic"))+
    scale_shape_manual(values = seq(0,11), name = "Genotype") +
    theme_bw() +
    theme(panel.background = element_blank()) 

}


#### For motif plot in presentation ###

plot_counts_motif <- function(candidate_gene = "AT5G46950") {
  
  # Give nitrate levels
  nitratelvl <- c("LN","HN")
  
  # Retreive the count data
  geneCounts <- plotCounts(dds, gene = candidate_gene, 
                           intgroup = c("nitrate","plasticity","genotype"),
                           returnData = TRUE)
  
  # Plot the counts
  ggplot(geneCounts, aes(x = nitrate, y = count, colour = plasticity, shape = genotype)) +
    #geom_point(position = position_dodge(0.1), size = 2) +
    geom_point(stat = "summary", fun.data = "mean_cl_boot", aes(group = plasticity, size = 1.5)) +
    geom_line(stat = "summary", fun.y = "mean", aes(group = plasticity),
              size = 1.5) +
    scale_x_discrete(limits = nitratelvl) +
    #Aesthetics...
    labs(title = candidate_gene, x = "Nitrate", y = "Counts") +
    scale_colour_manual(values=c("deepskyblue4", "red"), 
                        name  ="Plasticity", 
                        labels=c("Low Plasticity", "High Plasticity"))+
    scale_shape_manual(values = seq(0,11), name = "Genotype") +
    theme_bw() +
    theme(panel.background = element_blank()) 
  
}
