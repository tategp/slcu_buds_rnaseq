##############     clustering_samples.R     ##############

# A program that applies hierachical clustering to the Euclidean distance of the vst counts.

################
#### Set up ####
################

# Libraries required
library(tidyverse)
library(DESeq2)
library(ggrepel)

# Set working directory
setwd("//slcu.cam.ac.uk/data/TeamOL/personal_folders/gregt/bud_rnaseq/")
setwd("~/Part_III_Project/")

# Load R objects
vst_blind <- read_rds("scripts/R/objects/vst_blind.rds")


#######################################
#### Clustering of sample distance ####
#######################################

# Euclidean distances
dists <- dist(t(assay(vst_blind)))

# Apply clustering
cluster <- hclust(dists)




# use ggdendro package to extract "data" from hclust
dend_data <- ggdendro::dendro_data(cluster)

# Add sample information  
dend_data$labels <- dend_data$labels %>% 
  separate(label, c("genotype", "plasticity", "nitrate", "rep"), 
            sep = "_", remove = FALSE) %>% 
  mutate(colour = ifelse(nitrate == "LN" & plasticity == "lowplast","deepskyblue",NA)) %>% 
  mutate(colour = ifelse(nitrate == "LN" & plasticity == "highplast","magenta",colour)) %>% 
  mutate(colour = ifelse(nitrate == "HN" & plasticity == "lowplast","deepskyblue4",colour)) %>% 
  mutate(colour = ifelse(nitrate == "HN" & plasticity == "highplast","red",colour)) %>% 
  mutate(y = -4) 



# Make a plot
ggplot() + 
  geom_segment(data = dend_data$segments, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_point(data = dend_data$labels, aes(x = x, y = y, colour = colour, shape = genotype), size = 2.5) +
  scale_shape_manual(values = seq(0,11), name = "Genotype") +
  scale_colour_manual(values=c("deepskyblue", "deepskyblue4", "magenta", "red")) +
  coord_flip() + 
  scale_y_reverse() +
  theme_bw() +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        legend.position="none",
        axis.title.y=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank()) +
  labs(y = "Height")
  
