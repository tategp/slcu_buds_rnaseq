##############     MDS.R     ##############

# A program that creates MDS plots for the samples.


################
#### Set up ####
################

# Libraries required
library(tidyverse)
library(DESeq2)
library(ggrepel)

# Set working directory
setwd("//slcu.cam.ac.uk/data/TeamOL/personal_folders/gregt/bud_rnaseq/")

# Load R objects
read_rds("scripts/R/objects/vds.rds")


###################
#### Euclidean ####
###################

# Create distance matrix
dists <- dist(t(assay(vsd)))

# Create mds data frame
mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(dists))

# Plot mds
mds %>% ggplot(aes(`1`, `2`, colour = plasticity, 
           shape = genotype, alpha = nitrate)) + 
  geom_point(size = 3) +
  scale_shape_manual(values = seq(0,11)) +
  scale_color_manual(values = c("red", "blue")) +
  scale_alpha_manual(values = c(0.3, 1))


#################
#### Poisson ####
#################

# Create distance matrix
p_dists <- PoissonDistance(t(counts(dds_counts_filtered)))

# Create mds data frame
mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(p_dists$dd))

# Plot mds
mds %>% ggplot(aes(`1`, `2`, colour = plasticity, 
                   shape = genotype, alpha = nitrate)) + 
  geom_point(size = 3) +
  scale_shape_manual(values = seq(0,11)) +
  scale_color_manual(values = c("red", "blue")) +
  scale_alpha_manual(values = c(0.3, 1))
