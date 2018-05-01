# Applying vst and zscores to dds

############
## Set up ##
############

library(tidyverse)
library(DESeq2)
library(dplyr)
library(readr)
library(tidyr)

# Load objects
dds <- readRDS("scripts/R/objects/dds.rds")

# Applying normalizing
vst <- varianceStabilizingTransformation(dds, blind = FALSE) %>% 
  assay()

# Combining replicate
vst_df <- as.data.frame(vst) %>% 
  mutate(gene = rownames(.)) %>% 
  gather(id, vst, -gene) %>% 
  separate(id, c("genotype", "plasticity", "nitrate", "rep"), remove = FALSE)

# Combining genotypes
vst_per_plas <- vst_df %>% 
  group_by(gene, plasticity, nitrate) %>% 
  summarise(vst_mean = mean(vst),
            vst_range = diff(range(vst)),
            n = n()) %>% 
  ungroup()

# For counts
vst_counts <- vst_per_plas %>% 
  mutate(id = paste(plasticity, nitrate, sep = "_")) %>% 
  select(id, gene, vst_mean) %>% 
  spread(id, vst_mean) %>% 
  as.data.frame()

rownames(vst_counts) <- vst_counts$gene
vst_counts <- vst_counts[, -1]

# For z-scores
vst_zscore <- vst_per_plas %>% 
  mutate(id = paste(plasticity, nitrate, sep = "_")) %>% 
  group_by(gene) %>% 
  mutate(vst_z = (vst_mean - mean(vst_mean))/sd(vst_mean)) %>% 
  select(id, gene, vst_z) %>% 
  spread(id, vst_z) %>% 
  as.data.frame()

rownames(vst_zscore) <- vst_zscore$gene
vst_zscore <- vst_zscore[, -1]


####################
## Saving objects ##
####################

write_rds(vst, "scripts/R/objects/vst.rds")
write_rds(vst_counts, "scripts/R/objects/vst_counts.rds")
write_rds(vst_zscore, "scripts/R/objects/vst_zscore.rds")

