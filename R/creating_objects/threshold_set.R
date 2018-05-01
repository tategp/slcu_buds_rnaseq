# Wald - On alternative design and with a threshold

############
## Set up ##
############

library(reshape2)
library(DESeq2)
library(readr)
library(dplyr)

setwd("~/Part_III_Project/")

dds_alt_wald <- readRDS("scripts/R/objects/dds_alt_wald.rds")
dds_ind_wald <- readRDS("scripts/R/objects/dds_ind_wald.rds")

Alpha <- 0.05
threshold <- 0.585  # Equivalent to a change of ~1.5


#############
## Results ##
#############

threshold_N_LP <- results(dds_ind_wald, name = "plasticitylowplast.nitrateHN" ,alpha = Alpha, lfcThreshold = threshold, tidy = TRUE)
threshold_P_LN <- results(dds_alt_wald, contrast = list("plasticityhighplast.nitrateLN", "plasticitylowplast.nitrateLN") ,alpha = Alpha, lfcThreshold = threshold, tidy = TRUE)
threshold_N_HP <- results(dds_ind_wald, name = "plasticityhighplast.nitrateHN" ,alpha = Alpha, lfcThreshold = threshold, tidy = TRUE)
threshold_P_HN <- results(dds_alt_wald, contrast = list("plasticityhighplast.nitrateHN", "plasticitylowplast.nitrateHN") ,alpha = Alpha, lfcThreshold = threshold, tidy = TRUE)

#######################
## Classifying genes ##
#######################

sig_N_LP <- mutate(threshold_N_LP, dir_N_LP = ifelse(log2FoldChange > 0, "up","down")) %>%
  mutate(dir_N_LP = ifelse(padj < Alpha, dir_N_LP, "none")) %>% 
  mutate(dir_N_LP = ifelse(is.na(dir_N_LP), "none", dir_N_LP)) %>% 
  select(gene = row,dir_N_LP)

sig_N_HP <- mutate(threshold_N_HP, dir_N_HP = ifelse(log2FoldChange > 0, "up","down")) %>% 
  mutate(dir_N_HP = ifelse(padj < Alpha, dir_N_HP, "none")) %>% 
  mutate(dir_N_HP = ifelse(is.na(dir_N_HP), "none", dir_N_HP)) %>% 
  select(gene = row,dir_N_HP)

sig_P_LN <- mutate(threshold_P_LN, dir_P_LN = ifelse(log2FoldChange > 0, "up","down")) %>%
  mutate(dir_P_LN = ifelse(padj < Alpha, dir_P_LN, "none")) %>% 
  mutate(dir_P_LN = ifelse(is.na(dir_P_LN), "none", dir_P_LN)) %>% 
  select(gene = row,dir_P_LN)

sig_P_HN <- mutate(threshold_P_HN, dir_P_HN = ifelse(log2FoldChange > 0, "up","down")) %>%
  mutate(dir_P_HN = ifelse(padj < Alpha, dir_P_HN, "none")) %>% 
  mutate(dir_P_HN = ifelse(is.na(dir_P_HN), "none", dir_P_HN)) %>% 
  select(gene = row,dir_P_HN)

threshold_sig_direction <- sig_N_LP %>% 
  full_join(sig_N_HP, by = "gene") %>% 
  full_join(sig_P_LN, by = "gene") %>% 
  full_join(sig_P_HN, by = "gene") %>% 
  mutate(full = paste(dir_N_LP,dir_N_HP,dir_P_LN,dir_P_HN, sep = "_")) 

remove(sig_N_LP,sig_N_HP,sig_P_LN,sig_P_HN)

# Create set of genes that are differentially expressed with the threshold
threshold_set <- as.data.frame(threshold_sig_direction[which(threshold_sig_direction$full != "none_none_none_none"),]$gene)
colnames(threshold_set) <- c("row")



####################
## Saving objects ##
####################

write_rds(threshold_sig_direction, "scripts/R/objects/threshold_sig_direction.rdS")
write_rds(threshold_set, "scripts/R/objects/threshold_set.rdS")
