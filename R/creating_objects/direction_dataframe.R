# Creating a directional dataframe and set of differentially expressed genes

############
## Set Up ## 
############

library(readr)
library(dplyr)

setwd("~/Part_III_Project/")

res_N_LP <- readRDS("scripts/R/objects/res_N_LP.rds")
res_N_HP <- readRDS("scripts/R/objects/res_N_HP.rds")
res_P_LN <- readRDS("scripts/R/objects/res_P_LN.rds")
res_P_HN <- readRDS("scripts/R/objects/res_P_HN.rds")

Alpha <- 0.05


# Just the direction of the lfc
N_LP <- mutate(res_N_LP, dir_N_LP = ifelse(log2FoldChange > 0, "up","down")) %>%
  select(gene = row,dir_N_LP)
N_HP <- mutate(res_N_HP, dir_N_HP = ifelse(log2FoldChange > 0, "up","down")) %>%
  select(gene = row,dir_N_HP)
P_LN <- mutate(res_P_LN, dir_P_LN = ifelse(log2FoldChange > 0, "up","down")) %>%
  select(gene = row,dir_P_LN)
P_HN <- mutate(res_P_HN, dir_P_HN = ifelse(log2FoldChange > 0, "up","down")) %>%
  select(gene = row,dir_P_HN)

direction <- N_LP %>% 
  full_join(N_HP, by = "gene") %>% 
  full_join(P_LN, by = "gene") %>% 
  full_join(P_HN, by = "gene") %>% 
  mutate(full = paste(dir_N_LP,dir_N_HP,dir_P_LN,dir_P_HN, sep = "_"))

rm(N_LP,N_HP,P_LN,P_HN)


# Direction of differentially expressed genes
sig_N_LP <- mutate(res_N_LP, dir_N_LP = ifelse(log2FoldChange > 0, "up","down")) %>%
  mutate(dir_N_LP = ifelse(padj < Alpha, dir_N_LP, "none")) %>% 
  mutate(dir_N_LP = ifelse(is.na(dir_N_LP), "none", dir_N_LP)) %>% 
  select(gene = row,dir_N_LP)

sig_N_HP <- mutate(res_N_HP, dir_N_HP = ifelse(log2FoldChange > 0, "up","down")) %>% 
  mutate(dir_N_HP = ifelse(padj < Alpha, dir_N_HP, "none")) %>% 
  mutate(dir_N_HP = ifelse(is.na(dir_N_HP), "none", dir_N_HP)) %>% 
  select(gene = row,dir_N_HP)

sig_P_LN <- mutate(res_P_LN, dir_P_LN = ifelse(log2FoldChange > 0, "up","down")) %>%
  mutate(dir_P_LN = ifelse(padj < Alpha, dir_P_LN, "none")) %>% 
  mutate(dir_P_LN = ifelse(is.na(dir_P_LN), "none", dir_P_LN)) %>% 
  select(gene = row,dir_P_LN)

sig_P_HN <- mutate(res_P_HN, dir_P_HN = ifelse(log2FoldChange > 0, "up","down")) %>%
  mutate(dir_P_HN = ifelse(padj < Alpha, dir_P_HN, "none")) %>% 
  mutate(dir_P_HN = ifelse(is.na(dir_P_HN), "none", dir_P_HN)) %>% 
  select(gene = row,dir_P_HN)

sig_direction <- sig_N_LP %>% 
  full_join(sig_N_HP, by = "gene") %>% 
  full_join(sig_P_LN, by = "gene") %>% 
  full_join(sig_P_HN, by = "gene") %>% 
  mutate(full = paste(dir_N_LP,dir_N_HP,dir_P_LN,dir_P_HN, sep = "_")) 

rm(sig_N_LP,sig_N_HP,sig_P_LN,sig_P_HN)

# Create set of genes which are differentially expressed
no_threshold_set <-  as.data.frame(sig_direction[which(sig_direction$full != "none_none_none_none"),]$gene)
colnames(no_threshold_set) <- c("row")


##################
## Save objects ##
##################

write_rds(direction, "scripts/R/objects/direction.rds")
write_rds(sig_direction, "scripts/R/objects/sig_direction.rds")
write_rds(no_threshold_set, "scripts/R/objects/no_threshold_set.rds")
