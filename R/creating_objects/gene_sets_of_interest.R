# Gene sets of interest

library(readr)
library(tidyverse)

setwd("~/Part_III_Project/")

# Parameters
clus <- 1

# Load objects
gene_clus <- readRDS("scripts/R/objects/gene_clus.rds")

# Cluster sets
clus_genes <- split(gene_clus, f = as.factor(as.numeric(gene_clus$clus)))


# Significant change set
sig_sd_change <- tibble::rownames_to_column(x,"gene") %>% 
  left_join(gene_clus) %>% 
  mutate(N_LP = ifelse(abs(lowplast_HN-lowplast_LN) < 1, 0, 1)) %>% 
  mutate(N_LP = ifelse(lowplast_HN-lowplast_LN < -1, -1, N_LP)) %>% 
  mutate(N_HP = ifelse(abs(highplast_HN-highplast_LN) < 1, 0, 1)) %>% 
  mutate(N_HP = ifelse(highplast_HN-highplast_LN < -1, -1, N_HP))%>% 
  mutate(P_LN = ifelse(abs(highplast_LN-lowplast_LN) < 1, 0, 1)) %>% 
  mutate(P_LN = ifelse(highplast_LN-lowplast_LN < -1, -1, P_LN))%>% 
  mutate(P_HN = ifelse(abs(highplast_HN-lowplast_HN) < 1, 0, 1)) %>% 
  mutate(P_HN = ifelse(highplast_HN-lowplast_HN < -1, -1, P_HN)) %>% 
  select(gene,N_LP,N_HP,P_LN,P_HN)
