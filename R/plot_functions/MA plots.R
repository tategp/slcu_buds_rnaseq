# MA plots

############
## Set up ##
############

library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)

setwd("~/Part_III_Project/")

res_N_LP <- readRDS("scripts/R/objects/res_N_LP.rds")
res_N_HP <- readRDS("scripts/R/objects/res_N_HP.rds")
res_P_LN <- readRDS("scripts/R/objects/res_P_LN.rds")
res_P_HN <- readRDS("scripts/R/objects/res_P_HN.rds")
res_int <- readRDS("scripts/R/objects/res_int.rds")

Alpha = 0.05

#############
## MA plot ##
#############

MA = function(set){
  set %>% 
    mutate(sig = ifelse(padj < Alpha, log2FoldChange, NA)) %>% 
    mutate(above = ifelse(log2FoldChange > 5, 5, NA)) %>% 
    mutate(below = ifelse(log2FoldChange < -5, -5, NA)) %>%
    mutate(sig_above = ifelse(padj < Alpha, above, NA)) %>% 
    mutate(sig_below = ifelse(padj < Alpha, below, NA)) %>%  
    ggplot(aes(baseMean, log2FoldChange)) +
    geom_point(size = 0.5) +
    geom_point(aes(y = above), shape = 2) +
    #geom_point(aes(y = below), shape = 6) +
    geom_point(aes(y = sig_above), shape = 2, colour = "red3") +
    #geom_point(aes(y = sig_below), shape = 6, colour = "red3") +
    geom_point(aes(y = sig), colour = "red3", size = 0.5) +
    geom_hline(yintercept = 0, colour = "royalblue", size = 1) +
    scale_x_log10() +
    ylim(-5,5) +
    labs(title = "Change in nitrate for high plasticity", x = "Mean expression", y = "Log2 Fold Change") +
    theme_bw() +
    theme(panel.background = element_blank()) 
      
}



# For the actual plots

# Change in nitrate for low plasticity
res_N_LP %>% 
  mutate(sig = ifelse(padj < Alpha, log2FoldChange, NA)) %>% 
  mutate(above = ifelse(log2FoldChange > 5, 5, NA)) %>% 
  mutate(below = ifelse(log2FoldChange < -5, -5, NA)) %>%
  mutate(sig_above = ifelse(padj < Alpha, above, NA)) %>% 
  mutate(sig_below = ifelse(padj < Alpha, below, NA)) %>%  
  ggplot(aes(baseMean, log2FoldChange)) +
  geom_point(size = 0.5) +
  #geom_point(aes(y = above), shape = 2) +
  geom_point(aes(y = below), shape = 6) +
  #geom_point(aes(y = sig_above), shape = 2, colour = "red3") +
  #geom_point(aes(y = sig_below), shape = 6, colour = "red3") +
  #geom_point(aes(y = sig), colour = "red3", size = 0.5) +
  geom_hline(yintercept = 0, colour = "royalblue", size = 1) +
  scale_x_log10() +
  ylim(-5,5) +
  labs(title = "A. Change in nitrate for non-plastic lines", x = "Mean expression", y = "Log2 Fold Change") +
  theme_bw() +
  theme(panel.background = element_blank())

# Change in nitrate for high plasticity
res_N_HP %>% 
  mutate(sig = ifelse(padj < Alpha, log2FoldChange, NA)) %>% 
  mutate(above = ifelse(log2FoldChange > 5, 5, NA)) %>% 
  mutate(below = ifelse(log2FoldChange < -5, -5, NA)) %>%
  mutate(sig_above = ifelse(padj < Alpha, above, NA)) %>% 
  mutate(sig_below = ifelse(padj < Alpha, below, NA)) %>%  
  ggplot(aes(baseMean, log2FoldChange)) +
  geom_point(size = 0.5) +
  geom_point(aes(y = above), shape = 2) +
  #geom_point(aes(y = below), shape = 6) +
  geom_point(aes(y = sig_above), shape = 2, colour = "red3") +
  #geom_point(aes(y = sig_below), shape = 6, colour = "red3") +
  geom_point(aes(y = sig), colour = "red3", size = 0.5) +
  geom_hline(yintercept = 0, colour = "royalblue", size = 1) +
  scale_x_log10() +
  ylim(-5,5) +
  labs(title = "B. Change in nitrate for plastic", x = "Mean expression", y = "Log2 Fold Change") +
  theme_bw() +
  theme(panel.background = element_blank())

# Different plasticity group for low nitrate
res_P_LN %>% 
  mutate(sig = ifelse(padj < Alpha, log2FoldChange, NA)) %>% 
  mutate(above = ifelse(log2FoldChange > 5, 5, NA)) %>% 
  mutate(below = ifelse(log2FoldChange < -5, -5, NA)) %>%
  mutate(sig_above = ifelse(padj < Alpha, above, NA)) %>% 
  mutate(sig_below = ifelse(padj < Alpha, below, NA)) %>% 
  ggplot(aes(baseMean, log2FoldChange)) +
  geom_point(size = 0.5) +
  geom_point(aes(y = above), shape = 2) +
  geom_point(aes(y = below), shape = 6) +
  geom_point(aes(y = sig_above), shape = 2, colour = "red3") +
  geom_point(aes(y = sig_below), shape = 6, colour = "red3") +
  geom_point(aes(y = sig), colour = "red3", size = 0.5) +
  geom_hline(yintercept = 0, colour = "royalblue", size = 1) +
  scale_x_log10() +
  ylim(-5,5) +
  labs(title = "C. Different plasticity groups at low nitrate", x = "Mean expression", y = "Log2 Fold Change") +
  theme_bw() +
  theme(panel.background = element_blank())

# Different plasticity group for high nitrate
res_P_HN %>% 
  mutate(sig = ifelse(padj < Alpha, log2FoldChange, NA)) %>% 
  mutate(above = ifelse(log2FoldChange > 5, 5, NA)) %>% 
  mutate(below = ifelse(log2FoldChange < -5, -5, NA)) %>%
  mutate(sig_above = ifelse(padj < Alpha, above, NA)) %>% 
  mutate(sig_below = ifelse(padj < Alpha, below, NA)) %>% 
  ggplot(aes(baseMean, log2FoldChange)) +
  geom_point(size = 0.5) +
  geom_point(aes(y = above), shape = 2) +
  geom_point(aes(y = below), shape = 6) +
  geom_point(aes(y = sig_above), shape = 2, colour = "red3") +
  geom_point(aes(y = sig_below), shape = 6, colour = "red3") +
  geom_point(aes(y = sig), colour = "red3", size = 0.5) +
  geom_hline(yintercept = 0, colour = "royalblue", size = 1) +
  scale_x_log10() +
  ylim(-5,5) +
  labs(title = "D. Different plasticity groups at high nitrate", x = "Mean expression", y = "Log2 Fold Change") +
  theme_bw() +
  theme(panel.background = element_blank())
