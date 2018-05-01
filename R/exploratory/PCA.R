##############     PCA.R     ##############

# A program that creates PCA plots for the samples.

################
#### Set up ####
################

# Libraries required
library(tidyverse)
library(DESeq2)
library(ggrepel)
library(RColorBrewer)

# Set working directory
setwd("~/Part_III_Project/")
setwd("//slcu.cam.ac.uk/data/TeamOL/personal_folders/gregt/bud_rnaseq/")

# Load R objects
dds_full <- readRDS("scripts/R/objects/dds_full.rds")

# Variance stabilisation 
vst_blind <- varianceStabilizingTransformation(dds_full, blind = TRUE)


##############################################
######## Adaption of Hugo's PCA code #########
##############################################

pcaDESeq <- function(object, ...){
  
  # Get the data for all genes
  selected_expr <- assay(object)
  
  # perform a PCA on the data in assay(x) for the selected genes
  ## Need to transpose the matrix as prcomp clusters by rows
  pca <- prcomp(t(selected_expr), ...)
  
  #### Eigen vectors table ####
  # Get sample information from DESeq object
  # and bind the PCA vectors
  eigen_vectors <- colData(object) %>% 
    as.data.frame() %>% 
    bind_cols(as.data.frame(pca$x))
  
  #### Eigen values table ####
  eigen_values <- data.frame(PC = colnames(pca$x), stdev = pca$sdev) %>% 
    mutate(var = stdev^2, var_pct = var/sum(var), cum_var = cumsum(var_pct))
  
  #### Factor loadings table ####
  factor_loadings <- pca$rotation %>% 
    as.data.frame() %>% 
    mutate(variable = row.names(.)) %>% 
    dplyr::select(variable, everything())
  
  #### Convert the original data to a data.frame ####
  selected_expr <- selected_expr %>% 
    as.data.frame() %>% 
    rename_all(funs(paste0("sample", .))) %>% 
    mutate(variable = rownames(.)) %>% 
    dplyr::select(variable, everything())
  
  # Return a list with each of these objects
  return(list(vectors = eigen_vectors, 
              values = eigen_values, 
              loadings = factor_loadings,
              original = selected_expr))
  
}


##############################################
#### Applying the program to our data set ####
##############################################

pca <- pcaDESeq(vst_blind)

# PCA plot
pca$vectors %>% 
  ggplot(aes(PC1, PC2, colour = interaction(plasticity,nitrate), 
             shape = genotype)) + 
  geom_point(size = 5) +
  scale_colour_manual(values=c("deepskyblue", "magenta", "deepskyblue4", "red"), 
                      name  ="Group", 
                      labels=c("Low Nitrate, Non-Plastic", "Low Nitrate, Plastic", 
                               "High Nitrate, Non-Plastic", "High Nitrate, Plastic")) +
  scale_shape_manual(values = seq(0,11), name = "Genotype") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  xlim(-50,50) +
  ylim(-50,50) 



###############
# Save object #
############### 

write_rds(vst_blind,"scripts/R/objects/vst_blind.rds")

