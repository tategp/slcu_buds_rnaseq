##############     PCA.R     ##############

# A program that creates PCA plots for the samples.


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
read_rds("scripts/R/objects/dds_filtered.rds")

# Variance stabilisation 
vsd <- varianceStabilizingTransformation(dds_counts_filtered, blind = FALSE)

write_rds(vds, "scripts/R/objects")


##############################################
######## Hugo's code for running PCA #########
##############################################

pcaDESeq <- function(object, ntop = 500, ...){
  
  # calculate the variance for each gene
  rv <- genefilter::rowVars(assay(object))
  
  # select the ntop genes by variance
  selected_genes <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  
  # Get the data for those genes
  selected_expr <- assay(object)[selected_genes,]
  
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
    select(variable, everything())
  
  #### Convert the original data to a data.frame ####
  selected_expr <- selected_expr %>% 
    as.data.frame() %>% 
    rename_all(funs(paste0("sample", .))) %>% 
    mutate(variable = rownames(.)) %>% 
    select(variable, everything())
  
  # Return a list with each of these objects
  return(list(vectors = eigen_vectors, 
              values = eigen_values, 
              loadings = factor_loadings,
              original = selected_expr))
  
}


##############################################
#### Applying the program to our data set ####
##############################################

pca <- pcaDESeq(vds)

# Cumulative plot
position <- pca$values$PC
pca$values %>% 
  ggplot(aes(PC, cum_var, group = 1)) +
  geom_line() +
  scale_x_discrete(limits = position)


# PCA plot
pca$vectors %>% 
  ggplot(aes(PC1, PC2, colour = nitrate, 
             shape = genotype, alpha = plasticity)) + 
  geom_point(size = 3) +
  scale_shape_manual(values = seq(0,11)) +
  scale_color_manual(values = c("red", "blue")) +
  scale_alpha_manual(values = c(0.3, 1))


# Plotting the main contributing genes
pca$loadings %>% 
  filter(PC1^2 > quantile(PC1^2, 0.98) | PC2^2 > quantile(PC2^2, 0.98)) %>% 
  ggplot(aes(PC1, PC2)) +
  geom_hline(yintercept = 0, colour = "grey") + geom_vline(xintercept = 0, colour = "grey") +
  geom_segment(data = pca$loadings, aes(xend = 0, yend = 0),
               arrow = arrow(ends = "first", length = unit(0.2, "cm")),
               colour = "grey") +
  geom_segment(aes(xend = 0, yend = 0),
               arrow = arrow(ends = "first", length = unit(0.2, "cm")),
               colour = "red3") +
  geom_text_repel(aes(label = variable)) +
  theme_classic()
