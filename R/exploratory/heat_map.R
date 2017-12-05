##############     heat_map.R     ##############

# A program that creates heat maps for the samples.


################
#### Set up ####
################

# Libraries required
library(tidyverse)
library(DESeq2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(PoiClaClu)

# Set working directory
setwd("//slcu.cam.ac.uk/data/TeamOL/personal_folders/gregt/bud_rnaseq/")

# Load R objects
read_rds("scripts/R/objects/vds.rds")


######################################
#### Heat map for sample distance ####
######################################

# Euclidean distances
dists <- dist(t(assay(vsd_counts)))

dists_matrix <- as.matrix(dists)
rownames(dists_matrix) <- paste(vsd_counts$genotype, vsd_counts$plasticity, vsd_counts$nitrate, vsd_counts$rep, sep = "_")
colnames(dists_matrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

# Plot heat map
pheatmap(dists_matrix,
         clustering_distance_rows = dists,
         clustering_distance_cols = dists,
         col = colors)


# Poisson distances
p_dists <- PoissonDistance(t(counts(dds_counts_filtered)))

p_dists_matrix <- as.matrix(p_dists$dd)
rownames(p_dists_matrix) <- paste(vsd_counts$genotype, vsd_counts$plasticity, vsd_counts$nitrate, vsd_counts$rep, sep = "_")
colnames(p_dists_matrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

# Plot heat map
pheatmap(dists_matrix,
         clustering_distance_rows = p_dists$dd,
         clustering_distance_cols = p_dists$dd,
         col = colors)


###################################
#### Heat map for count matrix ####
###################################

select <- order(rowMeans(counts(dds_counts_filtered,normalized=FALSE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds_counts_filtered)[,c("nitrate","plasticity")])
pheatmap(assay(vsd_counts)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

