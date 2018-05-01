# Plotting clusters

# Set up

library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cluster)

setwd("~/Part_III_Project/")

############################################
# Function to plot the motifs with zscores #
############################################

motifs_zscore <- function(gene_set = "no_threshold_set", metric = c("euclidean","correlation"),  clustering = c("hierarchical","bruteforce"), num_clus = 8 ) {
  
  ############## 
  # Load files #
  ##############
  
  split_set <- strsplit(gene_set, "_")
  split_set <- split_set[[1]]
  
  if (split_set[1] == "go") {
    gene_set <- read.table(paste0("scripts/R/objects/",gene_set,".txt"), header = TRUE, sep = "\t")
    gene_set <- as.data.frame(gene_set[,1])
    colnames(gene_set) <- c("row")
  } else {
    gene_set <- readRDS(paste0("scripts/R/objects/",gene_set,".rds"))
    colnames(gene_set) <- c("row")
  }
    
  x <- readRDS("scripts/R/objects/vst_zscore.rds")
  x <- x[rownames(x) %in% gene_set$row,]
  
  
  # Separate case when you want to do brute force clustering
  
  if (clustering == "bruteforce") {
    
    brute_cluster <- readRDS("scripts/R/objects/brute_cluster.rds")
    
    gene_clus <- left_join(gene_set, brute_cluster, by = "row") %>% 
      mutate(clus = as.numeric(as.factor(cluster))) %>% 
      select(row,clus)
    colnames(gene_clus) <- c("gene","clus")
    
  } else {
    
    ########################
    # Dissimilarity matrix # 
    ########################
    
    if (metric == "euclidean") {
      clustering_matrix <- dist(x)
    } else if (metric == "correlation") {
      # Calculate correlation
      x_cor <- cor(t(x))
      # Express as (1 - correlation) such that
      ## r = 1 becomes 0 (very similar)
      ## r = 0 becomes 1 (quite dissimilar)
      ## r = -1 becomes 2 (very dissimilar)
      x_dis <- 1 - x_cor
      # Return as a distance object
      clustering_matrix <- as.dist(x_dis)
      rm(x_cor,x_dis)
    } else {
      stop()
    }
    
    if (clustering == "hierarchical") {
    
      # Get the clustering
      cluster <- hclust(clustering_matrix)
    
      # Cut the tree and lable genes
      gene_clus <- as.data.frame(cutree(cluster, num_clus))
      gene_clus <- as.data.frame(cbind(row.names(gene_clus),gene_clus[,1]))
      colnames(gene_clus) <- c("gene","clus")
    
    }  else {
      stop()
    }

  }
  
  #################
  # Plot clusters #
  #################

  # Get the zscore data frame in the right form 
  y <- tibble::rownames_to_column(x,"gene") %>% 
    left_join(gene_clus) %>% 
    gather(nitplas,zscore,-gene,-clus) %>% 
    mutate(nitrate = ifelse(nitplas == "lowplast_HN" | nitplas == "highplast_HN", "High","Low")) %>% 
    mutate(plasticity = ifelse(nitplas == "lowplast_LN" | nitplas == "lowplast_HN","LowPlas","HighPlas"))
  
  y$nitrate <- factor(y$nitrate, levels = c("Low","High"),ordered = TRUE)
  y$plasticity <- as.factor(y$plasticity)
  
  m <- y %>% 
    group_by(clus,nitrate, plasticity) %>% 
    summarise(mz = mean(zscore))
  
  label_text <- gene_clus %>% 
    group_by(clus) %>% 
    summarise(total = n(), nitrate = n()/n(), z = n()) %>% 
    mutate(total = paste(total,"genes",sep = " "))
  
  # Plot
  ggplot(y, aes(x = nitrate, y = zscore, 
                group = interaction(gene,plasticity) , 
                colour = plasticity)) +
    geom_point() +
    geom_line() +
    geom_point(data = m, aes(x = nitrate, y = mz, group = interaction(clus,plasticity)), colour = "black") +
    geom_line(data = m, aes(x = nitrate, y = mz, group = interaction(clus,plasticity)), colour = "black") +
    #geom_text(data = label_text, aes(x = nitrate, y = z, label = total)) +
    lims(y = c(-1.6,1.6)) +
    facet_wrap( ~ as.factor(as.numeric(clus)), ncol = 3) +
    #Aesthetics...
    labs(x = "Nitrate", y = "Z-score") +
    scale_colour_manual(values=c("tomato", "royalblue1"), 
                        name  ="Plasticity", 
                        labels=c("High Plasticity", "Low Plasticity"))+
    theme_bw() +
    theme(panel.background = element_blank())  
  
}






##################################
# Volin plot for the interaction #
##################################

motifs_violin <- function(gene_set = "threshold_set", metric = c("euclidean","correlation"),  clustering = c("hierarchical","bruteforce"), num_clus = 5 ) {
  
  ############## 
  # Load files #
  ##############
  
  split_set <- strsplit(gene_set, "_")
  split_set <- split_set[[1]]
  
  if (split_set[1] == "go") {
    gene_set <- read.table(paste0("scripts/R/objects/",gene_set,".txt"), header = TRUE, sep = "\t")
    gene_set <- as.data.frame(gene_set[,1])
    colnames(gene_set) <- c("row")
  } else {
    gene_set <- readRDS(paste0("scripts/R/objects/",gene_set,".rds"))
    colnames(gene_set) <- c("row")
  }
  
  x <- readRDS("scripts/R/objects/vst_zscore.rds")
  x <- x[rownames(x) %in% gene_set$row,]
  
  
  # Separate case when you want to do brute force clustering
  
  if (clustering == "bruteforce") {
    
    brute_cluster <- readRDS("scripts/R/objects/brute_cluster.rds")
    
    gene_clus <- left_join(gene_set, brute_cluster, by = "row") %>% 
      mutate(clus = as.numeric(as.factor(cluster))) %>% 
      select(row,clus)
    colnames(gene_clus) <- c("gene","clus")
    
  } else {
    
    ########################
    # Dissimilarity matrix # 
    ########################
    
    if (metric == "euclidean") {
      clustering_matrix <- dist(x)
    } else if (metric == "correlation") {
      # Calculate correlation
      x_cor <- cor(t(x))
      # Express as (1 - correlation) such that
      ## r = 1 becomes 0 (very similar)
      ## r = 0 becomes 1 (quite dissimilar)
      ## r = -1 becomes 2 (very dissimilar)
      x_dis <- 1 - x_cor
      # Return as a distance object
      clustering_matrix <- as.dist(x_dis)
      rm(x_cor,x_dis)
    } else {
      stop()
    }
    
    if (clustering == "hierarchical") {
      
      # Get the clustering
      cluster <- hclust(clustering_matrix)
      
      # Cut the tree and lable genes
      gene_clus <- as.data.frame(cutree(cluster, num_clus))
      gene_clus <- as.data.frame(cbind(row.names(gene_clus),gene_clus[,1]))
      colnames(gene_clus) <- c("gene","clus")
      
    }  else {
      stop()
    }
    
  }
  
  
  ###############
  # Violin plot #
  ###############
  
  z <- tibble::rownames_to_column(x,"gene") %>% 
    left_join(gene_clus) %>%
    mutate(int = (highplast_HN - highplast_LN) - (lowplast_HN - lowplast_LN))
  
  ggplot(z, aes(x = as.factor(as.numeric(clus)), y = int)) +
    geom_violin() +
    labs(x = "Cluster", y = "Interaction")
  
}



###################################################
# For K-means to determine the number of clusters #
###################################################

# Plotting the elbow graph
elbow_plot <- function(gene_set = "threshold_set", metric = c("euclidean","correlation"), nclus = 25 ) {
  
  ############## 
  # Load files #
  ##############
  
  split_set <- strsplit(gene_set, "_")
  split_set <- split_set[[1]]
  
  if (split_set[1] == "go") {
    gene_set <- read.table(paste0("scripts/R/objects/",gene_set,".txt"), header = TRUE, sep = "\t")
    gene_set <- as.data.frame(gene_set[,1])
    colnames(gene_set) <- c("row")
  } else {
    gene_set <- readRDS(paste0("scripts/R/objects/",gene_set,".rds"))
    colnames(gene_set) <- c("row")
  }
  
  x <- readRDS("scripts/R/objects/vst_zscore.rds")
  x <- x[rownames(x) %in% gene_set$row,]
  
  # Clustering matrix
  
  if (metric == "euclidean") {
    clustering_matrix <- dist(x)
  } else if (metric == "correlation") {
    # Calculate correlation
    x_cor <- cor(t(x))
    # Express as (1 - correlation) such that
    ## r = 1 becomes 0 (very similar)
    ## r = 0 becomes 1 (quite dissimilar)
    ## r = -1 becomes 2 (very dissimilar)
    x_dis <- 1 - x_cor
    # Return as a distance object
    clustering_matrix <- as.dist(x_dis)
    rm(x_cor,x_dis)
  } else {
    stop()
  }
  
  # K-means clustering
  elbow <- matrix(nrow = (nclus - 1), ncol = 2)
  
  for (k in 2:nclus) {
    elbow[k-1,1] <- k
    elbow[k-1,2] <- kmeans(as.matrix(clustering_matrix),k,iter.max = 10, nstart = 4)$tot.withinss
  }
  
  elbow <- as.data.frame(elbow)
  colnames(elbow) <- c("k","wss")
  
  # Elbow plot
  
  ggplot(elbow, aes(x = k, y = wss)) +
    geom_point() +
    geom_line()
  
}

# Calualting the gap statistic for optimal cluster
gap_statistic <- clusGap(as.matrix(clustering_matrix), FUN = kmeans, iter.max = 10, nstart = 25, K.max = 15, B = 500)
