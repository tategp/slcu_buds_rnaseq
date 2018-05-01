# Plotting clusters

# Set up

library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cluster)

setwd("~/Part_III_Project/")

# Parameters
gene_set = "no_threshold_set"
num_clus = 16

# Load objects
gene_set <- readRDS(paste0("scripts/R/objects/",gene_set,".rds"))
colnames(gene_set) <- c("row")
x <- readRDS("scripts/R/objects/vst_zscore.rds")
x <- x[rownames(x) %in% gene_set$row,]

# Calculate correlation
x_cor <- cor(t(x))
# Express as (1 - correlation) such that
## r = 1 becomes 0 (very similar)
## r = 0 becomes 1 (quite dissimilar)
## r = -1 becomes 2 (very dissimilar)
x_dis <- 1 - x_cor
# Return as a distance object
clustering_matrix <- as.dist(x_dis)


# Get the clustering
cluster <- hclust(clustering_matrix)

# Cut the tree and lable genes
gene_clus <- as.data.frame(cutree(cluster, num_clus))
gene_clus <- as.data.frame(cbind(row.names(gene_clus),gene_clus[,1]))
colnames(gene_clus) <- c("gene","clus")
gene_clus$clus <- as.factor(as.numeric(gene_clus$clus))

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
  facet_wrap( ~ as.factor(as.numeric(clus)), ncol = 4) +
  #Aesthetics...
  labs(x = "Nitrate level", y = "Z-score") +
  scale_colour_manual(values=c("tomato", "royalblue1"), 
                      name  ="Plasticity", 
                      labels=c("Non-Plastic", "Plastic"))+
  theme_bw() +
  theme(panel.background = element_blank(),legend.position="bottom")  




#####################################

# K-means clustering
nclus <- 18

clusters <- vector("list", nclus-1)
centers <- vector("list", nclus-1)
totss <- vector("list", nclus-1)
withinss <- vector("list", nclus-1)
tot.withinss <- vector("list", nclus-1)
betweenss <- vector("list", nclus-1)
size <- vector("list", nclus-1)
iter <- vector("list", nclus-1)
ifault <- vector("list", nclus-1)

for (i in 2:nclus) {
  
  a <- kmeans(as.matrix(clustering_matrix),i,iter.max = 10, nstart = 5)
  clusters[[i-1]] <- a$cluster
  centers[[i-1]] <- a$centers
  totss[[i-1]] <- a$totss
  withinss[[i-1]] <- a$withinss
  tot.withinss[[i-1]] <- a$tot.withinss
  betweenss[[i-1]] <- a$betweenss
  size[[i-1]] <- a$size
  iter[[i-1]] <- a$iter
  ifault[[i-1]] <- a$ifault
  print(i)
}

setNames(clusters, paste0("cluster", 2:nclus))
setNames(centers, paste0("cluster", 2:nclus))
setNames(totss, paste0("cluster", 2:nclus))
setNames(withinss, paste0("cluster", 2:nclus))
setNames(tot.withinss, paste0("cluster", 2:nclus))
setNames(betweenss, paste0("cluster", 2:nclus))
setNames(size, paste0("cluster", 2:nclus))
setNames(iter, paste0("cluster", 2:nclus))
setNames(ifault, paste0("cluster", 2:nclus))




elbow <- matrix(nrow = (nclus - 1), ncol = 2)

for (k in 2:nclus) {
  elbow[k-1,1] <- k
  elbow[k-1,2] <- tot.withinss[[k]]
}


elbow <- as.data.frame(elbow)
colnames(elbow) <- c("k","wss")

# Elbow plot
ggplot(elbow, aes(x = k, y = wss)) +
  geom_point() +
  geom_line()


####################
### Save objects ###
####################

write_rds(gene_clus,"scripts/R/objects/gene_clus.rds")
