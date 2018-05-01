# Comparing data to external set

no_threshold_set <- readRDS("scripts/R/objects/no_threshold_set.rds")
full_set <- readRDS("scripts/R/objects/full_set.rds")
gene_clus <- readRDS("scripts/R/objects/gene_clus.rds")

# Load gene set
root_set <- read.csv("scripts/R/external_data/root_set.csv", header = TRUE, sep = ",")

# Remove last row
root_set <- root_set[1:(nrow(root_set)-1),]

# Turn to characters
root_set$id <- as.character(root_set$id)
# Extract gene id
id <- sapply(strsplit(root_set$id, split = '|', fixed = TRUE),'[[',3)

root_set <- root_set %>% 
  dplyr::select(Description,Nitrate.induced,Nitrate.repressed)

root_set <- cbind(row = id,gene = id,root_set)





# Overlap between root set and full set of expressed genes
expressed_root_set <- inner_join(root_set,full_set)

nrow(expressed_root_set)/nrow(root_set)

# Intersect with the differentially expressed genes
diff_root_set <- inner_join(expressed_root_set,no_threshold_set)

nrow(diff_root_set)/nrow(expressed_root_set)

# See what clusters they're in
clus_diff_root <- left_join(diff_root_set,gene_clus)

table(clus_diff_root$clus) %>% as.data.frame() %>% mutate(Freq = Freq*6183/nrow(diff_root_set)) %>% View()



#Induced genes
induced_root_set <- diff_root_set %>% 
  filter(Nitrate.induced == 1)

clus_induced_root <- left_join(induced_root_set,gene_clus)

table(clus_induced_root$clus) %>% as.data.frame() %>% 
  #mutate(Freq = Freq*6183/nrow(induced_root_set)) %>% 
  View()

#Reduced genes
reduced_root_set <- diff_root_set %>% 
  filter(Nitrate.repressed == 1)

clus_reduced_root <- left_join(reduced_root_set,gene_clus)

table(clus_reduced_root$clus) %>% as.data.frame() %>% 
  #mutate(Freq = Freq*6183/nrow(reduced_root_set)) %>% 
  View()
