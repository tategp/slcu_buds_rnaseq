# GO enrichment

# Set up

library(tidyverse)
library(goseq)

setwd("~/Part_III_Project/")


# Call R objects
go <- readRDS("scripts/R/objects/go.rds")
full_set <- readRDS("scripts/R/objects/full_set.rds")

#Get length vector
go_length <- go %>% 
  filter(go_id != "" & tair_locus != "") %>% 
  group_by(tair_locus) %>% 
  summarise(length = mean(end_position) - mean(start_position))

#Only consider set of genes we have reads for
go_length <- go_length[which(go_length$tair_locus %in% full_set$row),]

#Set of genes which protein coding and we have reads for
universe_set <- go_length$tair_locus

#Length vector
go_length <- as.integer(go_length$length)


###############################
### Define gene set to test ###
###############################

#gene_set <- sig_sd_change %>% 
#  filter(N_LP == 0) %>% 
#  select(row = gene)

gene_set <- as.data.frame(clus_genes[[2]]) %>% 
  dplyr::select(row = gene)

###############################

#Differentially expressed genes
gene <- ifelse(universe_set %in% gene_set$row,1,0)
names(gene) <- universe_set


#GO identiviers
go_identifier <- go %>% 
  dplyr::select(tair_locus,go_id) %>% 
  filter(tair_locus %in% universe_set)



#GO enrichment
pwf <- nullp(gene,bias.data = go_length)

go_wall <- goseq(pwf, gene2cat = go_identifier, method = "Wallenius") %>% 
  filter(ontology == "BP")

#Adjusted pvalues
go_enriched <- go_wall %>% 
  mutate(over_padj = p.adjust(go_wall$over_represented_pvalue, method="BH")) %>% 
  mutate(under_padj = p.adjust(go_wall$under_represented_pvalue, method="BH"))

