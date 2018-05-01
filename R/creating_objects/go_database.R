# Creating the GO object for use in other tests

library(tidyverse)
library(biomaRt)

setwd("~/Part_III_Project/")

# Specify the dataset we want to use - Arabidopsis Thaliana gene set 
m <- useMart("plants_mart", host="plants.ensembl.org", dataset="athaliana_eg_gene")

# Get attributes from bioMart
go <- getBM(attributes= c("tair_locus", "tair_symbol", "external_gene_name", "gene_biotype",
                          "chromosome_name", "start_position", "end_position",
                          "go_id", "name_1006", "definition_1006",  "description"), mart=m)

# Filter data for protein coding genes and ones with loci
go <- go %>% 
  filter(tair_locus != "" & gene_biotype == "protein_coding")

# Save object
write_rds(go,"scripts/R/objects/go.rds")
