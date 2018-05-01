##############     01_import_data.R     ##############

# A program that imports the Salmon quant.sf files for all the samples and converts them to a form that is accesible by DESeq2

################
#### Set up ####
################

# Libraries required
library(tidyverse)
library(stringr)
library(tximport)
library(readr)
library(DESeq2)

# Set working directory
setwd("//slcu.cam.ac.uk/data/TeamOL/personal_folders/gregt/bud_rnaseq/")

# Load information from the HPC
# All samples
samples <- read.csv("sample_information.csv", stringsAsFactors=FALSE)

# Removed low quality sample (LMC4220_lowplast_LN_rep1)
samples <- samples[!grepl("LMC4220_lowplast_LN_rep1", samples$sample_name),] 



#####################################
#### Create the tx2gene database ####
#####################################

# Set up path
dir <- "gene_quantification"
files <- file.path(dir, samples$sample_name, "quant.sf")

# Read the quatn.sf files and bind them together using the sample id
quant <- map(files, read_tsv) 
names(quant) <- samples$sample_name

quant <- bind_rows(quant, .id = "sample") %>% 
  separate(sample, c("genotype", "plasticity", "nitrate", "rep"), 
           sep = "_", remove = FALSE) %>%
  mutate(gene_id = gsub("\\..*", "", Name)) %>%
  mutate(gene_id = gsub(".*:","",gene_id)) %>%
  rename(transcript = Name)

# Create tx2gene database, one column of transcipt ID and the other of the corresponding gene ID
tx2gene <- quant %>% 
  distinct(transcript, gene_id)


####################################################
##### Create the txi file from the tx2gene file ####
####################################################

# Import using tximport
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

colnames(txi$abundance) <- str_extract(files, "/.*/") %>% str_replace_all("/", "")
colnames(txi$counts) <- str_extract(files, "/.*/") %>% str_replace_all("/", "")
colnames(txi$length) <- str_extract(files, "/.*/") %>% str_replace_all("/", "")


################################################
#### Create DESeq2 object from the txi file ####
################################################

# Create colData, a dataframe containing info for each sample
coldata <- quant %>% 
  distinct(genotype, plasticity, nitrate, rep) %>% 
  mutate(genotype = factor(genotype),
         plasticity = factor(plasticity, levels = c("lowplast", "highplast")),
         nitrate = factor(nitrate, levels = c("LN", "HN")), replicate = factor(rep, levels = c("rep1","rep2")))

rownames(coldata) <- str_extract(files, "/.*/") %>% str_replace_all("/", "")

# Add individual ID per plasticity group
coldata <- coldata %>% 
  group_by(plasticity) %>% 
  mutate(ind_n = factor(as.numeric(factor(genotype))))


# Create the DESeq object - using full txi with counts and avgTxLength
dds <- DESeqDataSetFromTximport(txi, colData = coldata, 
                                     design = ~ plasticity + nitrate + plasticity:nitrate)

# Create set of genes that have at least one read present
dds_full <- dds[rowSums(counts(dds)) > 0, ]

full_set <- as.data.frame(rownames(dds_full))
colnames(full_set) <- ("row")

# Filter out lowly expressed genes
dds <- dds_full[ rowSums((counts(dds_full) > 10)) > 6, ]

filtered_set <- as.data.frame(rownames(dds))
colnames(filtered_set) <- ("row")


##################
# Saving objects #
##################

write_rds(txi, "scripts/R/objects/txi.rds")
write_rds(quant, "scripts/R/objects/quant.rds")
write_rds(coldata, "scripts/R/objects/coldata.rds")
write_rds(full_set, "scripts/R/objects/full_set.rds")
write_rds(filtered_set, "scripts/R/objects/filtered_set.rds")
write_rds(dds_full, "scripts/R/objects/dds_full.rds")
write_rds(dds, "scripts/R/objects/dds.rds")
