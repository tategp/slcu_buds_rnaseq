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

dir <- "gene_quantification"
files <- file.path(dir, samples$sample_name, "quant.sf")


#####################################
#### Create the tx2gene database ####
#####################################

# Read the quatn.sf files and bind them together using the sample id
quant <- map(files, read_tsv) 
names(quant) <- samples$sample_name

quant <- bind_rows(quant, .id = "sample") %>% 
  separate(sample, c("genotype", "plasticity", "nitrate", "rep"), 
           sep = "_", remove = FALSE) %>%
  mutate(gene_id = gsub("\\..*", "", Name)) %>%
  mutate(gene_id = gsub(".*:","",gene_id)) %>%
  rename(Name = "transcript")

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

# Saving txi for later use in DESeq2
write_rds(txi, "scripts/R/objects")
write_rds(quant, "scripts/R/objects")


################################################
#### Create DESeq2 object from the txi file ####
################################################

# Create colData, a dataframe containing info for each sample
colData <- quant %>% 
  distinct(genotype, plasticity, nitrate, rep) 

rownames(colData) <- str_extract(files, "/.*/") %>% str_replace_all("/", "")

# Create the DESeq object - using full txi with counts and avgTxLength
dds <- DESeqDataSetFromTximport(txi, colData = colData, design = ~ 1)

# Counts only - note assays(dds)$counts and counts(dds) are identical, use when program takes the wrong assay (ie when applying vsd())
dds_counts <- DESeqDataSetFromMatrix(round(counts(dds)), colData = colData, design = ~ 1)

# Filter out lowly expressed genes
dds_filtered <- dds[ rowSums((counts(dds) > 10)) > 5, ]
dds_counts_filtered <- dds_counts[ rowSums((counts(dds) > 10)) > 5, ]

# Saving objects
write_rds(dds_filtered, "scripts/R/objects")
write_rds(dds_counts_filtered, "scripts/R/objects")




##########################################################################
#### Comparing counts from tximport to summing transcripts from quant ####
##########################################################################

salmon_counts <- txi$counts %>% 
  as.data.frame() %>% 
  dplyr::mutate(gene_id = rownames(.)) %>% 
  gather(sample, counts_txi, LMA3945_lowplast_HN_rep1:acc86_lowplast_LN_rep2)

raw_counts <- quant %>%
  group_by(sample,gene_id) %>%
  summarise(couts_raw = sum(NumReads))

compare_counts <- full_join(salmon_counts, raw_counts)
compare_counts$difference <- compare_counts$counts_txi - compare_counts$couts_raw

compare_counts <- compare_counts %>%
  distinct(difference)

max_dif <- max(compare_counts$difference)







# To go direct from the gff file

# library(GenomicFeatures)
# # Extracting the information from the gff3 file used for Salmon annotation
# gff <- "annotation/Arabidopsis_thaliana.TAIR10.37.gff3"
# txdb <- makeTxDbFromGFF(gff, format = "gff3")
# k <- keys(txdb, keytype = "TXNAME")
# df <- select(txdb, keys = k, keytype = "TXNAME", columns = "GENEID")
# df <- df %>% 
#   mutate(TXNAME = ifelse(str_detect(TXNAME, "transcript:"), 
#                          TXNAME, paste0("transcrit:", TXNAME)),
#          GENEID = str_replace()) ### Need to update here

# tx2gene <- df[, 2:1]

