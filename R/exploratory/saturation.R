##############     saturation.R     ##############

# A program that creates saturation plots for transcripts and genes. It also has histogram plots for read count and unique transcripts and genes.

################
#### Set up ####
################

# Libraries required
library(tidyverse)
library(ggplot2)

# Set working directory
setwd("//slcu.cam.ac.uk/data/TeamOL/personal_folders/gregt/bud_rnaseq/")

# Load information from the HPC - use all samples
samples_sat <- read.csv("sample_information.csv", stringsAsFactors=FALSE)

dir <- "gene_quantification"
files <- file.path(dir, samples_sat$sample_name, "quant.sf")

# Read the quatn.sf files and bind them together using the sample id
quant_sat <- map(files, read_tsv) 
names(quant_sat) <- samples_sat$sample_name

quant_sat <- bind_rows(quant_sat, .id = "sample") %>% 
  separate(sample, c("genotype", "plasticity", "nitrate", "rep"), 
           sep = "_", remove = FALSE) %>%
  mutate(gene_id = gsub("\\..*", "", Name)) %>%
  mutate(gene_id = gsub(".*:","",gene_id)) %>%
  rename(transcript = Name)


#################################
##### Set up initial values #####
#################################

# Input values
n_freq <- 15   # Number of points to plot
start <- 0     # Start value

# Creating probabilities table
prob_transcript <- quant_sat %>% 
  group_by(sample) %>% 
  mutate(probability = NumReads/sum(NumReads)) %>% 
  ungroup()


###############################################
############# Unique transcripts ##############
###############################################

tot_transcripts <- round(nrow(prob_transcript)/47)

# Creating reads which we are mapping at
reads_transcript <- prob_transcript %>% 
  group_by(sample) %>%
  summarise(total = sum(NumReads)) %>%
  ungroup()

reads_transcript <- reads_transcript[rep(seq_len(nrow(reads_transcript)), n_freq), ] %>%
  group_by(sample) %>%
  mutate(i = 1:n_freq) %>%
  ungroup() %>%
  mutate(num_read = round(start + (total-start)*(i-1)/(n_freq-1))) %>%
  mutate(unique_transcripts = 0) %>%
  mutate(percent_transcripts = 0)

# Counting the unique transcripts

for(i in 1:nrow(reads_transcript)){
  
  i_sample <- reads_transcript[[i,"sample"]]
  i_reads <- reads_transcript[[i,"num_read"]]
  
  x <- prob_transcript %>%
    filter(sample == i_sample) %>%
    mutate(bin = i_reads*probability >= 10) %>% 
    summarise(n = sum(bin))
  
  reads_transcript[i,"unique_transcripts"] = x$n
  reads_transcript[i,"percent_transcripts"] = reads_transcript[i,"unique_transcripts"]*100/tot_transcripts
}


###############################################
############### Unique genes ##################
###############################################

prob_gene <- prob_transcript %>%
  group_by(sample,gene_id) %>%
  summarise(probability = sum(probability), NumReads = sum(NumReads))

tot_genes <- round(nrow(prob_gene)/47)

reads_gene <- prob_gene %>%
  group_by(sample) %>%
  summarise(total = sum(NumReads)) %>%
  ungroup()

reads_gene <- reads_gene[rep(seq_len(nrow(reads_gene)), n_freq), ] %>%
  group_by(sample) %>%
  mutate(i = 1:n_freq) %>%
  ungroup() %>%
  mutate(num_read = round(start + (total-start)*(i-1)/(n_freq-1))) %>%
  mutate(unique_genes = 0) %>%
  mutate(percent_genes = 0)


# Counting the unique genes

for(i in 1:nrow(reads_gene)){
  
  i_sample <- reads_gene[[i,"sample"]]
  i_reads <- reads_gene[[i,"num_read"]]
  
  x <- prob_gene %>%
    filter(sample == i_sample) %>%
    mutate(bin = i_reads*probability >= 10) %>% 
    summarise(n = sum(bin))
  
  reads_gene[i,"unique_genes"] = x$n
  reads_gene[i,"percent_genes"] = reads_gene[i,"unique_genes"]*100/tot_genes
}


#############################################
############ Plotting results ###############
#############################################

# Saturation curve of unique transcripts
ggplot(data=reads_transcript, 
       aes(num_read, percent_transcripts, group=sample, colour=sample)) +
  geom_line() +
  geom_point() +
  geom_point(data=subset(reads_transcript, i %in% c(n_freq)), aes(num_read, percent_transcripts), shape = 21, fill = 'white', size=2) +
  theme(legend.position="none") +
  labs(x = "Reads", y = "Unique transcripts") 

# Saturation curve of percent of unique transcripts mapped
ggplot(data=reads_transcript, 
       aes(num_read, unique_transcripts, group=sample, colour=sample)) +
  geom_line() +
  geom_point() +
  geom_point(data=subset(reads_transcript, i %in% c(n_freq)), aes(num_read, unique_transcripts), shape = 21, fill = 'white', size=2) +
  theme(legend.position="none") +
  labs(x = "Reads", y = "Percentage of unique transcripts") 

# Histogram of unique transcripts
ggplot(subset(reads_transcript, i %in% c(n_freq))) +
  geom_histogram(aes(unique_transcripts)) 



reads_gene <- reads_gene %>% 
  separate(sample, c("genotype", "plasticity", "nitrate", "rep"), 
           sep = "_", remove = FALSE)

# Saturation curve of unique genes
ggplot(data=reads_gene, 
       aes(num_read, unique_genes, group=sample, fill = interaction(plasticity,nitrate))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_point(data=subset(reads_gene, i %in% c(n_freq)), aes(num_read, unique_genes), shape = 21, size=4) +
  theme(legend.position="none") +
  labs(x = "Reads", y = "Unique genes") +
  scale_colour_manual(values=c("black"))+
  scale_fill_manual(values=c("red", "deepskyblue4", "magenta", "deepskyblue"), 
                    name  ="Group", 
                    labels=c("High Nitrate, Plastic","High Nitrate, Non-Plastic",  
                             "Low Nitrate, Plastic", "Low Nitrate, Non-Plastic")) +
  theme_bw() +
  theme(panel.background = element_blank(), legend.position = c(0.9, 0.8)) 



# Saturation curve of percent of unique genes mapped
ggplot(data=reads_gene, 
       aes(num_read, percent_genes, group=sample, colour=sample)) +
  geom_line() +
  geom_point() +
  geom_point(data=subset(reads_gene, i %in% c(n_freq)), aes(num_read, percent_genes), shape = 21, fill = 'white', size=2) +
  theme(legend.position="none") +
  labs(x = "Reads", y = "Percentage of genes")

# Histogram of unique genes
ggplot(subset(reads_gene, i %in% c(n_freq))) +
  geom_histogram(aes(unique_genes)) 



# Histogram of number of reads
ggplot(subset(reads_transcript, i %in% c(n_freq))) +
  geom_histogram(aes(num_read)) +
  scale_x_log10()



##################
# Saving objects #
##################

write_rds(quant_sat, "scripts/R/objects/quant_sat.rds")
