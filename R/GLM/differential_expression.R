# Finding differentially expressed genes using DESeq2

library(DESeq2)
library(tidyverse)
library(readr)

# Set up
setwd("~/Part_III_Project/")

txi <- readRDS("scripts/R/objects/txi.rds")
coldata <- readRDS("scripts/R/objects/coldata.rds")

Alpha <- 0.05


################################################
#### Change in nitrate for given plasticity ####
################################################

# Create the DESeq object
dds_ind <- DESeqDataSetFromTximport(txi, colData = coldata, 
                                    design = ~ plasticity + plasticity:ind_n + 
                                      plasticity:nitrate)
# Filter data
dds_ind <- dds_ind[ rowSums((counts(dds_ind) > 10)) > 6, ]

# Apply DESeq
dds_ind_wald <- DESeq(dds_ind)


###########################################
### Change in nitate for low plasticity ###
###########################################

res_N_LP <- results(dds_ind_wald, name = "plasticitylowplast.nitrateHN", tidy = TRUE, alpha = Alpha)

############################################
### Change in nitate for high plasticity ###
############################################

res_N_HP <- results(dds_ind_wald, name = "plasticityhighplast.nitrateHN", tidy = TRUE, alpha = Alpha)

#############################################
### Change in nitate for interaction term ###
#############################################

res_int <- results(dds_ind_wald, contrast = c(list("plasticityhighplast.nitrateHN","plasticitylowplast.nitrateHN")), tidy = TRUE, alpha = Alpha)



######################################################
#### Difference in plasticity for a given nitrate ####
######################################################

# Create the DESeq object
dds_alt <- DESeqDataSetFromTximport(txi, colData = coldata, 
                                    design = ~ 0 + plasticity:nitrate)

# Filter data
dds_alt <- dds_alt[ rowSums((counts(dds_alt) > 10)) > 6, ]

# Apply DEseq
dds_alt_wald <- DESeq(dds_alt)


################################################
### Difference in plasticity for low nitrate ###
################################################

res_P_LN <- results(dds_alt_wald, contrast = list("plasticityhighplast.nitrateLN", "plasticitylowplast.nitrateLN"), tidy = TRUE ,alpha = Alpha)

#################################################
### Difference in plasticity for high nitrate ###
#################################################

res_P_HN <- results(dds_alt_wald, contrast = list("plasticityhighplast.nitrateHN", "plasticitylowplast.nitrateHN") ,alpha = Alpha, tidy = TRUE)


##################
# Saving objects #
##################

write_rds(dds_ind, "scripts/R/objects/dds_ind.rds")
write_rds(dds_ind_wald, "scripts/R/objects/dds_ind_wald.rds")
write_rds(dds_alt, "scripts/R/objects/dds_alt.rds")
write_rds(dds_alt_wald, "scripts/R/objects/dds_alt_wald.rds")
write_rds(res_N_LP, "scripts/R/objects/res_N_LP.rds")
write_rds(res_N_HP, "scripts/R/objects/res_N_HP.rds")
write_rds(res_P_LN, "scripts/R/objects/res_P_LN.rds")
write_rds(res_P_HN, "scripts/R/objects/res_P_HN.rdS")
write_rds(res_int, "scripts/R/objects/res_int.rdS")
