# genetic_rna_sequence_features
# extracts all GENETIC RNA sequence features

library(here)
library(tidymodels)
library(tidyverse)
library(conflicted)

conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)

#----------------------------#
# ADDING BASIC FEATURES
basic <- readRDS(file = here("4_processed_data/keratinocytes/basic_rna_features.RDS"))

# ADDING CODON FREQ
genetic.codonfreq <- readRDS(file = here("4_processed_data/keratinocytes/genetic_codonfreq.RDS"))

# ADDING KMERS
genetic.kmers <- readRDS(file = here("4_processed_data/keratinocytes/genetic_kmers_3UTR.RDS"))

# ADDING TARGETSCAN
genetic.miRNA <- readRDS(file = here("4_processed_data/keratinocytes/genetic_miRNA.RDS"))

# ADDING SEQWEAVER 
genetic.seqw <- readRDS(file = "4_processed_data/keratinocytes/genetic_seqweaver.RDS")

#----------------------------#

basic.genetic <- inner_join(basic, genetic.codonfreq, by = "ensembl_gene_id") %>% 
  inner_join(genetic.kmers, by = "ensembl_gene_id") %>% 
  inner_join(genetic.miRNA, by = "ensembl_gene_id") %>% 
  inner_join(genetic.seqw, by = "ensembl_gene_id")

saveRDS(basic.genetic, file = here("4_processed_data/keratinocytes/basic.genetic.RDS"))
