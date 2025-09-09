# clean version of basic_rna_sequence_features
# will only extract BASIC features

library(here)
library(tidymodels)
library(tidyverse)
library(BiocManager)
# to download Biostrings: 
# BiocManager::install("Biostrings")
library(Biostrings)
library(biomaRt)
library(biomartr)
library(conflicted)

conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)

identifier <- "20250619"

# loading dependencies
load(here("1_raw_data/Mulder/Normalised_GSdata.RData"))
features_keratinocyte <- read.csv(here("1_raw_data/Mulder/Features_keratinocyte_HLdata.csv"))
source(here("2_code/scripts/20240417_make_hist_from_df_v2.R"))

#----------------------------#
# FILTERING DATASET

ATG_start_features_kera <- features_keratinocyte %>% 
  filter(substr(coding, 1, 3) == "ATG")


#----------------------------#
# ADDING EXON JUNCTION DENSITY
ensembl106 <- useEnsembl(biomart = 'genes',
                         dataset = 'hsapiens_gene_ensembl',
                         version = 106)

exon_sequences <- getBM(attributes = c("ensembl_transcript_id","gene_exon"),
                        filters = "ensembl_transcript_id",
                        values = ATG_start_features_kera$ensembl_transcript_id,
                        mart = ensembl106)

length(unique(exon_sequences$ensembl_transcript_id))


grouped_exon_sequences <- exon_sequences %>%
  group_by(ensembl_transcript_id) %>%   # only the number of exons are needed
  summarize(gene_exon = list(gene_exon)) %>%  # combines the exon sequences into a list per ensembl_transcript_id 
  mutate(
    exon_count = lengths(gene_exon), #adds the number of exons
    exon_junction_count = case_when( # adds the number of exon junctions
      exon_count == 1 ~ exon_count, # ignore exon counts of 1 to exclude 0 values
      TRUE ~ exon_count - 1))

# get a df with only the exon junctions
ejc_exon_counts <- grouped_exon_sequences %>% 
  select(-gene_exon)

# merge in order to get df with codon frequences and the exon junction density
features_keratinocyte_ejc_exon_count <- inner_join(features_keratinocyte, ejc_exon_counts, by = "ensembl_transcript_id")

features_keratinocyte_exon_ejc_count_density <- features_keratinocyte_ejc_exon_count %>% 
  mutate(
    ejc_density = features_keratinocyte_ejc_exon_count$exon_junction_count/(features_keratinocyte_ejc_exon_count$cds_length/1000)
  )

basic <- features_keratinocyte_exon_ejc_count_density

#----------------------------#
# ADD IDENTIFIER TO BASIC FEATURES
basic_features <- c("percentage_gene_gc_content", "transcript_length", "cds_length"
                    , "UTR5_length", "UTR3_length","cds_GC", "UTR5_GC", "UTR3_GC", "ejc_density")
added_identifiers <- paste("basic", basic_features, sep = ".")
colnames(basic)[colnames(basic) %in% basic_features] <- added_identifiers

colnames(basic)

saveRDS(basic, file = here("4_processed_data/keratinocytes/basic_rna_features.RDS"))

#----------------------------# 
# log transformations can be done when making a recipe

