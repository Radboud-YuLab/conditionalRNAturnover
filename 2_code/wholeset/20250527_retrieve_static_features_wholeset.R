# script for generating all features
# NOTE: this was done for the wholeset 

library(tidyverse)
library(here)
library(biomaRt)
library(here)
library(Biostrings)
library(conflicted)

conflict_prefer(name = "select", winner = "dplyr")
conflict_prefer(name = "filter", winner = "dplyr")
# extract samplenames
hl_df <- readRDS(here("4_processed_data/20250516_hl_4sU.RDS"))
samples <- row.names(hl_df)

# 1 adding Biomart features to the dataset
# biomart feature list
# df generated in 2_code/Seq_features.R
# includes sequences, GC content, lengths 
features_all.genes <- read.csv("5_supplementary_data/Features_genes.csv")
df_biomart <- features_all.genes[features_all.genes$ensembl_gene_id %in% samples,] # three genes lost

colnames(df_biomart)

# remove genes that do not start with ATG
df_biomart <- df_biomart %>% 
  filter(substr(coding, 1, 3) == "ATG")

# 2. calculate and add exon junction density to the df
# retrieve exon sequences
ensembl106 <- useEnsembl(biomart = 'genes',
                         dataset = 'hsapiens_gene_ensembl',
                         version = 106)

exon_sequences <- getBM(attributes = c("ensembl_transcript_id","gene_exon"),
                        filters = "ensembl_transcript_id",
                        values = df_biomart$ensembl_transcript_id,
                        mart = ensembl106)
saveRDS(exon_sequences, here("5_supplementary_data/wholeset_exon_sequences_biomart106.RDS")
exon_sequences <- readRDS("5_supplementary_data/wholeset_exon_sequences_biomart106.RDS") # use above script to retrieve, but takes some time
length(unique(exon_sequences$ensembl_transcript_id))

# retrieve number of exons per transcript
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
# no genes lost
df_biomart_ejc <- inner_join(df_biomart, ejc_exon_counts, by = "ensembl_transcript_id")

# calculate density
df_biomart_ejc.density <- df_biomart_ejc %>% 
  mutate(
    ejc_density = df_biomart_ejc$exon_junction_count/(df_biomart_ejc$cds_length/1000)
  )

# add basic identifier
basic_features <- c("percentage_gene_gc_content", "transcript_length", "cds_length"
                    , "UTR5_length", "UTR3_length","cds_GC", "UTR5_GC", "UTR3_GC", "ejc_density")
added_identifiers <- paste("basic", basic_features, sep = ".")
colnames(df_biomart_ejc.density)[colnames(df_biomart_ejc.density) %in% basic_features] <- added_identifiers

# 3. add codonfrequency

df_biomart_ejc.density_codonfreq <- cbind(df_biomart_ejc.density, do.call(rbind, lapply(df_biomart_ejc.density$coding, function(x){
  if (length(x) > 0) {
    y = oligonucleotideFrequency(DNAStringSet(x), 3, step=3)
    colnames(y)=paste0("genetic.codon.",colnames(y)) 
    y/sum(y)
  }
})))

stop_codons <- c("TAA","TAG","TGA")
stop_codons <- paste("genetic.codon", stop_codons, sep = ".")
removed_stop_codons <-  df_biomart_ejc.density_codonfreq %>% # make dataset that don't take stop codon frequencies into account
  select(-all_of(stop_codons))

first_col = which(colnames(removed_stop_codons) == "genetic.codon.AAA")
last_col = which(colnames(removed_stop_codons) == "genetic.codon.TTT")
removed_stop_codons$sums <- rowSums(removed_stop_codons[, first_col:last_col])


# normalization step
removed_stop_codons <- removed_stop_codons %>% 
  mutate(across("genetic.codon.AAA":"genetic.codon.TTT", ~ .x/sums))

basic.genetic <- removed_stop_codons %>% select(-sums)

# check for NA:
colnames(basic.genetic)[colSums(is.na(basic.genetic)) > 0] # some genes do not have a UTR3 and UTR5
NA_containing_data <- basic.genetic[rowSums(is.na(basic.genetic))>0,]
# turn those NA to 0
basic.genetic[is.na(basic.genetic)] <- 0
colSums(is.na(basic.genetic))

saveRDS(basic.genetic, file = here("4_processed_data/20250527_wholeset_basic.genetic.RDS"))
