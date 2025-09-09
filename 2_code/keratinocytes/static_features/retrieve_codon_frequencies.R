# Rscripts that retrieves the codon frequencies

library(here)
library(tidymodels)
library(tidyverse)
library(BiocManager)
library(conflicted)
library(Biostrings)


conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)

identifier <- "20250619"

# load filtered features
ATG_start_features_kera <- readRDS(here("4_processed_data/keratinocytes/ATG_start_data.RDS"))

#----------------------------#
# ADDING CODON FREQ
features_kera_codon_freq <- cbind(ATG_start_features_kera, do.call(rbind, lapply(ATG_start_features_kera$coding, function(x){
  if (length(x) > 0) {
    y = oligonucleotideFrequency(DNAStringSet(x), 3, step=3)
    colnames(y)=paste0("genetic.codon.",colnames(y)) 
    y/sum(y)
  }
})))

stop_codons <- c("TAA","TAG","TGA")
stop_codons <- paste("genetic.codon", stop_codons, sep = ".")
removed_stop_codons <-  features_kera_codon_freq %>% # make dataset that don't take stop codon frequencies into account
  select(-all_of(stop_codons))

first_col = which(colnames(removed_stop_codons) == "genetic.codon.AAA")
last_col = which(colnames(removed_stop_codons) == "genetic.codon.TTT")
removed_stop_codons$sums <- rowSums(removed_stop_codons[, first_col:last_col]) # calculate the sum of all codon frequencies excluding stop codons

# normalization step
removed_stop_codons <- removed_stop_codons %>% 
  mutate(across("genetic.codon.AAA":"genetic.codon.TTT", ~ .x/sums)) 

genetic.codonfreq <- removed_stop_codons %>% select(-sums)


# only selecting geneids and codonfreq in order to inner join when compiling the data
genetic.codonfreq.selection <- genetic.codonfreq %>% select(ensembl_gene_id, contains("genetic"))

# save to file 
saveRDS(genetic.codonfreq.selection, file = here("4_processed_data/keratinocytes/genetic_codonfreq.RDS"))
