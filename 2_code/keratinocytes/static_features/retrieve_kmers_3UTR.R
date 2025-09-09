library(here)
library(tidymodels)
library(tidyverse)
library(BiocManager)
library(conflicted)
library(Biostrings)

conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)

identifier <- "20250624"

# load filtered features
ATG_start_features_kera <- readRDS("4_processed_data/keratinocytes/ATG_start_data.RDS")

#------------------------#
# ADDING 1TO7 KMERS OF 3UTR

kmer_len <- 7
genetic_kmers <- ATG_start_features_kera  # Initialize with original dataframe

for (k in 1:kmer_len) {
  # Calculate k-mer frequencies for each DNA sequence in the 'coding' column
  kmer_frequencies <- lapply(ATG_start_features_kera$coding, function(x) {
    if (length(x) > 0) {
      y <- oligonucleotideFrequency(DNAStringSet(x), k, step = 1)
      colnames(y) <- paste0("genetic.kmer_3UTR.", colnames(y))
      y / sum(y)
    }
  })
  
  # Combine k-mer frequencies into a single dataframe
  combined_kmers <- do.call(rbind, kmer_frequencies)
  
  # Merge k-mer features with existing dataframe
  genetic_kmers <- cbind(genetic_kmers, combined_kmers)
}


# only select identifier and features
genetic.kmers_3UTR.selection <- genetic_kmers %>% select(ensembl_gene_id, contains("genetic"))

# save to file 
saveRDS(genetic.kmers_3UTR.selection, here("4_processed_data/keratinocytes/genetic_kmers_3UTR.RDS"))
