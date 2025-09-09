library(tidyverse)
library(tidymodels)
library(here)

ATG_features_kera <- readRDS(here("4_processed_data/keratinocytes/ATG_start_data.RDS"))


seqweaver_3UTR <- read.csv(here("1_raw_data/Agarwal/3pUTR_avg.txt"), sep = "\t")
seqweaver_5UTR <- read.csv(here("1_raw_data/Agarwal/5pUTR_avg.txt"), sep = "\t")
seqweaver_orf <- read.csv(here("1_raw_data/Agarwal/ORF_avg.txt"), sep = "\t")

# adding identifiers to seqweaver and renaming identifier
seqweaver_3UTR_renamed <- seqweaver_3UTR
seqweaver_5UTR_renamed <- seqweaver_5UTR
seqweaver_orf_renamed <- seqweaver_orf

colnames(seqweaver_3UTR_renamed)[1] <- "ensembl_gene_id"
colnames(seqweaver_5UTR_renamed)[1] <- "ensembl_gene_id"
colnames(seqweaver_orf_renamed)[1] <- "ensembl_gene_id"

colnames(seqweaver_3UTR_renamed)[-1] <- paste("genetic.seqweaver_3UTR", colnames(seqweaver_3UTR_renamed)[-1], sep = ".")
colnames(seqweaver_5UTR_renamed)[-1] <- paste("genetic.seqweaver_5UTR", colnames(seqweaver_5UTR_renamed)[-1], sep = ".")
colnames(seqweaver_orf_renamed)[-1] <- paste("genetic.seqweaver_orf", colnames(seqweaver_orf_renamed)[-1], sep = ".")

# joining with keratinocyte data
genetic_seqw <- inner_join(ATG_features_kera, seqweaver_3UTR_renamed, by = "ensembl_gene_id") %>% 
  inner_join(seqweaver_5UTR_renamed, by = "ensembl_gene_id") %>% 
  inner_join(seqweaver_orf_renamed, by = "ensembl_gene_id")

# only select the identifier and features
genetic_seqw.selection <- genetic_seqw %>% select(ensembl_gene_id, contains("genetic"))
  
# save to file 
saveRDS(genetic_seqw.selection, file = here("4_processed_data/keratinocytes/genetic_seqweaver.RDS"))

