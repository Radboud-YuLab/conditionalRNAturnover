library(here)
library(tidyverse)
library(tidymodels)

# ADDING TARGETSCAN DATA

ATG_features_kera <- readRDS(here("4_processed_data/keratinocytes/ATG_start_data.RDS"))

targetscan_predictions <- read.csv(here("1_raw_data/Agarwal/CWCS.txt"), sep = "\t")

# change column names by adding identifiers
targetscan_predictions_with_ident <- targetscan_predictions
colnames(targetscan_predictions_with_ident) <- gsub("\\.", "_", colnames(targetscan_predictions))
colnames(targetscan_predictions_with_ident)[-1] <- paste("genetic.miRNA",colnames(targetscan_predictions_with_ident)[-1], sep = ".") 
colnames(targetscan_predictions_with_ident)[1] <- "ensembl_gene_id"

# negating the values 
targetscan_predictions_with_ident_abs <- mutate_if(targetscan_predictions_with_ident, is.numeric, abs)

# joining the data 
genetic_miRNA <- inner_join(ATG_features_kera, targetscan_predictions_with_ident_abs, by = "ensembl_gene_id")

# selecting only identifier and faetures
genetic_miRNA_selection <- genetic_miRNA %>% select(ensembl_gene_id, contains("genetic"))

saveRDS(genetic_miRNA_selection, file = here("4_processed_data/keratinocytes/genetic_miRNA.RDS"))
