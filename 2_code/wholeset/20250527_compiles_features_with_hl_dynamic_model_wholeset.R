# this script compiles, the samples of the dynamic model with all the features (basic RNA sequence features and codonfreq, progeny)
# WHOLESET version

library(here)
library(tidyverse)


# dynamic half life value with identifiers dataset
hl <- readRDS(file = here("4_processed_data/wholeset/20250527_dynamic_hl_4sU_cutoff6.RDS"))
samples_hl <- unique(hl$ensembl_gene_id)

# append the static features
static.features <- readRDS("4_processed_data/static.features/20250530_basic.genetic.RDS") 

hl.static <- merge(x = hl, y = static.features, all = FALSE) #some genes are lost as they did not start with ATG
colnames(hl.static)
length(unique(hl.static$ensembl_gene_id)) 


# append progeny scores
progeny_scores <- readRDS("4_processed_data/wholeset/progeny/unscaled/20240910_progeny_results_all_mean_per_condition_unscaled.RDS")
progeny_scores_df <- as.data.frame(progeny_scores)

# add prefix 
colnames(progeny_scores_df) <- paste0("dynamic.progeny.", colnames(progeny_scores_df))

progeny_scores_df <- rownames_to_column(progeny_scores_df, var = "group")
progeny_scores_df$group <- stringr::str_replace(string = progeny_scores_df$group, pattern = "GM12878", replacement = "LCL")

hl.static.progeny <- merge(hl.static, progeny_scores_df, all = TRUE, by = "group")
colnames(hl.static.progeny)

saveRDS(hl.static.progeny, file = here("4_processed_data/ML_input/20250527_dynamicHL_basic.genetic.RDS"))



