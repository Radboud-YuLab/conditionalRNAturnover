# goal of this script is to compile the hl datasets with their static features.
# first determine whether threshold cutoff reduces the dataset too much

library(here)
library(tidyverse)

## 1. Determine reduction in number of genes when using a threshshold cutoff
# half life value with identifiers dataset
hl_no.cutoff <- readRDS(file = here("4_processed_data/wholeset/no_cutoff/20250530_hl_4sU_mean.per.condition.RDS"))
hl <- hl_no.cutoff
# extract vals mean per condition 
hl_vals <- hl[,-1]

# calculate the mean NOTE: there is no cutoff here
hl_mean <- rowMeans(hl_vals, na.rm = T)
hl_idents <- hl[,1]

hl_mean_idents <- data.frame(
  identifier = hl_idents,
  hl_mean = hl_mean
)

b.g.t.s <- readRDS("4_processed_data/static.features/20250530_basic.genetic.ts.seqw.RDS")
no_cutoff  <-  merge(x = hl_mean_idents, y = b.g.t.s, all = FALSE)
length(unique(no_cutoff$ensembl_gene_id)) # 10461

hl_w.cutoff <- readRDS(file = here("4_processed_data/wholeset/w_cutoff/20250527_static_hl_4sU_cutoff6.RDS"))

w_cutoff  <-  merge(x = hl_w.cutoff, y = b.g.t.s, all = FALSE)
length(unique(w_cutoff$ensembl_gene_id)) # 6936

# too high of a loss in genes, continuing with no threshold

## 2. create basic.codonfreq hl dataframe, static.all already created as "no_cutoff"
saveRDS(no_cutoff, file = here("4_processed_data/wholeset/no_cutoff/20250618_b.g.t.s.HL_ML.input.RDS"))

# append the static features
hl <- readRDS(here("4_processed_data/wholeset/no_cutoff/20250530_hl_4sU_mean.per.condition.RDS"))
hl_vals <- hl[,-1]

# calculate the mean NOTE: there is no cutoff here
hl_mean <- rowMeans(hl_vals, na.rm = T)
hl_idents <- hl[,1]

hl_mean_idents <- data.frame(
  identifier = hl_idents,
  hl_mean = hl_mean
)

b.g <- readRDS(here("4_processed_data/wholeset/static.features/20250530_basic.genetic.RDS"))

# create basic.codon df 
basic.codon <- merge(x = hl_mean_idents, y = b.g, all = FALSE) #some genes are lost as they did not start with ATG
colnames(basic.codon)
length(unique(basic.codon$ensembl_gene_id))

# needs to be the same number of genes as no_cutoff in order to make a fair comparison
basic.codon_final <- basic.codon[basic.codon$ensembl_gene_id %in% no_cutoff$ensembl_gene_id,] 
saveRDS(basic.codon_final, here("4_processed_data/wholeset/no_cutoff/20250618_b.g.HL_ML.input.RDS"))

setdiff(basic.codon_final$ensembl_gene_id, b.g.t.s$ensembl_gene_id)
colnames(basic.codon_final)
colnames(no_cutoff)

length(colnames(basic.codon_final))
length(colnames(no_cutoff))
## 3. After evaluation (optimization step 2), b.g. is better compared to static all. Make new df with a cutoff 

hl <- readRDS(here("4_processed_data/wholeset/w_cutoff/20250527_static_hl_4sU_cutoff6.RDS"))
hl
b.g <- readRDS(here("4_processed_data/wholeset/static.features/20250530_basic.genetic.RDS"))
b.g.w_cutoff  <-  merge(x = hl, y = b.g, all = FALSE) # some genes lost as coding sequence  did not start with ATG
saveRDS(b.g.w_cutoff, file = here("4_processed_data/ML_input/20250703_b.g.w.cutoff_static.4sU_ML.input.RDS"))



