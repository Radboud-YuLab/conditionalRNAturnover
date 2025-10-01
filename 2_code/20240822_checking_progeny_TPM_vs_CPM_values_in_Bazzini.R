library(tidyverse)
library(tidymodels)
library(here)

library(stringr)
library(biomaRt)
library(progeny)

conflicted::conflict_prefer("select", winner = "dplyr")
conflicted::conflict_prefer("filter", winner = "dplyr")

# get latest ensembl for the hgnc symbols
ensembl <- useEnsembl(biomart = 'genes', 
                      dataset = 'hsapiens_gene_ensembl',
                      version = 112)

attributes <- listAttributes(ensembl)

bazzini <- read.csv("data_KC_A//raw_data/Bazzini/GSE126522_time_course_slam_CPM.csv")
features_all.genes <- read.csv("data_KC_A/raw_data/raw_data_from_janou/Features_genes.csv")

bazzini_biomart <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                          filters = "ensembl_gene_id",
                          values = bazzini$Gene_ID,
                          mart = ensembl)
# only take the 0 hours
bazzini <- bazzini %>% select(Gene_ID, k562_0h_A, k562_0h_B, k562_0h_C)

# renaming gene id
colnames(bazzini)[1] <- "ensembl_gene_id"
bazzini <- column_to_rownames(bazzini,var = "ensembl_gene_id")

# merge dataset
features_all.genes_transcript.length <- features_all.genes %>% select(ensembl_gene_id,transcript_length)
bazzini <- rownames_to_column(.data = bazzini,"ensembl_gene_id")
bazzini_w.transcript.length <- merge(bazzini, features_all.genes_transcript.length, by = "ensembl_gene_id")

bazzini_length <- bazzini_w.transcript.length %>% select(transcript_length)
bazzini_names <- bazzini_w.transcript.length %>% select(ensembl_gene_id)
bazzini_vals <- bazzini_w.transcript.length[, c(2,3,4)]

# create as vector (same dimensions)
bazzini_length <- as.vector(bazzini_length$transcript_length)
bazzini_vals <- as.matrix(bazzini_vals)

length_normalized_counts <- bazzini_vals / bazzini_length
tpm <- sweep(length_normalized_counts, 2, colSums(length_normalized_counts), FUN = "/") * 1e6

colSums(tpm)

tpm <- as.data.frame(tpm)
bazzini_tpm <- cbind(bazzini_names, tpm)

# checking difference in progeny scores 
bazzini_tpm_hgnc <- merge(bazzini_tpm, features_all.genes %>% select(ensembl_gene_id, external_gene_name))

# filter out replicate and empty values in "external_gene_name"

if (length(unique(bazzini_tpm_hgnc$external_gene_name)) < nrow(bazzini_tpm_hgnc)) {
  print("There are duplicate names.")
  duplicates <- bazzini_tpm_hgnc[duplicated(bazzini_tpm_hgnc$external_gene_name) | duplicated(bazzini_tpm_hgnc$external_gene_name, fromLast = TRUE), ]
} else {
  print("there are no duplicate names")
}

bazzini_tpm_hgnc<-bazzini_tpm_hgnc[!(bazzini_tpm_hgnc$external_gene_name %in% unique(duplicates$external_gene_name)),]
rownames(bazzini_tpm_hgnc) <- NULL
bazzini_tpm_hgnc <- bazzini_tpm_hgnc %>% select(-ensembl_gene_id)
bazzini_tpm_hgnc <- column_to_rownames(bazzini_tpm_hgnc, var = "external_gene_name")

bazzini_tpm_hgnc <- as.matrix(bazzini_tpm_hgnc)
bazzini_tpm_hgnc
progeny_result_bazzini_tpm_hgnc <- progeny(bazzini_tpm_hgnc)

# cpm data bazzini
bazzini_cpm <- bazzini
bazzini_cpm <- bazzini_cpm %>% rownames_to_column("ensembl_gene_id")
bazzini_cpm_hgnc <- merge(bazzini_cpm, features_all.genes %>% select(ensembl_gene_id,external_gene_name))
bazzini_cpm_hgnc <- bazzini_cpm_hgnc %>% select(-ensembl_gene_id)


if (length(unique(bazzini_cpm_hgnc$external_gene_name)) < nrow(bazzini_cpm_hgnc)) {
  print("There are duplicate names.")
  duplicates <- bazzini_cpm_hgnc[duplicated(bazzini_cpm_hgnc$external_gene_name) | duplicated(bazzini_cpm_hgnc$external_gene_name, fromLast = TRUE), ]
} else {
  print("there are no duplicate names")
}
bazzini_cpm_hgnc<-bazzini_cpm_hgnc[!(bazzini_cpm_hgnc$external_gene_name %in% unique(duplicates$external_gene_name)),]

rownames(bazzini_cpm_hgnc) <- NULL
bazzini_cpm_hgnc <- column_to_rownames(bazzini_cpm_hgnc, var = "external_gene_name")


bazzini_cpm_hgnc <- as.matrix(bazzini_cpm_hgnc)
bazzini_cpm_hgnc
progeny_result_bazzini_cpm_hgnc <- progeny(bazzini_cpm_hgnc)

saveRDS(here("/4_processed_data/wholeset/progeny/bazzini/20240822_Bazzini_progeny_cpm.RDS"))
saveRDS(here("/4_processed_data/wholeset/progeny/bazzini/20240822_Bazzini_progeny_tpm.RDS"))

# not that much of a difference in cpm vs tpm values
