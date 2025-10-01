library(tidyverse)
library(tidymodels)
library(here)
library(stringr)
library(biomaRt)
library(progeny)

conflicted::conflict_prefer("select", winner = "dplyr")
conflicted::conflict_prefer("filter", winner = "dplyr")



remove_duplicates_df <- function(df, name_column) {
  # Check if there are duplicate names
  if(length(unique(df[[name_column]])) < nrow(df)) {
    print("There are duplicate names.")
    
    # Find duplicates
    duplicates <- df[duplicated(df[[name_column]]) | duplicated(df[[name_column]], fromLast = TRUE), ]
    
    # Remove rows with duplicate names
    df <- df[!(df[[name_column]] %in% unique(duplicates[[name_column]])), ]
  } else {
    print("There are no duplicate names.")
  }
  
  return(df)
}

# biomart database in order to get the hgnc symbol

ensembl <- useEnsembl(biomart = 'genes', 
                      dataset = 'hsapiens_gene_ensembl',
                      version = 106)

attributes <- listAttributes(ensembl)


# retrieved from "/2_code/20240822_checking_progeny_TPM_vs_CPM_values_in_Bazzini.R"
bazzini <- readRDS("4_processed_data/wholeset/progeny/bazzini/20240822_Bazzini_cpm_data.RDS")
# needs to be log transformed 
bazzini_log10 <- log10(bazzini+0.1)
bazzini_log10_progeny <- progeny(bazzini_log10, scale = F, verbose =T, organism = "Human", top = 100)
# saveRDS(bazzini_log10_progeny, file = "4_processed_data/wholeset/progeny/bazzini/20240822_Bazzini_progeny_log10_cpm.RDS")

## cramer1 
cramer1 <- read_tsv("1_raw_data/count_data/Cramer1/GSE75792_norm_counts_TPM_GRCh38.p13_NCBI.tsv")
cramer1_annotation <- read_tsv("1_raw_data/count_data/Cramer1/Human.GRCh38.p13.annot.tsv")
cramer1_annotation <- cramer1_annotation %>% select(GeneID, Symbol)
cramer1_w.annotation <- merge(cramer1_annotation, cramer1,by = "GeneID")
cramer1_w.annotation <- cramer1_w.annotation %>% select(Symbol,GSM1967867,GSM1967868) # chosen total RNA fragmented because TTseq used total RNA fragmented for calculating half-life

cramer1_w.annotation <- remove_duplicates_df(cramer1_w.annotation, name_column = "Symbol")
row.names(cramer1_w.annotation) <- NULL
cramer1_w.annotation <- cramer1_w.annotation %>% column_to_rownames("Symbol")
cramer1_progeny <- as.matrix(cramer1_w.annotation)
cramer1_progeny_log10_2 <- log10(cramer1_progeny+0.1)
cramer1_progeny_output1 <- progeny(cramer1_progeny_log10_2, scale = F, verbose =T, organism = "Human", top = 100)

# saveRDS(cramer1_progeny_output1, file = "4_processed_data/wholeset/progeny/cramer1/20240822_cramer1_progeny_log10_tpm.RDS")

## dieterich

dieterich <- read_tsv("1_raw_data/count_data/Dieterich/GSE49831_norm_counts_TPM.gz")
# dieterich annotation
dieterich_annotation <- read_tsv("1_raw_data/count_data/Dieterich/Human.GRCh38.p13.annot.tsv.gz")
dieterich_annotation <- dieterich_annotation %>% select(GeneID, EnsemblGeneID,Symbol)
dieterich_w.annotation <- merge(dieterich_annotation, dieterich, by = "GeneID")
dieterich_w.annotation <- dieterich_w.annotation %>% select(Symbol,GSM1269466, GSM1269469)

dieterich_w.annotation <- remove_duplicates_df(dieterich_w.annotation, name_column = "Symbol")
row.names(dieterich_w.annotation) <- NULL
dieterich_w.annotation <- dieterich_w.annotation %>% column_to_rownames("Symbol")
dieterich_progeny <- as.matrix(dieterich_w.annotation)
dieterich_progeny_log10 <- log10(dieterich_progeny+0.1)
dieterich_progeny_output <- progeny(dieterich_progeny_log10, scale =F, verbose =T, organism = "Human", top = 100)
# saveRDS(dieterich_progeny_output, file = "4_processed_data/wholeset/progeny/dieterich/20240822_dieterich_progeny_log10_tpm.RDS")

## mortavazi
mortavazi <- read_tsv("1_raw_data/count_data/Mortavazi/Mortazavi_TPM_first_col.txt")
mortavazi <- remove_duplicates_df(mortavazi, name_column = "gene")
mortavazi <- mortavazi %>% column_to_rownames(var = "gene")
mortavazi_matrix <- as.matrix(mortavazi)
mortavazi_matrix_log <- log10(mortavazi_matrix + 0.1)
mortavazi_progeny_output <- progeny(mortavazi_matrix_log, scale =F, verbose =T, organism = "Human", top = 100)

# saveRDS(mortavazi_progeny_output, file = "4_processed_data/wholeset/progeny/mortavazi/20240822_mortavazi_progeny_log10_tpm.RDS")

## rissland 
rissland <- readRDS("1_raw_data/count_data/Rissland/rissland_count_data.RDS")
rissland <- rissland %>% select(Symbol, rep2, rep3)
rissland <- remove_duplicates_df(rissland, name_column = "Symbol")
row.names(rissland) <- NULL
rissland <- column_to_rownames(rissland, var = "Symbol")
rissland_matrix <- as.matrix(rissland)
rissland_matrix_log10 <- log10(rissland_matrix + 0.1)

rissland_progeny_output <- progeny(rissland_matrix_log10, scale =F, verbose =T, organism = "Human", top = 100)
# saveRDS(rissland_progeny_output, file = "20240822_rissland_progeny_log10_tpm.RDS")

## Mulder
mulder <- read_tsv("1_raw_data/Mulder/GRCh38.p13-TPM.tsv")
colnames(mulder)[1] <- "ensembl_gene_id"
head(mulder)

mulder_biomart <- getBM(
  attributes = c("ensembl_gene_id","hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = mulder$ensembl_gene_id,
  mart = ensembl
)

mulder_w.annotation <- merge(mulder, mulder_biomart, by = "ensembl_gene_id")
head(mulder_w.annotation)
mulder_w.annotation <- mulder_w.annotation %>% select(-ensembl_gene_id)


mulder_w.annotation <- remove_duplicates_df(mulder_w.annotation, name_column = "hgnc_symbol")
row.names(mulder_w.annotation) <- NULL
mulder_w.annotation <- mulder_w.annotation %>% column_to_rownames("hgnc_symbol")
mulder_matrix <- as.matrix(mulder_w.annotation)
mulder_matrix_log10 <- log10(mulder_matrix + 0.1)
mulder_progeny_output <- progeny(mulder_matrix_log10, scale =F, verbose =T, organism = "Human", top = 100)

# saveRDS(mulder_progeny_output, file = "4_processed_data/wholeset/progeny/mulder/20240822_mulder_progeny_log10_tpm.RDS")

# check count data from here and make the plot as a QC figure 

