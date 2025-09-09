# Extracting a table of available features that might relate to transcript half-life
## Note: First building a table by transcript ID, but GRANDSLAM half-life data
## is annotated by ensembl gene ID. This script also creates a table where for each
## ensembl gene ID, the canonical transcript is kept. 

setwd("C:/Users/Janou Roubroeks/surfdrive/Teaching/Internships/2024_Moelong/Data")

library(biomaRt)
library(biomartr)
library(seqinr)

# using version 106 = GRCh38.p13 (same as used for GRANDSLAM)
ensembl106 <- useEnsembl(biomart = 'genes', 
                        dataset = 'hsapiens_gene_ensembl',
                        version = 106)

# command for downgrading otherwise biomart doesnt work
devtools::install_version("dbplyr", version = "2.3.4")


# Extract attributes:
# ID variables: transcript_is_canonical,external_gene_name, ensembl_gene_id,
#   ensembl_transcript_id, chromosome_name
# Sequences: 5utr and 3utr, coding
# Length: transcript_length, cds_length
# GC content: percentage_gene_gc_content, per sequence section (calculate)

attributes = listAttributes(ensembl106) # Overview of all available features
write.csv(attributes, "Available_attributes_biomaRt.csv")
# Note: you cannot extract all features at once

# Start building by transcript ID, extract all available 
features <- getBM(attributes = c("ensembl_transcript_id","transcript_is_canonical",
                                 "ensembl_gene_id", "external_gene_name", 
                                 "chromosome_name","percentage_gene_gc_content"),
                  mart = ensembl106)


# Extract sequences (too much to download at once, need to filter by transcript ID)
coding <- getBM(attributes = c("ensembl_transcript_id",
                                  "coding"),
              filters = "ensembl_transcript_id",
              values = features$ensembl_transcript_id,
              mart = ensembl106)
utr5 <- getBM(attributes = c("ensembl_transcript_id",
                               "5utr"),
                filters = "ensembl_transcript_id",
                values = features$ensembl_transcript_id,
                mart = ensembl106)
utr3 <- getBM(attributes = c("ensembl_transcript_id",
                               "3utr"),
                filters = "ensembl_transcript_id",
                values = features$ensembl_transcript_id,
                mart = ensembl106)


sequences <- merge(coding, utr5, by = "ensembl_transcript_id")
sequences <- merge(sequences, utr3, by = "ensembl_transcript_id")

# Extract feature lengths
transcr_lengths <- getBM(attributes = c("ensembl_transcript_id",
                                  "transcript_length", "cds_length"),
                 filters = "ensembl_transcript_id",
                 values = features$ensembl_transcript_id,
                 mart = ensembl106)


# NOTE: when 5'and 3'UTR data is unavailable, the length of the transcript is 
#   equal to the length of the coding sequence

# Merge by transcript ID
feature_table <- merge(features, transcr_lengths, by = "ensembl_transcript_id")
feature_table <- merge(feature_table, transcr_lengths, by = "ensembl_transcript_id")

# add 5'and 3'UTR lengths based on sequences
# set to NA when the sequence is unavailable
feature_table$UTR5_length <- ifelse(feature_table$'5utr'=="Sequence unavailable", NA,nchar(feature_table$'5utr'))
feature_table$UTR3_length <- ifelse(feature_table$'3utr'=="Sequence unavailable", NA,nchar(feature_table$'3utr'))

# Add GC content for each of the sequences
feature_table$cds_GC <- NA
feature_table$UTR5_GC <- NA
feature_table$UTR3_GC <- NA
for(i in 1:nrow(feature_table)){
  feature_table[i,"cds_GC"] <- ifelse(feature_table[i,"coding"]=="Sequence unavailable", NA,round(GC(s2c(feature_table[i,"coding"]))*100,2))
  feature_table[i,"UTR5_GC"] <- ifelse(feature_table[i,"5utr"]=="Sequence unavailable", NA,round(GC(s2c(feature_table[i,"5utr"]))*100,2))
  feature_table[i,"UTR3_GC"] <- ifelse(feature_table[i,"3utr"]=="Sequence unavailable", NA,round(GC(s2c(feature_table[i,"3utr"]))*100,2))
}

# Save full table
write.csv(feature_table, "Features_all_transcripts.csv")

# Filter to canonical genes only (for each unique ens gene ID there is 1 canonical transcript)
feature_table_genes <- feature_table[!is.na(feature_table$transcript_is_canonical),]

write.csv(feature_table_genes, "Features_genes.csv")

# Filter to hl data
feature_table_ker <- feature_table_genes[ feature_table_genes$ensembl_gene_id %in% rownames(hl), ]

write.csv(feature_table_ker, "Features_keratinocyte_HLdata.csv")

