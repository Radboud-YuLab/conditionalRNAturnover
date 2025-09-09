library(tidyverse)
library(here)
library(biomaRt)
library(here)
library(Biostrings)
library(conflicted)

conflict_prefer(name = "select", winner = "dplyr")
conflict_prefer(name = "filter", winner = "dplyr")


# extract samplenames
# need some geneids otherwise extraction from biomart will take too long
hl_df <- readRDS(here("4_processed_data/wholeset/no_cutoff/20250530_hl_4sU_mean.per.condition.RDS"))
samples <- hl_df$ensembl_gene_id

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
saveRDS(basic.genetic, file = here("4_processed_data/static.features/20250530_basic.genetic.RDS"))

# 4. add targetscan
# targetscan database
targetscan_predictions <- read.csv(here("data_crowdsourcing/data/raw_data/target_scan/CWCS.txt"), sep = "\t")
targetscan_predictions_with_ident <- targetscan_predictions
colnames(targetscan_predictions_with_ident) <- gsub("\\.", "_", colnames(targetscan_predictions))
colnames(targetscan_predictions_with_ident)[-1] <- paste("genetic.miRNA",colnames(targetscan_predictions_with_ident)[-1], sep = ".") 
colnames(targetscan_predictions_with_ident)[1] <- "ensembl_gene_id"
colnames(targetscan_predictions_with_ident)

targetscan_predictions_with_ident_abs <- mutate_if(targetscan_predictions_with_ident, is.numeric, abs)

basic.genetic.ts <- inner_join(basic.genetic, targetscan_predictions_with_ident, by = "ensembl_gene_id") # removes 6000 genes
colnames(basic.genetic.ts)


# 5. add seqweaver
seqweaver_3UTR <- read.csv(here("data_crowdsourcing/data/raw_data/seq_weaver/3pUTR_avg.txt"), sep = "\t")
seqweaver_5UTR <- read.csv(here("data_crowdsourcing/data/raw_data/seq_weaver/5pUTR_avg.txt"), sep = "\t")
seqweaver_orf <- read.csv(here("data_crowdsourcing/data/raw_data/seq_weaver/ORF_avg.txt"), sep = "\t")


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

basic.genetic.ts.seqw <- inner_join(basic.genetic.ts, seqweaver_3UTR_renamed, by = "ensembl_gene_id") %>% 
  inner_join(seqweaver_5UTR_renamed, by = "ensembl_gene_id") %>% 
  inner_join(seqweaver_orf_renamed, by = "ensembl_gene_id") 

saveRDS(basic.genetic.ts.seqw, file = here("4_processed_data/static.features/20250530_basic.genetic.ts.seqw.RDS"))








