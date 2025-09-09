# data processing of agarwal's data
# only selecing half-life datasets calculated using 4sU labelling methods
# data processing of Mulder data
# QC using histograms for both
# compile agarwal and mulder hl (raw)


library(here)
library(tidyverse)

a.studies <- read.csv(here("1_raw_data/Agarwal/all_HLs_human.txt.gz"), sep = "\t")
a.idents <- a.studies[,c(1,2)]

# 1. create a sample dataframe
names.replicates_studies <- colnames(a.studies)[-c(1,2)]
samples_a <- data.frame(names.replicates_studies)
colnames(samples_a) <- "sample"
samples_a$group <- stringr::str_sub(samples_a$sample, end = -3L)

# add whether a study is 4sU or not
samples_4sU <- samples_a[grep("4sU", samples_a$group),]

# add cell line column
for (i in seq(length(samples_4sU$group))){
  samples_4sU$cell.line[i] <- str_split(samples_4sU$group,pattern = "_")[[i]][3]
}

# finalized sample list
samples_4sU <- samples_4sU %>%
  mutate(cell.line_2 = ifelse(grepl("GM", cell.line), "LCL", cell.line))

# 2. check replicates

# select datasets with replicates
counts_group <- table(samples_4sU$group)
group_names_with_replicates <- names(counts_group)[counts_group > 1] 

# selection of datasets that have replicates 
group_with_replicates <- a.studies[, grep(paste(group_names_with_replicates, collapse = "|"), colnames(a.studies))]
group_with_replicates <- cbind(a.idents, group_with_replicates)
group_with_replicates <- column_to_rownames(group_with_replicates, colnames(group_with_replicates)[1])
group_with_replicates


# all 4sU studies
a.studies.4sU <- a.studies[, grep(paste(samples_4sU$group, collapse = "|"), colnames(a.studies))]
a.studies.4sU <- cbind(a.idents, a.studies.4sU)
a.studies.4sU <- column_to_rownames(a.studies.4sU, colnames(a.studies.4sU)[1])
colnames(a.studies.4sU)

# make histograms to check
source(here("2_code/scripts/20240417_make_hist_from_df_v2.R"))
make_hist_from_df(df = a.studies.4sU, pdfname = here("6_results/QC/20250516_Agarwal.4sU_histograms.pdf"))

# Simon_4sU_K562_2 shows a large peak in the histogram, however it is no replicative value as shown here
counts_4su_studies <-  as.data.frame(table(a.studies.4sU$Simon_4sU_K562_2))
# Mortazavi_4sU_GM12878_1 has a lot of 0s, which might be artifacts
counts_4su_studies <-  as.data.frame(table(a.studies.4sU$Mortazavi_4sU_GM12878_1))

# include only the datasets with readcounts available
a.studies.4sU_select <- a.studies.4sU[,grep("Bazzini|Cramer_4sU|Dieterich|Rissland_4sU_HEK293|Mortazavi",colnames(a.studies.4sU))]

# replace GM12878 with LCL (Lymphoblastoid Cell Line)
colnames(a.studies.4sU_select) <- stringr::str_replace(string = colnames(a.studies.4sU_select), pattern = "GM12878", replacement = "LCL")
colnames(a.studies.4sU_select)

# function that replaces the minumum value with NA
replace_min <- function(x){
  value <- min(x, na.rm = T)
  ifelse(value == x, yes = NA, no = x)
}

# replace repetitive 0's to NA
a.studies.4sU_select$Mortazavi_4sU_LCL_1 <- replace_min(a.studies.4sU_select$Mortazavi_4sU_LCL_1)

# create histogram to see difference 
make_hist_from_df(df = a.studies.4sU_select, pdfname = here("6_results/QC/20250516_Agarwal.4sU_histograms.pdf"))

# 3. process keratinocyte data

# load keratinocyte HL data
load(here("1_raw_data/Mulder/Normalised_GSdata.RData"))
samples$descriptive_name[c(6,7,8)] <- c("Vehicle_1", "Vehicle_2", "Vehicle_3")
samples$new_group_names <- stringr::str_sub(samples$descriptive_name, end = -3L)
samples$new_group_names

# change column names
hl_values <- hl[,-1]

kera_treatment_names <- paste("KC", samples$descriptive_name, sep = ".")
colnames(hl_values) <- paste("Mulder","4sU",kera_treatment_names, sep = "_")
samples$column_name <- paste("Mulder","4sU",kera_treatment_names, sep = "_")
samples$column_groupname <- stringr::str_sub(samples$column_name, end = -3L)

make_hist_from_df(df = hl, pdfname = here("6_results/QC/20250516_Mulder_histograms.pdf"))
# Mulder has possible artifact, multiple 0's

# removed repetitive minimum values
hl_values_removed_min <- as.data.frame(lapply(hl_values, replace_min))

# log10 transformation of HL as was done in Agarwal's dataset
hl_values_removed_min_log10 <- log10(hl_values_removed_min)

# create histograms of transformed KC data
make_hist_from_df(df = hl_values_removed_min_log10, pdfname = here("6_results/QC/20250516_Mulder_removedmin_log10_histograms.pdf"))



# add back identifiers to df
hl_removed_min_log10 <- cbind(as.data.frame(rownames(hl)), hl_values_removed_min_log10)
hl_removed_min_log10 <- tibble::column_to_rownames(hl_removed_min_log10, var = colnames(hl_removed_min_log10)[1])


# 4. compile agarwal and mulder dataset
# NOTE: two parts, one is a subset, the other is complete
# most analysis done on subset

# subset
merged_df_a.studies_mulder <- merge(hl_removed_min_log10, a.studies.4sU_select, by = 0)
merged_df_a.studies_mulder <- tibble::column_to_rownames(merged_df_a.studies_mulder, var = "Row.names")
saveRDS(merged_df_a.studies_mulder, here("4_processed_data/subset/20250516_hl_4sU_subset.RDS"))

# complete
merged_df_a.studies_mulder_complete <- merge(hl_removed_min_log10, a.studies.4sU_select, by = 0, all = TRUE)
merged_df_a.studies_mulder_complete <- tibble::column_to_rownames(merged_df_a.studies_mulder_complete, var = "Row.names")
saveRDS(merged_df_a.studies_mulder_complete, here("4_processed_data/wholeset/no_cutoff/20250516_hl_4sU.RDS"))


# check naming LCL
colnames(merged_df_a.studies_mulder)
colnames(merged_df_a.studies_mulder_complete)
# check subset is inner joined
str(merged_df_a.studies_mulder)

# check complete is outer joined
str(merged_df_a.studies_mulder_complete)
# check no repetitive minimum values in histogram





