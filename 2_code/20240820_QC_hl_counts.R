library(tidyverse)
library(here)
library(stringr)

x <- readRDS("4_processed_data/wholeset/no_cutoff/20250516_hl_4sU.RDS")

  
  
source("data_KC_A/scripts/20240723_make_scatterplots.R")
source("data_KC_A/scripts/20240417_make_hist_from_df_v2.R")


# only mulder, dieterich and rissland have replicates

conditions_w_replicates <- unique(samples$conditions[grep(pattern = "Mulder|Rissland|Dieterich",samples$conditions)])
get_scatterplots(df = counts, combis = conditions_w_replicates, title = "", output_file = "data_KC_A/derived_data/QC/half.life_count.data/scatterplots_readcounts_replicates.png")
get_scatterplots(df = hl, combis = conditions_w_replicates, title = "", output_file = "data_KC_A/derived_data/QC/half.life_count.data/scatterplots_halflife_replicates.png")

make_hist_from_df(df = counts,pdfname = "data_KC_A/derived_data/QC/half.life_count.data/histogram_counts_individual_samples.pdf", number_of_bins = 150)
make_hist_from_df(df = hl,pdfname = "data_KC_A/derived_data/QC/half.life_count.data/histogram_halflife_individual_samples_histograms.pdf", number_of_bins = 150)

# 0's need to be removed in mortavazi
replace_min <- function(x){
  value <- min(x, na.rm = TRUE) # checks minimum value
  ifelse(value == x, yes = NA, no = x) # replaces minimum value with NA
}

hl$Mortazavi_4sU_GM12878_1<- replace_min(hl$Mortazavi_4sU_GM12878_1)
make_hist_from_df(df = hl,pdfname = "data_KC_A/derived_data/QC/half.life_count.data/histogram_halflife_individual_samples_histograms_QC.pdf", number_of_bins = 150)


save.image("data_KC_A/raw_data/compiled_data/20240822_hl_counts_selected_studies.RData")
