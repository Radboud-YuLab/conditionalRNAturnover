library(here)
library(tidymodels)
library(tidyverse)

identifier <- "20250624"


# load dataframe with basic and genetic features
basic.genetic <- readRDS(here("4_processed_data/keratinocytes/basic.genetic.RDS"))


# load half-lives
load(here("1_raw_data/Mulder/Normalised_GSdata.RData"))

# rename the columns to match with feature dataframe
head(hl)
colnames(hl)[1] <- "external_gene_name"
colnames(hl)[-1] <- paste("HL",colnames(hl)[-1], sep = ".")

source("2_code/scripts/20240417_make_hist_from_df_v2.R")

# function that replaces the minimum value with NA
replace_min <- function(x){
  value <- min(x, na.rm = TRUE) # checks minimum value
  ifelse(value == x, yes = NA, no = x) # replaces minimum value with NA
}

hl <- tibble::rownames_to_column(hl,var = "ensembl_gene_id")
head(hl)
hl_idents <- hl[,c(1,2)]
hl_values <- hl[,-c(1,2)]
hl_values_removed.min  <- data.frame(lapply(hl_values,replace_min))
hl.df_rmv.min <- cbind(hl_idents, hl_values_removed.min)


# expand the data to include mean of group A,B,V and the mean of all groups 
hl_with_means <- hl.df_rmv.min %>% 
  mutate(mean_HL.A = rowMeans(select(.,c(HL.A1,HL.A2)), na.rm = TRUE)) %>% 
  mutate(mean_HL.B = rowMeans(select(.,c(HL.B1,HL.B2, HL.B3)), na.rm = TRUE)) %>% 
  mutate(mean_HL.V = rowMeans(select(.,c(HL.V1,HL.V2,HL.V3)), na.rm = TRUE)) %>% 
  mutate(mean_HL.all = rowMeans(select(.,c(HL.A1,HL.A2,
                                           HL.B1,HL.B2, HL.B3,
                                           HL.V1,HL.V2,HL.V3)),na.rm = TRUE))

# change Nan to NA
hl_with_means <- hl_with_means %>% mutate_all(~ifelse(is.nan(.), NA, .))

# check missing values
sum(is.na(hl_with_means$mean_HL.all)) # 19 missing values, due to the minimum values that have been removed

# removing those 19 genes
filt_hl.means <- hl_with_means %>% filter(!is.na(mean_HL.all))

# make rownames a column in order to inner join
basic.genetic.hl <- inner_join(basic.genetic, filt_hl.means, by = "ensembl_gene_id")

basic.genetic.hl <- basic.genetic.hl %>% 
  relocate(colnames(basic.genetic.hl)[grep("HL.",colnames(basic.genetic.hl))], .after = chromosome_name) # relocate the joined hl columns to view the half lives easier

basic.genetic.hl$ensembl_gene_id[is.na(basic.genetic.hl$basic.UTR3_length)]

# change NA to0 when UTR3 and UTR5 are missing
clean.basic.genetic.hl <- basic.genetic.hl %>%
  replace_na(list(basic.UTR5_length = 0,
                  basic.UTR3_length = 0,
                  basic.UTR3_GC = 0,
                  basic.UTR5_GC = 0))

# checking how many missing values there are
NA_counts <- colSums(is.na(basic.genetic.hl))
NA_counts <- data.frame(NA_counts)

# save to RDS file 
saveRDS(clean.basic.genetic.hl, file = here(sprintf("4_processed_data/ML_input/%s_clean_kera_HL.basic_genetic.RDS",
                                                    identifier)))

