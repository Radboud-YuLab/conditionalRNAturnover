library(here)
library(dplyr)
#----------------------------#
# FILTERING DATASET

# only taking in ids where the coding sequence start with ATG
features_keratinocyte <- read.csv(here("1_raw_data/Mulder/Features_keratinocyte_HLdata.csv"))


ATG_start_features_kera <- features_keratinocyte %>% 
  filter(substr(coding, 1, 3) == "ATG")

saveRDS(ATG_start_features_kera, file = here("4_processed_data/keratinocytes/ATG_start_data.RDS"))
