# from the heatmaps and scatterplots seen in: 

# "6_results/20250527_Rissland_HL_scatterplots.png"
# "6_results/20250527_heatmap_selected_4sU.studies.png"

# Rissland sample 1 and 4 will be removed from the dataset
# WHOLESET

library(here)
library(tidyverse)


hl <- readRDS(here("4_processed_data/wholeset/20250516_hl_4sU.RDS"))

# removal of rissland sample 1 and 4
hl <- hl %>% select(-any_of(c("Rissland_4sU_HEK293_1", "Rissland_4sU_HEK293_4")))

# removes genes with only NA
hl_removedNA <- hl[rowSums(is.na(hl)) != ncol(hl), ]

hl_removedNA  <- tibble::rownames_to_column(hl_removedNA, var = "ensembl_gene_id")
hl_long <- pivot_longer(hl_removedNA, cols = colnames(hl_removedNA)[-1])
# add the groupnames to it
hl_long_w.group <- hl_long %>% mutate(group = stringr::str_sub(name, end = -3L))
hl_long_mean <- hl_long_w.group %>% group_by(ensembl_gene_id,group) %>% 
  summarise(mean_hl = mean(value, na.rm = T))

hl_wide_mean <- pivot_wider(data = hl_long_mean, names_from = group, values_from = mean_hl)
saveRDS(hl_wide_mean, file = here("4_processed_data/wholeset/no_cutoff/20250530_hl_4sU_mean.per.condition.RDS"))
# here a cutoff can be made (like minimum number of values in order to get a mean value)
# NOTE: if a cutoff is made then only these selected genes can be used for making the dynamic model as it will be compared to



# choose the cutoff here
for (i in seq(1,9)){
  print(c(i,nrow(hl_wide_mean[rowSums(!is.na(hl_wide_mean[,-1])) >= i,])))
}

# chosen cutoff = 6
hl_wide_mean_cutoff <- hl_wide_mean[rowSums(!is.na(hl_wide_mean[,-1])) >= 6,]
samples <- hl_wide_mean_cutoff$ensembl_gene_id

# create mean hl value for the static model
hl_static_average <- apply(X = hl_wide_mean_cutoff[, -1], MARGIN = 1, FUN = function(x) mean(x, na.rm = TRUE))
hl_static_average_w.samples <- data.frame(samples, hl_static_average)
colnames(hl_static_average_w.samples) <- c("ensembl_gene_id", "mean_hl")

hl_static_average_w.samples
ggplot(data = hl_static_average_w.samples) + 
  geom_histogram(aes(x = mean_hl), bins = 100)
ggsave(filename = here("6_results/wholeset/20250528_cutoff6_genes.per.dataset.png"), dpi = 400, height = 5, width = 5)
saveRDS(hl_static_average_w.samples, here("4_processed_data/20250527_static_hl_4sU_cutoff6.RDS"))

# create HL dataframe for dynamic model
# this needs to to be in long format
hl_long_mean_cutoff <- hl_long_mean[hl_long_mean$ensembl_gene_id %in% samples,]

freq <- colSums(!is.na(pivot_wider(hl_long_mean_cutoff, names_from = group, values_from = mean_hl)))
freq_df <- data.frame(Dataset = names(freq), non_NA_values = freq)
freq_df[-1,]
ggplot(data = freq_df[-1,], aes(x = Dataset, y = non_NA_values)) +
  geom_bar(stat = "identity") +
  theme(
    axis.text.x.bottom =element_text(angle = 90)
  )

hl_long_mean_cutoff_clean <- na.omit(hl_long_mean_cutoff, target.colnames = mean_hl)

saveRDS(hl_long_mean_cutoff_clean, here("4_processed_data/20250527_dynamic_hl_4sU_cutoff6.RDS"))

# without a cutoff, distribution of genes
freq <- colSums(!is.na(pivot_wider(hl_long_mean, names_from = group, values_from = mean_hl)))
freq_df <- data.frame(Dataset = names(freq), non_NA_values = freq)
freq_df[-1,]
ggplot(data = freq_df[-1,], aes(x = Dataset, y = non_NA_values)) +
  geom_bar(stat = "identity") +
  theme(
    axis.text.x.bottom =element_text(angle = 90)
  )
ggsave(filename = here("6_results/wholeset/20250528_no.cutoff_genes.per.dataset.png"), dpi = 400, height = 5, width = 5)


