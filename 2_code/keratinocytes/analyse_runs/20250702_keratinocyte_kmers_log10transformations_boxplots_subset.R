# script for creating a boxplot of run1-3 from keratinocyte subset ML run
# this run determines how transformations of the dataset and the addition of 3'UTR kmers result in a better performance

library(tidyverse)
library(tidymodels)
library(tools)
library(here)
library(patchwork)
library(stringr)
library(scales)
library(ggplot2)
library(ggpubr)
library(ggbreak)

source(here("2_code/scripts/20240528_readmetrics_compiled.R"))


# random forest
metrics_dir <- "6_results/keratinocytes/tune_results/metrics/subset/"
metric_tables <- list()

# loop through each file in the list
for (file in list.files(metrics_dir)) {
  print(file)
  # remove the file extension from the file name
  file_name <- file_path_sans_ext(file)
  
  # split the file name into components using "_" as the delimiter
  file_split <- strsplit(file_name, "_")
  print(file_split)
  # extract date, model, and run information from the split parts
  date_metric <- file_split[[1]][1]
  model_metric <- file_split[[1]][2]
  run_metric <- file_split[[1]][4]
  
  
  # construct the full path to the CSV file
  file_path <- here(paste0(metrics_dir, file))
  print(file_path)
  # read the CSV file into a data frame
  metric_table <- readRDS(here(file_path))
  
  # print the data frame for debugging purposes
  
  # Add the metrics
  metric_table$date <- date_metric
  metric_table$model <- model_metric
  metric_table$run <- run_metric
  metric_table$file_name <- file
  
  metric_tables <- append(metric_tables, list(metric_table))
}

file_split
combined_metric_table <- bind_rows(metric_tables)

colnames(combined_metric_table)



names_runs <- unique(combined_metric_table$run)

predictions_dir <- list_files_with_exts(here("6_results/keratinocytes/tune_results/predictions/subset/"),exts = "RDS")
print(predictions_dir)
list_predictions <- list()

for (name_run in names_runs) {
  
  top_config <- combined_metric_table %>%
    filter(.metric == "rsq", run == name_run) %>%
    arrange(desc(mean)) %>%
    slice_head(n = 1) %>%
    pull(.config)
  print(top_config) # best configuration
  # get the predictions 
  prediction_RDS <- readRDS(predictions_dir[grep(pattern = name_run, predictions_dir)]) %>% 
    filter(.config == top_config) %>% 
    group_by(id) %>%
    summarise(rsq = rsq_vec(truth = mean_HL.all, estimate = .pred),
              rmse = rmse_vec(truth = mean_HL.all, estimate = .pred)) %>% 
    ungroup() %>%
    mutate(run = name_run)
  
  list_predictions[[name_run]] <- prediction_RDS
}
all_best_predictions_rf <- bind_rows(list_predictions)
all_best_predictions_rf$algorithm <- "Random Forest"


# link run numbers to used dataset
all_best_predictions_rf$description <- ifelse(
  all_best_predictions_rf$run == "run1", "+ 3'UTR-kmers\n - transformations",
  ifelse(all_best_predictions_rf$run == "run2", "- 3'UTR-kmers\n - transformations",
         ifelse(all_best_predictions_rf$run == "run3", "- 3'UTR-kmers\n + transformations", NA)
  )
)
plot_rsq <- ggplot(all_best_predictions_rf, aes(x = description, y = rsq)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2) + 
  theme(
    axis.title.x = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.text.y.right = element_blank()
  ) +
  ylab("R-squared\n(higher is better)") +
  scale_y_continuous(
    limits = c(0,0.40),
    breaks = c(0, seq(0.10, 0.40, by = 0.05))) +
  scale_y_break(c(0.001, 0.12),scales = 30)
print(plot_rsq)

ggsave(plot = plot_rsq, filename = here("6_results/keratinocytes/analysis/graphs/20250702_subset_keratinocytes_transformations_kmers_selection.png"),
       dpi = 500,width = 6, height = 4)

