
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
metrics_dir <- "6_results/keratinocytes/tune_results/metrics/rf/"
metric_tables <- list()

# loop through each file in the list
for (file in list.files(metrics_dir)) {
  print(file)
  # remove the file extension from the file name
  file_name <- file_path_sans_ext(file)
  
  # split the file name into components using "_" as the delimiter
  file_split <- strsplit(file_name, "_")
  print(file_split)
  #extract date, model, and run information from the split parts
  date_metric <- file_split[[1]][1]
  model_metric <- file_split[[1]][2]
  run_metric <- file_split[[1]][4]
  
  
  # construct the full path to the CSV file
  file_path <- here(paste0(metrics_dir, file))
  print(file_path)
  # read the CSV file into a data frame
  metric_table <- readRDS(file_path)
  
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

predictions_dir <- list_files_with_exts(here("6_results/keratinocytes/tune_results/predictions/rf/"),exts = "RDS")
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

combined_metric_table %>% filter(.metric == "rsq", run == "run2") %>% arrange(desc(mean)) 
# selecting the best config by rsq 


# svm (all)

metrics_dir <- "6_results/keratinocytes/tune_results/metrics/svm/"
metric_tables <- list()
# loop through each file in the list
for (file in list.files(metrics_dir)) {
  print(file)
  # remove the file extension from the file name
  file_name <- file_path_sans_ext(file)
  
  # split the file name into components using "_" as the delimiter
  file_split <- strsplit(file_name, "_")
  print(file_split)
  #extract date, model, and run information from the split parts
  date_metric <- file_split[[1]][1]
  model_metric <- file_split[[1]][2]
  run_metric <- file_split[[1]][4]
  
  # construct the full path to the CSV file
  file_path <- here(paste0(metrics_dir, file))
  print(file_path)
  # read the CSV file into a data frame
  metric_table <- readRDS(file_path)
  
  # print the data frame for debugging purposes
  
  # Add the metrics
  metric_table$date <- date_metric
  metric_table$model <- model_metric
  metric_table$run <- run_metric
  metric_table$file_name <- file
  
  metric_tables <- append(metric_tables, list(metric_table))
}
combined_metric_table_svm <- bind_rows(metric_tables)


predictions_dir <- list_files_with_exts(here("6_results/keratinocytes/tune_results/predictions/svm/"),exts = "RDS")
predictions_dir

recipe_models <- unique(combined_metric_table_svm$wflow_id)
run_numbers <- unique(combined_metric_table_svm$run)
list_predictions <- list()
for (recipe_model in recipe_models) {
  for (run_number in run_numbers) {
    top_config <- combined_metric_table_svm %>% 
      filter(wflow_id == recipe_model,
             .metric == "rsq",
             run == run_number) %>% 
      arrange(desc(mean)) %>%
      slice_head(n = 1) %>%
      pull(.config)
    print(top_config)
    
    # open the predicitons RDS
    predictions_RDS <- readRDS(predictions_dir[grep(pattern = run_number, predictions_dir)]) %>% #specifies predictions 
      filter(.config == top_config, wflow_id == recipe_model) %>% 
      group_by(id) %>%
      summarise(rsq = rsq_vec(truth = mean_HL.all, estimate = .pred),
                rmse = rmse_vec(truth = mean_HL.all, estimate = .pred)) %>%
      ungroup() %>%
      mutate(run = run_number,
             algorithm = recipe_model)
    
    print(predictions_RDS)
    name_data <- paste(recipe_model, run_number, sep = "_")
    
    list_predictions[[name_data]] <- predictions_RDS
  }
}
all_best_predictions_svm <- bind_rows(list_predictions)

# all_best_predictions_svm and all_best_predictions_rf contains the rsq for each fold based on the best performing hyperparameter space based on rsq

# renaming
all_best_predictions_svm$algorithm <- gsub("recipe_svm_lin", "SVM linear", all_best_predictions_svm$algorithm)
all_best_predictions_svm$algorithm <- gsub("recipe_svm_rbf", "SVM radial basis function", all_best_predictions_svm$algorithm)
all_best_predictions_svm$algorithm <- gsub("recipe_svm_poly", "SVM polynomial", all_best_predictions_svm$algorithm)

# combine the performances of rf and svm 
all_best_predictions <- rbind(all_best_predictions_rf, all_best_predictions_svm)

# link run numbers to used dataset
all_best_predictions$run <- gsub("run1", "Static All", all_best_predictions$run)
all_best_predictions$run <- gsub("run2", "Basic", all_best_predictions$run)
all_best_predictions$run <- gsub("run3", "Basic.Codon", all_best_predictions$run)
all_best_predictions$run <- gsub("run4", "Basic.miRNA", all_best_predictions$run)
all_best_predictions$run <- gsub("run5", "Basic.Seqweaver", all_best_predictions$run)
all_best_predictions$run <- gsub("run6", "Basic.Codon.miRNA", all_best_predictions$run)
all_best_predictions$run <- gsub("run7", "Basic.Codon.Seqweaver", all_best_predictions$run)
all_best_predictions$run <- gsub("run8", "Basic.miRNA.Seqweaver", all_best_predictions$run)

# add levels, to create order when making the boxplot
all_best_predictions$run <- factor(all_best_predictions$run, levels = c("Basic", "Basic.Codon","Basic.miRNA", "Basic.Seqweaver","Basic.Codon.miRNA","Basic.Codon.Seqweaver","Basic.miRNA.Seqweaver", "Static All"))
head(all_best_predictions)

all_best_predictions <- all_best_predictions %>%
  mutate(algorithm_formatted = ifelse(grepl("^SVM", algorithm), paste0(algorithm, "\n"), algorithm))

saveRDS(all_best_predictions, file = here("6_results/keratinocytes/analysis/best_performance_keratinocytes_rf_svm.RDS"))




highlighted_factors <- c("Static All", "Basic.Codon")

# custom color palette
custom_palette <- c("Static All" = "#1b9e77",  # darker shade 
                    "Basic.Codon" = "#d95f02",  # darker shade 
                    "Other" = "#7570b3")    # lighter shade

# generate all other colors 
all_factors <- levels(all_best_predictions$run)
non_highlighted_factors <- setdiff(all_factors, highlighted_factors)
non_highlighted_colors <- RColorBrewer::brewer.pal(length(non_highlighted_factors), "Pastel1")

color_mapping <- c(setNames(non_highlighted_colors, non_highlighted_factors), custom_palette)

# Top plot (plot_rsq) with custom colors
plot_rsq <- ggplot(all_best_predictions, aes(x =algorithm_formatted , y = rsq, fill = run)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.85), linewidth = 0.3) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.85), size = 0.4) +
  theme(
    legend.background = element_rect(fill = alpha('white', 0.5)),
    legend.title = element_blank(),
    axis.text.x = element_blank(),   
    axis.ticks.x = element_blank(),  
    plot.margin = margin(t = 3, b = -6),
    text = element_text(size = 15),
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank()
  ) +
  labs(x = "",
       y = "R-squared\n(higher is better)") +
  scale_y_continuous(
    limits = c(0,0.45),
    breaks = c(0, seq(0.20, 0.45, by = 0.05))) +
  scale_y_break(c(0.001, 0.15), scales =16) +
  scale_fill_manual(values = color_mapping)

print(plot_rsq)

# Bottom plot (plot_rmse) with custom colors and without legend
plot_rmse <- ggplot(all_best_predictions, aes(x = algorithm, y = rmse, fill = run)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.85), linewidth = 0.3) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.85), size = 0.4) +
  scale_x_discrete(labels = c("Random Forest", "SVM\nlinear", "SVM\npolynomial", "SVM\nradial basis function")) +
  theme(
    legend.background = element_rect(fill = alpha('white', 0.5)),
    legend.title = element_blank(),         # Remove legend
    plot.margin = margin(t = -5, b = 3),
    text = element_text(size = 15),    # Reduce top and bottom margins
    axis.text.y.right = element_blank(),
    axis.ticks.y.right = element_blank()
  ) +
  labs(x = "",
       y = "Root mean squared error\n(lower is better)") +
  scale_y_continuous(
    limits = c(0,0.36),
    breaks = c(0, seq(0.20, 0.40, by = 0.02))) +
  scale_y_break(c(0.001, 0.30), scales = 20) +
  scale_fill_manual(values = color_mapping)


combined_plot <- plot_rsq / plot_rmse + 
  plot_layout(guides = "collect") & theme(legend.position = "right")

print(combined_plot)

ggsave(file = "6_results/keratinocytes/analysis/graphs/2025101_KC_boxplots_all.models_feature.selection.png",
       plot = combined_plot, width = 12, height = 8, dpi = 500)




