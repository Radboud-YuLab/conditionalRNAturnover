# script for showing the difference 

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

# 20250702 CHANGE HERE
run1_metrics <- readRDS("6_results/ML/no_cutoff_threshhold/tune_results/metrics/20250618_all.models_metrics_run1.RDS")
run2_metrics <- readRDS("6_results/ML/no_cutoff_threshhold/tune_results/metrics/20250618_all.models_metrics_run2.RDS")

run1_metrics$description <- "Basic.Codon"
run2_metrics$description <- "Static all"

run1_metrics %>%
  select(model, .metric, mean, std_err, description, .config) %>%
  filter(.metric == "rsq") %>%
  group_by(model) %>%
  arrange(desc(mean)) %>%
  slice_head(n=1)

run2_metrics %>%
  select(model, .metric, mean, std_err, description, .config) %>%
  filter(.metric == "rsq") %>%
  group_by(model) %>%
  arrange(desc(mean)) %>%
  slice_head(n=1)
# here it shows that run1 with only basic.codon is performing as good or even better 

metrics_dir <- "6_results/ML/no_cutoff_threshhold/tune_results/metrics/"
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
  # read the RDS file into a data frame
  metric_table <- readRDS(file_path)
  
  # add the metrics
  metric_table$date <- date_metric
  metric_table$model <- model_metric
  metric_table$run <- run_metric
  metric_table$file_name <- file
  
  metric_tables <- append(metric_tables, list(metric_table))
}

file_split # should be date, models, metrics and run
combined_metric_table <- bind_rows(metric_tables)

colnames(combined_metric_table)

names_runs <- unique(combined_metric_table$run)
predictions_dir <- list_files_with_exts(here("6_results/ML/no_cutoff_threshhold/tune_results/predictions/"),exts = "RDS")

recipe_models <- unique(combined_metric_table$wflow_id)
print(recipe_models)
run_numbers <- unique(combined_metric_table$run)
print(run_numbers)
list_predictions <- list()
for (recipe_model in recipe_models) {
  for (run_number in run_numbers) {
    top_config <- combined_metric_table %>% 
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
      summarise(rsq = rsq_vec(truth = hl_mean, estimate = .pred),
                rmse = rmse_vec(truth = hl_mean, estimate = .pred)) %>%
      ungroup() %>%
      mutate(run = run_number,
             algorithm = recipe_model)
    
    print(predictions_RDS)
    name_data <- paste(recipe_model, run_number, sep = "_")
    
    list_predictions[[name_data]] <- predictions_RDS
  }
}
all_best_predictions <- bind_rows(list_predictions)
best_metrics <- all_best_predictions
best_metrics$algorithm <- stringr::str_replace(string = best_metrics$algorithm, pattern = "recipe_kknn", replacement = "Nearest Neighbor")
best_metrics$algorithm <- stringr::str_replace(string = best_metrics$algorithm, pattern = "recipe_rf", replacement = "Random Forest")
best_metrics$algorithm <- stringr::str_replace(string = best_metrics$algorithm, pattern = "recipe_svm_lin", replacement = "SVM\nlinear")
best_metrics$algorithm <- stringr::str_replace(string = best_metrics$algorithm, pattern = "recipe_svm_poly", replacement = "SVM\npolynomial")
best_metrics$algorithm <- stringr::str_replace(string = best_metrics$algorithm, pattern = "recipe_svm_rbf", replacement = "SVM\nradial basis function")

best_metrics$description <- ifelse(
  best_metrics$run == "run1",
  "Basic.Codon",
  "Static All"
)

best_metrics <- best_metrics %>% filter(algorithm != "Nearest Neighbor") # nearest neighbors not included due to bad performance
custom_palette <- c("Static All" = "#1b9e77",  
                    "Basic.Codon" = "#d95f02")

plot_rsq <- ggplot(best_metrics, aes(x = algorithm, y = rsq, fill = description)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.85), linewidth = 0.3) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.85), size = 1) +
  theme(
    legend.background = element_rect(fill = alpha('white', 1)),
    legend.title = element_blank(),
    legend.key.size = unit(1,units = "cm"),
    axis.text.x = element_blank(),   
    axis.ticks.x = element_blank(),  
    plot.margin = margin(t = 3, b = -10),
    axis.text.y.left = element_text(size = 10),
    axis.title.y.left = element_text(margin = margin(r =5, l=5), size =13),
    axis.text.y.right = element_blank(),     # Remove right y-axis text
    axis.ticks.y.right = element_blank(),    # Remove right y-axis ticks
    axis.line.y.right = element_blank()      # Remove right y-axis line
  ) +
  labs(x = "",
       y = "R-squared\n(higher is better)") +
  scale_y_continuous(
    limits = c(0,0.35),
    breaks = c(0, seq(0.20, 0.40,by = 0.05))) +
  scale_y_break(c(0.001, 0.20), scales =16) +
  scale_fill_manual(values = custom_palette)

print(plot_rsq)

plot_rmse <- ggplot(best_metrics, aes(x = algorithm, y = rmse, fill = description)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(0.85), linewidth = 0.3) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.85), size = 1) +
  scale_y_continuous(
    limits = c(0,0.475),
    breaks = c(0, seq(0.35, 0.475, by = 0.05))) +
  scale_y_break(c(0.001, 0.35), scales =16) +
  labs(x = "",
       y = "Root mean squared error\n(lower is better)") +
  scale_fill_manual(values = custom_palette) +
  theme(
    legend.background = element_rect(fill = alpha('white', 1)),
    legend.title = element_blank(),
    legend.key.size = unit(1,units = "cm"),
    plot.margin = margin(t = -10, b = 3),
    axis.text.y.left = element_text(size = 10),
    axis.title.y.left  = element_text(margin = margin(r =5, l=5), size = 13),
    axis.text.x = element_text(size = 10),
    axis.text.y.right = element_blank(),     # Remove right y-axis text
    axis.ticks.y.right = element_blank(),    # Remove right y-axis ticks
    axis.line.y.right = element_blank()      # Remove right y-axis line
  )
print(plot_rmse)
combined_plot <- plot_rsq / plot_rmse + 
  plot_layout(guides = "collect") & theme(legend.position = "right")

print(combined_plot)

ggsave(filename = here("6_results/ML/no_cutoff_threshhold/tune_results/graphs/20250702_boxplots_comparison_static_run1vs2_no_KNN.png"),
                       plot = combined_plot, 
                       dpi = 600, 
                       width = 8, 
                       height = 7)

