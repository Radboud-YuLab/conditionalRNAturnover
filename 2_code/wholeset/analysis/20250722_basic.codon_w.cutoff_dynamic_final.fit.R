library(tidymodels)
library(here)


initial_data <- readRDS(here("4_processed_data/ML_input/20250527_dynamicHL_basic.genetic.RDS"))
ml_data <- initial_data[,grep("basic|genetic.codon|dynamic.progeny|mean_hl|ensembl_gene_id|group", x = colnames(initial_data))]
data_train <- ml_data %>%
  filter(!(group %in% c("Mulder_4sU_KC.AG1478", "Rissland_4sU_HEK293"))) # about a 70/30 split

tune_result <- readRDS(here("6_results/ML/w_cutoff_threshhold/dynamic_model/tune_results/20250716_rand.forest_tune.result_run1.RDS"))
# best configuration determined to "Preprocessor1_Model37"

best_config <- tune_result %>% collect_metrics() %>% filter(.config == "Preprocessor1_Model37", .metric == "rsq") %>% select(mtry, trees, min_n)

workflow <- readRDS("6_results/ML/w_cutoff_threshhold/dynamic_model/tune_results/workflows/20250716_rand.forest_workflow_run1.RDS")

finalized_workflow <- workflow %>% 
  finalize_workflow(parameters = best_config)

fit_on_training <- fit(finalized_workflow ,data = data_train)
saveRDS(fit_on_training,"6_results/ML/w_cutoff_threshhold/dynamic_model/fit/20250722_dynamic_model_fit.RDS") 

testing_data <- ml_data %>% filter(group %in% c("Mulder_4sU_KC.AG1478", "Rissland_4sU_HEK293"))

# predict on testing data
predictions_on_heldout <- testing_data %>% select(ensembl_gene_id, group, mean_hl) %>% 
  mutate(.pred = predict(fit_on_training, new_data = testing_data))

save.image(file = "6_results/ML/w_cutoff_threshhold/dynamic_model/fit/20250722_dynamic_model_fit.RData")