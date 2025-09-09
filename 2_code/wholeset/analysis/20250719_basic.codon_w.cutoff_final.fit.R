# create scatterplots from the basic.codon w. threshold cutoff static models

library(here)
library(tidymodels)
library(ggpmisc)


# 1. retrieve from this run, the testing data, training data, predictions, metrics and workflow
set.seed(1998)
initial_data <- readRDS(here("4_processed_data/ML_input/20250703_b.g.w.cutoff_static.4sU_ML.input.RDS"))
ml_data <- initial_data[,grep("basic|genetic.codon|mean_hl|ensembl_gene_id", x = colnames(initial_data))]
data_split <- initial_split(ml_data, prop = 0.7)
data_train <- training(data_split)
data_test <- testing(data_split)

# 2. extract best configurations
tune_results <- readRDS(here("6_results/ML/w_cutoff_threshhold/static_model/tune_results/20250719_all.models_tune.result_run2.RDS"))

wflow_ids <- tune_results$wflow_id

# top 10 by RMSE
top10_rmse <- map_dfr(wflow_ids, ~ {
  result <- tune_results %>% extract_workflow_set_result(.x)
  show_best(result, metric = "rmse", n = 10) %>%
    mutate(wflow_id = .x)
})

# top 10 by Rsq
top10_rsq <- map_dfr(wflow_ids, ~ {
  result <- tune_results %>% extract_workflow_set_result(.x)
  show_best(result, metric = "rsq", n = 10) %>%
    mutate(wflow_id = .x)
})

# generally, when using bayesion optimization the best configuration for rsq is the same for rmse as shown here
top10_rmse %>% group_by(wflow_id) %>% slice_head(n = 1) %>% select(mean, .config, std_err) %>% print(n = 100)
top10_rsq %>% group_by(wflow_id) %>% slice_head(n = 1) %>% print(n = 100)

# 3. apply the best configuration to the workflow for each model\
fit_best_config <- function(tune_obj, name_recipe, data_t) {
  
  # get tune_result from one model
  tune_result <- tune_obj %>% extract_workflow_set_result(name_recipe)
  
  # retrieve the best configuration
  best_config <- tune_result %>% select_best()
  
  # get workflow from tune_result
  wf <- tune_obj %>% extract_workflow(name_recipe)
  
  # enter best configuration hyperparameters into the workflow
  finalized_wf <- wf %>% finalize_workflow(best_config)
  print(finalized_wf) # for control
  print("fitting")
  
  # fit the finalized wf to training data
  fit_train <- fit(finalized_wf, data_t)
  
  return(fit_train)
}

# fit on the whole training data
rf_fit <- fit_best_config(tune_obj = tune_results,
                          name_recipe = "recipe_rf",
                          data_t = data_train)
svm_lin_fit <- fit_best_config(tune_obj = tune_results,
                               name_recipe = "recipe_svm_lin",
                               data_t = data_train)
svm_poly_fit <- fit_best_config(tune_obj = tune_results,
                                name_recipe = "recipe_svm_poly",
                                data_t = data_train)
svm_rbf_fit <- fit_best_config(tune_obj = tune_results,
                               name_recipe = "recipe_svm_rbf",
                               data_t = data_train)

# fit on all of the data for predictions in condition-specific manner
rf_fit_whole <- fit_best_config(tune_obj = tune_results,
                          name_recipe = "recipe_rf",
                          data_t = ml_data)
saveRDS(rf_fit_whole, file = here("6_results/ML/w_cutoff_threshhold/static_model/fit/20250722_static_rf.fit.whole.RDS"))

# 4. predict the testing data
predictions_on_heldout <- data_test %>% select(ensembl_gene_id, mean_hl)

predictions_on_heldout <- data_test %>%
  select(ensembl_gene_id, mean_hl) %>%
  mutate(
    rand_forest = predict(rf_fit, new_data = data_test) %>% pull(.pred),
    svm_linear = predict(svm_lin_fit, new_data = data_test) %>% pull(.pred),
    svm_polynomial = predict(svm_poly_fit, new_data = data_test) %>% pull(.pred),
    svm_rbf = predict(svm_rbf_fit, new_data = data_test) %>% pull(.pred)
  )
saveRDS(predictions_on_heldout, file = here("6_results/ML/w_cutoff_threshhold/static_model/fit/20250722_predictions_on_testing_all.models.RDS"))

predictions_long <- predictions_on_heldout %>% pivot_longer(cols = c("rand_forest", "svm_linear", "svm_polynomial", "svm_rbf"),
                                        names_to = "model",
                                        values_to = "hl_predictions")



ggplot(predictions_long, aes(x = mean_hl, y = hl_predictions)) +
  geom_point(alpha = 0.3) +
  scale_y_continuous(limits = c(-0.1, 2.5),
                     breaks = seq(0, 2.5, by = 1)) +
  scale_x_continuous(limits = c(-0.1, 2.5),
                     breaks = seq(0, 2.5, by = 1)) +
  theme(aspect.ratio = 1,
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.title.x = element_text(margin = margin(t = 8)),
        axis.title.y = element_text(margin = margin(r =8, l= 5)),
        strip.text = element_text(size = 12)) +
  stat_poly_eq(aes(label = after_stat(rr.label))) +
  labs(x = "Mean RNA half-life\n(log10-transformed)",
       y = "Static models' half-life predictions\n(log10-transformed)") +
  facet_wrap(~model)

ggsave(here("6_results/ML/w_cutoff_threshhold/static_model/tune_results/graphs/20250722_scatterplots_all.models_static_run2.png"),  dpi = 400, height = 8, width = 6)

save.image(file= here("6_results/ML/w_cutoff_threshhold/static_model/fit/20250722_static_fit.RData"))





