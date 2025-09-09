# dynamic model of KC_A 

# random_forest 
# run1

# using normal grid search on different initial data, now using raw progeny scores
# runs on basic and genetic.codon

# running random forest only

# libaries
library(here)
library(tidyverse)
library(tidymodels)
library(workflows)
library(tune)
library(workflowsets)
library(multilevelmod)
library(finetune)
library(future)
# depending on machine learning model
library(ranger)



## identifiers
date_of_run <- "20250716" 
number_cores <- 7L
algorithm <- "rand.forest"
run_number <- "1"
cv_number <- 10
path_project <- "6_results/ML/w_cutoff_threshhold/dynamic_model/"
seednumber <- 1998

# create results directory
if (!dir.exists(here(sprintf("%s/tune_results/",path_project)))) {
  dir.create(here(sprintf("%s/tune_results/",path_project)))
}

## take subset
take_subset = F
take_subset_value = 100

# start data 
initial_data <- readRDS(here("4_processed_data/ML_input/20250527_dynamicHL_basic.genetic.RDS"))

# taking subset if take_subset is true
if (take_subset == T) {
  set.seed(seednumber)
  initial_data <- initial_data[sample(nrow(initial_data), take_subset_value),]
  
  algorithm <- paste(algorithm, "subset", take_subset_value, sep = ".")
}

# checking metrics folder
if (file.exists(here(sprintf("%s/tune_results/metrics/%s_%s_metrics_run%s.RDS",
                             path_project,
                             date_of_run,
                             algorithm,
                             run_number)))){
  stop("tuning metrics file already exists")
}

# checking predictions folder
if (file.exists(here(sprintf("%s/tune_results/predictions/%s_%s_predictions_run%s.RDS",
                             path_project,
                             date_of_run,
                             algorithm,
                             run_number)))){
  stop("tuning metrics file already exists")
}

# only select features and identifiers
ml_data <- initial_data[,grep("basic|genetic.codon|dynamic.progeny|mean_hl|ensembl_gene_id|group", x = colnames(initial_data))]
print(setdiff(colnames(initial_data), colnames(ml_data)))

data_test <- ml_data %>% 
  filter(group %in% c("Mulder_4sU_KC.AG1478", "Rissland_4sU_HEK293"))
data_train <- ml_data %>%
  filter(!(group %in% c("Mulder_4sU_KC.AG1478", "Rissland_4sU_HEK293"))) # about a 70/30 split


# create leave one out cv
set.seed(seednumber)
data_cv <- group_vfold_cv(data = data_train, group = "group", v = length(unique(data_train$group)))
print(data_cv)
## input the recipe and model specification
# model specification
# svm_linear_kernlab_spec <-
#   svm_linear(cost = tune(), margin = tune()) %>%
#   set_engine('kernlab', importance = "impurity") %>%
#   set_mode('regression')
# 
# svm_poly_kernlab_spec <-
#   svm_poly(cost = tune(), degree = tune(), scale_factor = tune(), margin = tune()) %>%
#   set_engine('kernlab',importance = "impurity") %>%
#   set_mode('regression')
# 
# svm_rbf_kernlab_spec <-
#   svm_rbf(cost = tune(), rbf_sigma = tune(), margin = tune()) %>%
#   set_engine('kernlab',importance = "impurity") %>%
#   set_mode('regression')
# 
# nearest_neighbor_kknn_spec <-
#   nearest_neighbor(neighbors = tune(), weight_func = tune(), dist_power = tune()) %>%
#   set_engine('kknn',importance = "impurity") %>%
#   set_mode('regression')

rand_forest_ranger_spec <-
  rand_forest(mtry = tune(), min_n = tune(), trees = tune()) %>%
  set_engine('ranger',importance = "impurity") %>%
  set_mode('regression')

# recipe specification 
rec <- recipe(x = data_train, formula = mean_hl ~ .) %>%
  update_role(ensembl_gene_id, new_role = "identifier") %>%
  update_role(group, new_role = "identifier") %>% 
  step_log(
    (contains(c("ejc", "length"))),
    base = 10,
    offset = 0.1) %>%
  step_zv(all_predictors()) %>%
  step_normalize(all_predictors()) # z transformation


# save the data 
if (!dir.exists(here(sprintf("%s/tune_results/used_data",path_project)))) {
  dir.create(here(sprintf("%s/tune_results/used_data",path_project)))
}

prepped_rec <- rec %>% prep()
juiced_rec <- juice(prepped_rec)

saveRDS(juiced_rec, file = here(sprintf("%s/tune_results/used_data/%s_%s_prepped.data_run%s.RDS",
                                        path_project,
                                        date_of_run,
                                        algorithm,
                                        run_number)))

print(colnames(juiced_rec))


rf_tune_grid <- expand_grid(mtry = round(seq(20,60, length.out = 20)),
                            trees = c(500, 1000, 1500),
                            min_n = 1)
# make workflowset object
# ml_wfset <- workflow_set(
#   preproc = list(rec),
#   models = list(svm_lin = svm_linear_kernlab_spec,
#                 svm_poly = svm_poly_kernlab_spec, 
#                 svm_rbf = svm_rbf_kernlab_spec,
#                 rf = rand_forest_ranger_spec,
#                 kknn = nearest_neighbor_kknn_spec),cross = TRUE) %>% 
#   option_add(param_info = svm_lin_param, id = "recipe_svm_lin") %>% 
#   option_add(param_info = svm_rbf_param, id = "recipe_svm_rbf") %>% 
#   option_add(param_info = svm_poly_param, id = "recipe_svm_poly") %>% 
#   option_add(param_info = rf_param, id = "recipe_rf") %>% 
#   option_add(param_info = kknn_param, id = "recipe_kknn")

ml_wf <- workflow() %>% 
  add_recipe(rec) %>% 
  add_model(rand_forest_ranger_spec)

head(ml_wf)

# create directory for workflows (used for fitting)
if (!dir.exists(here(sprintf("%s/tune_results/workflows",path_project)))) {
  dir.create(here(sprintf("%s/tune_results/workflows",path_project)))
}

# save workflow
saveRDS(ml_wf, file = here(
  sprintf("%s/tune_results/workflows/%s_%s_workflow_run%s.RDS",
          path_project,
          date_of_run,
          algorithm,
          run_number)))


plan(multisession, workers = number_cores)
# hyperparamter tuning using normal grid search


print("tuning")
ml_tg <- ml_wf %>% tune_grid(
  resamples = data_cv,
  grid = rf_tune_grid,
  metrics = metric_set(rsq, rmse),
  control = control_grid(save_pred = TRUE, parallel_over = "resamples", verbose = TRUE)
)

print("tune complete")
saveRDS(ml_tg, file = here(
  sprintf("%s/tune_results/%s_%s_tune.result_run%s.RDS",
          path_project,
          date_of_run,
          algorithm,
          run_number)
))

# create metrics directory
if (!dir.exists(here(sprintf("%s/tune_results/metrics",path_project)))) {
  dir.create(here(sprintf("%s/tune_results/metrics",path_project)))
}
ml_metrics <- ml_tg %>% collect_metrics()
saveRDS(ml_metrics, file = here(
  sprintf("%s/tune_results/metrics/%s_%s_metrics_run%s.RDS",
          path_project,
          date_of_run,
          algorithm,
          run_number)))


# create predictions directory
if (!dir.exists(here(sprintf("%s/tune_results/predictions",path_project)))) {
  dir.create(here(sprintf("%s/tune_results/predictions",path_project)))
}
ml_predictions <- ml_tg %>% collect_predictions(summarize = F)
saveRDS(ml_predictions, file = here(
  sprintf("%s/tune_results/predictions/%s_%s_predictions_run%s.RDS",
          path_project,
          date_of_run,
          algorithm,
          run_number)
))

# create graphs directory
if (!dir.exists(here(sprintf("%s/tune_results/graphs",path_project)))) {
  dir.create(here(sprintf("%s/tune_results/graphs",path_project)))
}


plot_auto <- ml_tg %>% autoplot()
ggsave(filename = sprintf("%s/tune_results/graphs/%s_%s_autoplot_run%s.png",
                          path_project,
                          date_of_run,
                          algorithm,
                          run_number), plot = plot_auto, dpi = 500)
