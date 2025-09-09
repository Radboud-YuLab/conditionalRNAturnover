# this run will use the newly processed keratinocyte dataset where the min values have been removed which are most likely artifacts 
# first checking if it runs better than runs where the minimum value has not been removed

# run 1 considers all features, now only considering basic features
# run 2 = static.basic
# run 3 = static.basic.codon
# run 4 = static.basic.miRNA


library(here)
library(tidyverse)
library(tidymodels)
library(workflows)
library(tune)
library(workflowsets)
library(multilevelmod)
library(finetune)
library(ranger)
library(future)


# input
date_of_run <- "20240801"
number_cores <- 10L
algorithm <- "RF_meanHL"
run_number <- "4"
cv_number <- 10
path_project <- "Data/derived_data/projects/machine_learning/removed_min_learning"
seednumber <- 1998

if (file.exists(here(sprintf("%s/tune_results/%s_%s_tune.results_run%s.RDS",
                             path_project,
                             date_of_run,
                             algorithm,
                             run_number)))){
  stop("tuning result file already exists")
}

## take subset
take_subset = F
take_subset_value = 4000

# start data 
initial_data <- readRDS(here("Data/derived_data/projects/compiling_data/20240708_clean_kera_HL.basic_genetic.RDS"))

# taking subset
if (take_subset == T) {
  set.seed(seednumber)
  initial_data <- initial_data[sample(nrow(initial_data), take_subset_value),]
  
  algorithm <- paste(algorithm, "subset", take_subset_value, sep = ".")
}

# checking whether tune file already exists in the folder
if (file.exists(here(sprintf("%s/tune_results/%s_%s_tune.results_run%s.RDS",
                             path_project,
                             date_of_run,
                             algorithm,
                             run_number)))){
  stop("tuning result file already exists")
}

# prep initial data
ml_data <- initial_data[,grep("basic|mean_HL.all|genetic.miRNA|ensembl_gene_id", x = colnames(initial_data))]
setdiff(colnames(initial_data), colnames(ml_data)) # check which ones are missing

# remove kmers
ml_data <- ml_data %>% select(-contains("genetic.kmer"))


# splitting data 
set.seed(seednumber)
data_split <- initial_split(ml_data, prop = 0.7)
data_train <- training(data_split)
data_test <- testing(data_split)

# create cross validation dataset from training data in order to tune parameters
set.seed(seednumber)
data_cv <- vfold_cv(data_train, strata = mean_HL.all, v = cv_number)


print("data splitted into training and test")


# random forest model
rand_forest_ranger_spec <-
  rand_forest(mtry = tune(), min_n = tune(), trees = tune()) %>%
  set_engine('ranger', importance = "impurity") %>%
  set_mode("regression")


# recipe 
rec <- recipe(mean_HL.all ~ ., data = ml_data) %>%
  update_role(ensembl_gene_id, new_role = "identifier") %>% 
  step_zv(all_predictors()) %>% 
  step_corr(all_predictors()) %>% 
  step_log(
    contains(c("ejc","length")),
    base = 10,
    offset = 0.1
  ) %>%
  step_log(mean_HL.all,
           base = 10) %>% 
  step_normalize(all_numeric_predictors())


# save the prepped data 
if (!dir.exists(here(sprintf("%s/tune_results/used_data",path_project)))) {
  dir.create(here(sprintf("%s/tune_results/used_data",path_project)))
}
prepped_rec <- rec %>% prep()
juiced_rec <- juice(prepped_rec)
print(colnames(juiced_rec))
saveRDS(juiced_rec, file = here(sprintf("%s/tune_results/used_data/%s_%s_prepped.data_run%s.RDS",
                                        path_project,
                                        date_of_run,
                                        algorithm,
                                        run_number)))

# workflow
ml_wf <- workflow() %>% 
  add_recipe(rec) %>% 
  add_model(rand_forest_ranger_spec)

print("workflow specified")

saveRDS(ml_wf, file = here(sprintf("%s/tune_results/workflows/%s_%s_workflow_run%s.RDS",
                                   path_project,
                                   date_of_run,
                                   algorithm,
                                   run_number
)))

# create directory for workflows (used for fitting)
if (!dir.exists(here(sprintf("%s/tune_results/workflows",path_project)))) {
  dir.create(here(sprintf("%s/tune_results/workflows",path_project)))
}

# save workflow
saveRDS(ml_wf, file = here(
  sprintf("%s/tune_results/workflows/%s_%s_workflow_run%s",
          path_project,
          date_of_run,
          algorithm,
          run_number)))



# CHANGING HERE
ml_param <- ml_wf %>% 
  extract_parameter_set_dials() %>% 
  update(mtry = mtry(c(1, 300))) %>% 
  update(trees = trees(c(1, 3000))) %>% 
  update(min_n = min_n(c(2, 80)))


start_grid <- 
  ml_param %>%
  grid_regular(levels = 3)

print("tuning initial")
plan(multisession, workers = number_cores)
set.seed(seednumber)
ml_initial <- 
  ml_wf %>% 
  tune_grid(resamples = data_cv, 
            grid = start_grid, 
            metrics = metric_set(rsq),
            control = control_grid(save_pred = TRUE,parallel_over = "resamples"))

print("GP")  
ctrl <- control_bayes(verbose = TRUE, save_pred = TRUE, parallel_over = "resamples", save_workflow = TRUE)
set.seed(seednumber)
ml_bo <-
  ml_wf %>%
  tune_bayes(
    resamples = data_cv,
    metrics = metric_set(rsq, rmse, mae),
    initial = ml_initial,
    param_info = ml_param,
    iter = 25,
    control = ctrl
  )
saveRDS(ml_bo, file = here(sprintf("%s/tune_results/%s_%s_tune.results_run%s.RDS",
                                   path_project,
                                   date_of_run,
                                   algorithm,
                                   run_number)))

# create metrics directory
if (!dir.exists(here(sprintf("%s/tune_results/metrics",path_project)))) {
  dir.create(here(sprintf("%s/tune_results/metrics",path_project)))
}
ml_metrics <- ml_bo %>% collect_metrics()
saveRDS(ml_metrics, file = here(
  sprintf("%s/tune_results/metrics/%s_%s_metrics_run%s",
          path_project,
          date_of_run,
          algorithm,
          run_number)))


# create predictions directory
if (!dir.exists(here(sprintf("%s/tune_results/predictions",path_project)))) {
  dir.create(here(sprintf("%s/tune_results/predictions",path_project)))
}
ml_predictions <- ml_bo %>% collect_predictions()
saveRDS(ml_predictions, file = here(
  sprintf("%s/tune_results/predictions/%s_%s_predictions_run%s",
          path_project,
          date_of_run,
          algorithm,
          run_number)
))




