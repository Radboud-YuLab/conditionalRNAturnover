
# The goal of this run is to see the performance of rf model on raw data (nothing transformed of the basic features) with 1000 random observations

# run 3
# all features (excluding kmers)
# 1000 observations
# seed 1998
# with transformations
# changed to only rand forest as the performance of svm models are really bad with just 1000 observations

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

date_of_run <- "20250117"
number_cores <- 10L
algorithm <- "rand.forest"
run_number <- "3"
cv_number <- 5
path_project <- "6_results/keratinocytes/"
take_subset <- TRUE
take_subset_value <- 1000
seednumber <- 1998

if (file.exists(here(sprintf("%s/tune_results/%s_%s_tune.results_run%s.RDS",
                             path_project,
                             date_of_run,
                             algorithm,
                             run_number)))){
  stop("tuning result file already exists")
}


# prep input data
initial_data <- readRDS("4_processed_data/keratinocytes/20250624_clean_kera_HL.basic_genetic.RDS")

# taking subset if take_subset is true
if (take_subset == TRUE) {
  set.seed(seednumber)
  initial_data <- initial_data[sample(nrow(initial_data), take_subset_value),]
  
  algorithm <- paste(algorithm, "subset", take_subset_value, sep = ".")
}


ml_data <- initial_data[,grep("basic|genetic|mean_HL.all|ensembl_gene_id", x = colnames(initial_data))]
setdiff(colnames(initial_data), colnames(ml_data)) # check which ones are missing

# remove kmers
ml_data <- ml_data %>% select(-contains("kmer"))

# splitting data 
set.seed(seednumber)
data_split <- initial_split(ml_data, prop = 0.7)
data_train <- training(data_split)
data_test <- testing(data_split)

# create cross validation dataset from training data in order to tune parameters
set.seed(seednumber)
data_cv <- vfold_cv(data_train, strata = mean_HL.all, v = cv_number)

print("data splitted into training and test")

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

# recipe
rec <- recipe(x = data_train, formula = mean_HL.all ~ .) %>%
  update_role(ensembl_gene_id, new_role = "identifier") %>% 
  step_log(
    contains("mean_HL.all"),
    base = 10
  ) %>% 
  step_log(
    (contains(c("ejc", "length"))),
    base = 10,
    offset = 0.1) %>%
  step_zv(all_predictors())

if (!dir.exists(here(sprintf("%s/tune_results/used_data",path_project)))) {
  dir.create(here(sprintf("%s/tune_results/used_data",path_project)))
}

# look at recipe transformed data
prepped_rec <- rec %>% prep()
juiced_rec <- juice(prepped_rec)
saveRDS(juiced_rec, file = here(sprintf("%s/tune_results/used_data/%s_%s_prepped.data_run%s.RDS",
                                        path_project,
                                        date_of_run,
                                        algorithm,
                                        run_number)))

# update the range of hyperparameters here
# svm_lin_param <- svm_linear_kernlab_spec %>% 
#   extract_parameter_set_dials() %>%
#   update(cost = cost(c(-20,0)))
# 
# 
# svm_poly_param <- svm_poly_kernlab_spec %>% 
#   extract_parameter_set_dials() %>% 
#   update(cost = cost(c(0,10)))
# 
# svm_rbf_param <- svm_rbf_kernlab_spec %>% 
#   extract_parameter_set_dials() %>%
#   update(cost = cost(c(-5, 10)))

rf_param <- rand_forest_ranger_spec %>% 
  extract_parameter_set_dials() %>% 
  update(mtry = mtry(c(10, 500))) %>% 
  update(trees = trees(c(1000, 6000))) %>% 
  update(min_n = min_n(c(2, 80)))

# kknn_param <- nearest_neighbor_kknn_spec %>% 
#   extract_parameter_set_dials() %>% 
#   update(neighbors = neighbors(c(1,40)))


# wfset object for multiple modelling testing
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

saveRDS(ml_wf, file = here(
  sprintf("%s/tune_results/workflows/%s_%s_workflow_run%s.RDS",
          path_project,
          date_of_run,
          algorithm,
          run_number)))

# enable paralell processing
plan(multisession, workers = number_cores)
set.seed(number_cores)

# hyperparamter tuning uses bayesian optimization
ctrl <- control_bayes(verbose = TRUE, save_pred = TRUE, parallel_over = "resamples")
ml_bo <-
  ml_wf %>% tune_bayes(
    resamples = data_cv,
    metrics = metric_set(rsq, rmse, mae),
    initial = 10,
    iter = 25,
    control = ctrl, 
    param_info = rf_param
  )


# create metrics directory
if (!dir.exists(here(sprintf("%s/tune_results/metrics",path_project)))) {
  dir.create(here(sprintf("%s/tune_results/metrics",path_project)))
}
ml_metrics <- ml_bo %>% collect_metrics()
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
ml_predictions <- ml_bo %>% collect_predictions(summarize = F)
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


plot_auto <- ml_bo %>% autoplot()
ggsave(filename = sprintf("%s/tune_results/graphs/%s_%s_autoplot_run%s.png",
                          path_project,
                          date_of_run,
                          algorithm,
                          run_number), plot = plot_auto, dpi = 500)



saveRDS(ml_bo, file = here(
  sprintf("%s/tune_results/%s_%s_tune.result_run%s.RDS",
          path_project,
          date_of_run,
          algorithm,
          run_number)
))

