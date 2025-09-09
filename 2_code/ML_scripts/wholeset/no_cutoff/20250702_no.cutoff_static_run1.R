## Checking the influence on having all features or just basic.codon 

# run1

# runs on basic.codon

# running random forests, svm and KNN 

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
library(kernlab)
library(ranger)
library(kknn)


## identifiers
date_of_run <- "20250702" 
number_cores <- 10L
algorithm <- "all.models"
run_number <- "1"
cv_number <- 10
path_project <- "6_results/ML/no_cutoff_threshhold/"
set.seed(1998)
## take subset
take_subset = F
take_subset_value = 4000

# start data 
initial_data <- readRDS(here("4_processed_data/wholeset/no_cutoff/20250618_b.g.HL_ML.input.RDS"))

# taking subset if take_subset is true
if (take_subset == T) {
  initial_data <- initial_data[sample(nrow(initial_data), take_subset_value),]
  
  algorithm <- paste(algorithm, "subset", take_subset_value, sep = ".")
}

# # checking metrics folder
# if (file.exists(here(sprintf("%s/tune_results/metrics/%s_%s_metrics_run%s.RDS",
#                              path_project,
#                              date_of_run,
#                              algorithm,
#                              run_number)))){
#   stop("tuning metrics file already exists")
# }
# 
# # checking predictions folder
# if (file.exists(here(sprintf("%s/tune_results/predictions/%s_%s_predictions_run%s.RDS",
#                              path_project,
#                              date_of_run,
#                              algorithm,
#                              run_number)))){
#   stop("tuning metrics file already exists")
# }


# only select features and identifiers
ml_data <- initial_data[,grep("basic|genetic.codon|hl_mean|ensembl_gene_id", x = colnames(initial_data))]
print(setdiff(colnames(initial_data), colnames(ml_data)))
colnames(ml_data)


data_split <- initial_split(ml_data, prop = 0.7)
data_train <- training(data_split)
data_test <- testing(data_split)

# create cross validation dataset from training data in order to tune parameters
data_cv <- vfold_cv(data_train, strata = hl_mean, v = cv_number)


## input the recipe and model specification
# model specification
svm_linear_kernlab_spec <-
  svm_linear(cost = tune(), margin = tune()) %>%
  set_engine('kernlab', importance = "impurity") %>%
  set_mode('regression')

svm_poly_kernlab_spec <-
  svm_poly(cost = tune(), degree = tune(), scale_factor = tune(), margin = tune()) %>%
  set_engine('kernlab',importance = "impurity") %>%
  set_mode('regression')

svm_rbf_kernlab_spec <-
  svm_rbf(cost = tune(), rbf_sigma = tune(), margin = tune()) %>%
  set_engine('kernlab',importance = "impurity") %>%
  set_mode('regression')

nearest_neighbor_kknn_spec <-
  nearest_neighbor(neighbors = tune(), weight_func = tune(), dist_power = tune()) %>%
  set_engine('kknn',importance = "impurity") %>%
  set_mode('regression')

rand_forest_ranger_spec <-
  rand_forest(mtry = tune(), min_n = tune(), trees = tune()) %>%
  set_engine('ranger',importance = "impurity") %>%
  set_mode('regression')

# recipe specification 
rec <- recipe(x = data_train, formula = hl_mean ~ .) %>%
  update_role(ensembl_gene_id, new_role = "identifier") %>% 
  step_log(
    contains(c("ejc", "length")),
    base = 10,
    offset = 0.1) %>%
  step_log(
    contains("hl_mean"),
    base = 10) %>%
  step_zv(all_predictors()) %>%
  step_normalize(all_predictors()) # z transformation

print(rec)
# create tune_results directory
if (!dir.exists(here(sprintf("%s/tune_results/",path_project)))) {
  dir.create(here(sprintf("%s/tune_results/",path_project)))
}

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


# update the range of hyperparameters here
svm_lin_param <- svm_linear_kernlab_spec %>% 
  extract_parameter_set_dials() %>%
  update(cost = cost(c(-20,0)))

svm_poly_param <- svm_poly_kernlab_spec %>% 
  extract_parameter_set_dials() %>% 
  update(cost = cost(c(0,10)))

svm_rbf_param <- svm_rbf_kernlab_spec %>% 
  extract_parameter_set_dials() %>%
  update(cost = cost(c(-5, 10)))

rf_param <- rand_forest_ranger_spec %>% 
  extract_parameter_set_dials() %>% 
  update(mtry = mtry(c(1, 500))) %>% 
  update(trees = trees(c(500, 3000))) %>% 
  update(min_n = min_n(c(2, 80)))

kknn_param <- nearest_neighbor_kknn_spec %>% 
  extract_parameter_set_dials() %>% 
  update(neighbors = neighbors(c(1,40)))


# make workflowset object
ml_wfset <- workflow_set(
  preproc = list(rec),
  models = list(svm_lin = svm_linear_kernlab_spec,
                svm_poly = svm_poly_kernlab_spec, 
                svm_rbf = svm_rbf_kernlab_spec,
                rf = rand_forest_ranger_spec,
                kknn = nearest_neighbor_kknn_spec),cross = TRUE) %>% 
  option_add(param_info = svm_lin_param, id = "recipe_svm_lin") %>% 
  option_add(param_info = svm_rbf_param, id = "recipe_svm_rbf") %>% 
  option_add(param_info = svm_poly_param, id = "recipe_svm_poly") %>% 
  option_add(param_info = rf_param, id = "recipe_rf") %>% 
  option_add(param_info = kknn_param, id = "recipe_kknn")

# create directory for workflows (used for fitting)
if (!dir.exists(here(sprintf("%s/tune_results/workflows",path_project)))) {
  dir.create(here(sprintf("%s/tune_results/workflows",path_project)))
}
# save workflow
saveRDS(ml_wfset, file = here(
  sprintf("%s/tune_results/workflows/%s_%s_workflow_run%s.RDS",
          path_project,
          date_of_run,
          algorithm,
          run_number)))


# enable paralell processing
plan(multisession, workers = number_cores)

# hyperparamter tuning uses bayesian optimization
ctrl <- control_bayes(verbose = TRUE, save_pred = TRUE)
wfset_bo <-
  ml_wfset %>% workflow_map(
    fn = "tune_bayes",
    resamples = data_cv,
    metrics = metric_set(rsq, rmse, mae),
    initial = 10,
    iter = 25,
    control = ctrl
  )

# create metrics directory
if (!dir.exists(here(sprintf("%s/tune_results/metrics",path_project)))) {
  dir.create(here(sprintf("%s/tune_results/metrics",path_project)))
}
ml_metrics <- wfset_bo %>% collect_metrics()
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
ml_predictions <- wfset_bo %>% collect_predictions(summarize = F)
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
plot_auto <- wfset_bo %>% autoplot()
ggsave(filename = sprintf("%s/tune_results/graphs/%s_%s_autoplot_run%s.png",
                          path_project,
                          date_of_run,
                          algorithm,
                          run_number), plot = plot_auto, dpi = 500)



saveRDS(wfset_bo, file = here(
  sprintf("%s/tune_results/%s_%s_tune.result_run%s.RDS",
          path_project,
          date_of_run,
          algorithm,
          run_number)
))