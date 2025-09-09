# svm all flavours

# run 1 considers all features, now only considering basic features
# run 2 = static.basic
# run 3 = static.basic.codon
# run 4 = static.basic.miRNA
# run 5 = static.basic.seqweaver

# run 6 = static.basic.codon.miRNA
# run 7 = static.basic.codon.seqweaver
# run 8 = static.basic.seqweaver.miRNA

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


### inputs

## identifiers
date_of_run <- "20240807" 
number_cores <- 10L
algorithm <- "svm.all"
run_number <- "4"
cv_number <- 10
path_project <- "Data/derived_data/projects/machine_learning/removed_min_learning"
seednumber <- 1998

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


# only select features and identifiers
ml_data <- initial_data[,grep("basic|genetic.miRNA|mean_HL.all|ensembl_gene_id", x = colnames(initial_data))]

ml_data <- ml_data %>% 
  select(-contains("genetic.kmer"))

# omit any columns containing NA
ml_data <- ml_data %>% select_if(~ !any(is.na(.)))

set.seed(seednumber)
data_split <- initial_split(ml_data, prop = 0.7)
data_train <- training(data_split)
data_test <- testing(data_split)

# create cross validation dataset from training data in order to tune parameters
set.seed(seednumber)
data_cv <- vfold_cv(data_train, strata = mean_HL.all, v = cv_number)


## input the recipe and model specification
# model specification
svm_linear_kernlab_spec <-
  svm_linear(cost = tune(), margin = tune()) %>%
  set_engine('kernlab') %>%
  set_mode('regression')

svm_poly_kernlab_spec <-
  svm_poly(cost = tune(), degree = tune(), scale_factor = tune(), margin = tune()) %>%
  set_engine('kernlab') %>%
  set_mode('regression')

svm_rbf_kernlab_spec <-
  svm_rbf(cost = tune(), rbf_sigma = tune(), margin = tune()) %>%
  set_engine('kernlab') %>%
  set_mode('regression')


# recipe specification 
rec <- recipe(x = data_train, formula = mean_HL.all ~ .) %>%
  update_role(ensembl_gene_id, new_role = "identifier") %>% 
  step_log(
    (contains(c("ejc", "length"))),
    base = 10,
    offset = 0.1) %>% 
  step_log(mean_HL.all,
           base = 10) %>% 
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


# workflowsets object
ml_wfset <- workflow_set(
  preproc = list(rec),
  models = list(svm_lin = svm_linear_kernlab_spec, svm_poly = svm_poly_kernlab_spec, svm_rbf = svm_rbf_kernlab_spec),
  cross = TRUE) %>% 
  option_add(param_info = svm_lin_param, id = "recipe_svm_lin") %>% 
  option_add(param_info = svm_rbf_param, id = "recipe_svm_rbf") %>% 
  option_add(param_info = svm_poly_param, id = "recipe_svm_poly")


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
set.seed(number_cores)

# hyperparamter tuning uses bayesian optimization
ctrl <- control_bayes(save_workflow = TRUE, verbose = TRUE, save_pred = TRUE, parallel_over = "resamples")
wfset_bo <-
  ml_wfset %>% workflow_map(
    fn = "tune_bayes",
    resamples = data_cv,
    metrics = metric_set(rsq, rmse, mae),
    initial = 10,
    iter = 25,
    control = ctrl
  )


saveRDS(wfset_bo, file = here(sprintf("%s/tune_results/%s_%s_tune.results_run%s.RDS",
                                      path_project,
                                      date_of_run,
                                      algorithm,
                                      run_number)))

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
ml_predictions <- wfset_bo %>% collect_predictions()
saveRDS(ml_predictions, file = here(
  sprintf("%s/tune_results/predictions/%s_%s_predictions_run%s.RDS",
          path_project,
          date_of_run,
          algorithm,
          run_number)
))

