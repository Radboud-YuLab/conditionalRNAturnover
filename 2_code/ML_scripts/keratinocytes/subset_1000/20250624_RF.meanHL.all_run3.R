# The goal of this run is to see the performance of rf model on raw data (nothing transformed of the basic features) with 1000 random observations

# run 3
# all features (no kmers)
# 1000 observations
# seed 1998
# log10 transformation HL and length based features
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

date_of_run <- "20260624"
number_cores <- 5L
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
initial_data <- readRDS(here("4_processed_data/ML_input/20250624_clean_kera_HL.basic_genetic.RDS"))

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

# create cross validation dataset from training data in order to tune parameters
set.seed(seednumber)
data_cv <- vfold_cv(ml_data, strata = mean_HL.all, v = cv_number)

print("data splitted into training and test")

## input the recipe and model specification

rand_forest_ranger_spec <-
  rand_forest(mtry = tune(), min_n = tune(), trees = tune()) %>%
  set_engine('ranger',importance = "impurity") %>%
  set_mode('regression')

# recipe
rec <- recipe(x = ml_data, formula = mean_HL.all ~ .) %>%
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

# create tune_results directory
if (!dir.exists(here(sprintf("%s/tune_results",path_project)))) {
  dir.create(here(sprintf("%s/tune_results",path_project)))
}


# create used_data directory
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


rf_param <- rand_forest_ranger_spec %>% 
  extract_parameter_set_dials() %>% 
  update(mtry = mtry(c(10, 500))) %>% 
  update(trees = trees(c(1000, 6000))) %>% 
  update(min_n = min_n(c(2, 80)))


ml_wf <- workflow() %>% 
  add_recipe(rec) %>% 
  add_model(rand_forest_ranger_spec)

# create workflows directory
if (!dir.exists(here(sprintf("%s/tune_results/workflows",path_project)))) {
  dir.create(here(sprintf("%s/tune_results/workflows",path_project)))
}

saveRDS(ml_wf, file = here(
  sprintf("%s/tune_results/workflows/%s_%s_workflow_run%s.RDS",
          path_project,
          date_of_run,
          algorithm,
          run_number)))

# enable paralell processing
plan(multisession, workers = number_cores)

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

