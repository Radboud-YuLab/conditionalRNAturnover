# this run will use the newly processed keratinocyte dataset where the min values have been removed which are most likely artifacts 
# first checking if it runs better than runs where the minimum value has not been removed


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
date_of_run <- "20240722"
number_cores <- 10L
algorithm <- "RF_meanHL"
run_number <- "1"
cv_number <- 10
path_project <- "Data/derived_data/projects/machine_learning/removed_min_learning/"

if (file.exists(here(sprintf("%s/tune_results/%s_%s_tune.results_run%s.RDS",
                             path_project,
                             date_of_run,
                             algorithm,
                             run_number)))){
  stop("tuning result file already exists")
}


# prep input data
initial_data <- readRDS("Data/derived_data/projects/compiling_data/20240708_clean_kera_HL.basic_genetic.RDS")
ml_data <- initial_data[,grep("basic|genetic|mean_HL.all|ensembl_gene_id", x = colnames(initial_data))]
setdiff(colnames(initial_data), colnames(ml_data)) # check which ones are missing

# remove kmers
ml_data <- ml_data %>% select(-contains("kmer"))


# splitting data 
set.seed(1998)
data_split <- initial_split(ml_data, prop = 0.7)
data_train <- training(data_split)
data_test <- testing(data_split)

# create cross validation dataset from training data in order to tune parameters
set.seed(1998)
data_cv <- vfold_cv(data_train, strata = mean_HL.all, v = cv_number)


print("data splitted into training and test")


# random forest model
rand_forest_ranger_spec <-
  rand_forest(mtry = tune(), min_n = tune(), trees = tune()) %>%
  set_engine('ranger', importance = "impurity") %>%
  set_mode("regression")


# recipe 
base_recipe <- recipe(mean_HL.all ~ ., data = ml_data) %>%
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


# prepped_rec <- base_recipe %>% prep()
# baked_rec <- bake(prepped_rec, new_data = NULL)

# ggplot(baked_rec) +
#   geom_histogram(aes(mean_HL.all), bins = 200)
# ggsave(filename = "Data/derived_data/projects/compiling_data/graphs/20240708_log.scaled.hl_removed.min.pdf")

# workflow
ml_wf <- workflow() %>% 
  add_recipe(base_recipe) %>% 
  add_model(rand_forest_ranger_spec)

print("workflow specified")

saveRDS(ml_wf, file = here(sprintf("%s/tune_results/workflows/%s_%s_workflow_run%s.RDS",
                                   path_project,
                                   date_of_run,
                                   algorithm,
                                   run_number
)))

# CHANGING HERE
ml_param <- ml_wf %>% 
  extract_parameter_set_dials() %>% 
  update(mtry = mtry(c(10, 500))) %>% 
  update(trees = trees(c(1000, 6000))) %>% 
  update(min_n = min_n(c(2, 80)))


start_grid <- 
  ml_param %>%
  grid_regular(levels = 3)

print("tuning initial")
start_initial <- Sys.time()
plan(multisession, workers = number_cores)
set.seed(1998)
ml_initial <- 
  ml_wf %>% 
  tune_grid(resamples = data_cv, 
            grid = start_grid, 
            metrics = metric_set(rsq),
            control = control_grid(save_pred = TRUE,parallel_over = "resamples"))
end_initial <- Sys.time()

start_gp <- Sys.time()
print("GP")  
ctrl <- control_bayes(verbose = TRUE, save_pred = TRUE, parallel_over = "resamples", save_workflow = TRUE)
set.seed(1998)
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
end_gp <- Sys.time()

duration <- data.frame(
  process = c("tune initial", "bayesion optimization"),
  start_time = c(start_initial, start_gp),
  end_time = c(end_initial, end_gp)
)

log_file <- sprintf("%s_duration_%s_run%s.csv",
                    date_of_run,
                    algorithm,
                    run_number)

write.csv(duration, file = here(sprintf("%s/tune_results/%s",
                                        path_project,
                                        log_file)))



