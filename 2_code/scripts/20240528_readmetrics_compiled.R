library(here)
library(tidymodels)
library(tools)
library(tidyverse)

# function name: metrics_get(), predictions_get() and metrics_read()
# written by Moelong Yu
# date: 20240508

## USAGE ## 
# Use metrics_get() to extract all the metrics from the RDS files from the tune_results directory
# Use metrics_read() to create a dataframe that has all the metrics with identifiers (model, run_number, date, filename)

## CHANGELOG ##
# addition of predictions_get() 2024-07-22
# predictions_get(), which extracts all the predictions made when after tuning hyperparameters
# save predictions flag during tuning must be TRUE to get predictions


metrics_get <- function(filepath = "Data/derived_data/projects/machine_learning/background_jobs/tune_results/",
                        dir_to_write = "Data/derived_data/projects/machine_learning/background_jobs/tune_results/metrics/",
                        list_missing = FALSE) {
  # retrieve a list of RDS files in specified path
  list_rds_files <- list.files(path = here(filepath), pattern = "\\.RDS$")
  missing_files <- c()
  
  for (file in list_rds_files) {
    
    # get the name of the csv file from the RDS file
    mod_file_name <- gsub(x = file, pattern = "tune.results", replacement = "metrics") 
    mod_file_name <- gsub(x = mod_file_name, pattern = "RDS", replacement = "csv")
    
    # checks if there are any metric files missing from the RDS file
    if (list_missing == TRUE) {
      if (!file.exists(paste0(dir_to_write,mod_file_name))) {
        print(paste(dir_to_write,mod_file_name, sep = ""))
        missing_files <- c(missing_files, mod_file_name)
      }
    }
    
    # reads RDS file, collect metrics and write to csv to specified metric directory
    if (list_missing == FALSE) {
      print(paste("collecting metrics from:", file))
      if (file.exists(here(paste0(dir_to_write, mod_file_name)))) {
        print(paste(mod_file_name,"already exists"))
      } else {
        print(paste0(".............processing ",file,"............."))
        
        tune_results <- readRDS(here(paste0(filepath, file))) # saves tuning result in object
        tune_metrics <- tune_results %>% collect_metrics()
        
        
        write_csv(tune_metrics, here(paste0(dir_to_write, mod_file_name)))
      }}
  }
  
  # reads the missing file object and print out which csv files are missing from RDS directory
  if (list_missing == TRUE) {
    if (length(missing_files) > 0) {
      print(paste(missing_files, "does not exist"))
    } else {
      print("no metric file missing from tune results directory")
    }
  }
}

predictions_get <- function(filepath = "Data/derived_data/projects/machine_learning/background_jobs/tune_results/",
                        dir_to_write = "Data/derived_data/projects/machine_learning/background_jobs/tune_results/metrics/",
                        list_missing = FALSE) {
  # retrieve a list of RDS files in specified path
  list_rds_files <- list.files(path = here(filepath), pattern = "\\.RDS$")
  missing_files <- c()
  
  for (file in list_rds_files) {
    
    # get the name of predictions file
    mod_file_name <- gsub(x = file, pattern = "tune.results", replacement = "predictions") 
    
    
    # checks if there are any prediction files missing from the RDS file
    if (list_missing == TRUE) {
      if (!file.exists(paste0(dir_to_write,mod_file_name))) {
        print(paste(dir_to_write,mod_file_name, sep = ""))
        missing_files <- c(missing_files, mod_file_name)
      }
    }
    
    # reads RDS file, collect metrics and write to csv to specified metric directory
    if (list_missing == FALSE) {
      print(paste("collecting predictions from:", file))
      if (file.exists(here(paste0(dir_to_write, mod_file_name)))) {
        print(paste(mod_file_name,"already exists"))
      } else {
        print(paste0(".............processing ",file,"............."))
        
        tune_results <- readRDS(here(paste0(filepath, file))) # saves tuning result in object
        tune_predictions <- tune_results %>% collect_predictions()
        
        
        write_csv(tune_predictions, here(paste0(dir_to_write, mod_file_name)))
      }}
  }
  
  # reads the missing file object and print out which csv files are missing from RDS directory
  if (list_missing == TRUE) {
    if (length(missing_files) > 0) {
      print(paste(missing_files, "does not exist"))
    } else {
      print("no metric file missing from tune results directory")
    }
  }
}

p_m_get <- function(filepath = "Data/derived_data/projects/machine_learning/background_jobs/tune_results/",
                    filter = NA,
                    list_missing = FALSE) {
  # Retrieve a list of RDS files in the specified path
  list_rds_files <- list.files(path = here(filepath), pattern = "\\.RDS$")
  
  # If a filter is provided, filter the list of files
  if (!is.na(filter)) {
    list_rds_files <- list_rds_files[grep(filter, list_rds_files)]
  }
  
  # Create directories if they do not exist
  metrics_dir <- file.path(filepath, "metrics")
  predictions_dir <- file.path(filepath, "predictions")
  
  if (!dir.exists(metrics_dir)) {
    dir.create(metrics_dir, recursive = TRUE)
  }
  
  if (!dir.exists(predictions_dir)) {
    dir.create(predictions_dir, recursive = TRUE)
  }
  
  for (file in list_rds_files) {
    # Get the names of metrics and predictions files 
    pred_file_name <- gsub(x = file, pattern = "tune.results", replacement = "predictions") 
    metrics_file_name <- gsub(x = file, pattern = "tune.results", replacement = "metrics") 
    metrics_file_name <- gsub(x = metrics_file_name, pattern = "RDS", replacement = "csv")
    
    # Read the RDS file, collect metrics and predictions, and write to the specified directories
    print(paste("Collecting predictions and metrics from:", file))
    
    tune_results <- readRDS(here(filepath, file)) # Save tuning result in an object
    print("object created")
    tune_predictions <- tune_results %>% collect_predictions()
    tune_metrics <- tune_results %>% collect_metrics()
    
    print(paste0("Processing ", pred_file_name))
    saveRDS(tune_predictions, file = here(predictions_dir, pred_file_name))
    print(paste0("Processing ", metrics_file_name))
    write_csv(tune_metrics, here(metrics_dir, metrics_file_name))
  }
}

metrics_read <- function(metrics_dir = "Data/derived_data/projects/machine_learning/background_jobs/tune_results/metrics/") {
  # if no argument is given, the metric directory from the data keratinocyte project directory will be selected
  # reads the csv files collected from collect_metrics() from tidymodels package
  # only works when the filename is described correctly i.e.20240508_rf_metrics_run0
  
  metric_tables <- list()
  # loop through each file in the list
  for (file in list.files(path = here(metrics_dir), pattern = ".csv")) {
    
    # remove the file extension from the file name
    file_name <- file_path_sans_ext(file)
    
    # split the file name into components using "_" as the delimiter
    file_split <- strsplit(file_name, "_")
    
    # extract date, model, and run information from the split parts
    date_metric <- file_split[[1]][1]
    model_metric <- file_split[[1]][2]
    run_metric <- file_split[[1]][4]
    
    # construct the full path to the CSV file
    file_path <- here(paste(metrics_dir, file, sep = ""))
    
    # read the CSV file into a data frame
    metric_table <- read.csv(file_path)
    
    # print the data frame for debugging purposes
    
    # Add the metrics 
    metric_table$date <- date_metric
    metric_table$model <- model_metric
    metric_table$run <- run_metric
    metric_table$file_name <- file
    
    metric_tables <- append(metric_tables, list(metric_table))
  }
  combined_metric_table <- bind_rows(metric_tables)
  return(combined_metric_table)
}





