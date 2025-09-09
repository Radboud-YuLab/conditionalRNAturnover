library(here)
library(tidyverse)


run1_metrics <- readRDS("6_results/keratinocytes/tune_results/metrics/20250117_rand.forest.subset.1000_metrics_run1.RDS")
run2_metrics <- readRDS("6_results/keratinocytes/tune_results/metrics/20250116_all.models.subset.1000_metrics_run2.RDS")
run3_metrics <- readRDS("6_results/keratinocytes/tune_results/metrics/20250117_rand.forest.subset.1000_metrics_run3.RDS")

run1_metrics %>% select(.metric, mean, std_err, .config) %>% filter(.metric == "rsq") %>% arrange(desc(mean)) %>% slice_head(n=1)
run2_metrics %>% select(.metric, mean, std_err, .config) %>% filter(.metric == "rsq") %>% arrange(desc(mean)) %>% slice_head(n=1)
run3_metrics %>% select(.metric, mean, std_err, .config) %>% filter(.metric == "rsq") %>% arrange(desc(mean)) %>% slice_head(n=1)
