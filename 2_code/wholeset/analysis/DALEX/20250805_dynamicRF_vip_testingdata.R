library(tidymodels)
library(DALEXtra)
library(tidyverse)
library(here)
# library(vip)
library(ggpubr)
library(stringr)


# retrieve the training data
# this loads the training data for the dynamic model, change the training data for static to make an explainer object with the static model
initial_data <- readRDS(here("4_processed_data/ML_input/20250527_dynamicHL_basic.genetic.RDS"))
ml_data <- initial_data[,grep("basic|genetic.codon|dynamic.progeny|mean_hl|ensembl_gene_id|group", x = colnames(initial_data))]
data_test <- ml_data %>%
  filter(group == "Mulder_4sU_KC.AG1478") # about a 70/30 split


rf_fit <- readRDS(here("6_results/ML/w_cutoff_threshhold/dynamic_model/fit/20250722_dynamic_model_fit.RDS"))
# importance <- extract_fit_parsnip(fitted_run)$fit
print("explaining")
set.seed(1998)
explainer_rf <- 
  explain_tidymodels(
    rf_fit, 
    data = data_test, 
    y = data_test$mean_hl,
    label = "random forest",
    verbose = TRUE
  )

set.seed(1998)
print("model parts")
vip_rf_rmse <- DALEX::model_parts(explainer_rf, loss_function = loss_root_mean_square, N = NULL)
saveRDS(vip_rf_rmse, here("6_results/ML/w_cutoff_threshhold/dynamic_model/DALEX/20250805_dynamicRF_vip_test"))

ggplot_imp <- function(..., top_n = 20, name_xlab) {
  obj <- list(...)
  metric_lab <- name_xlab
  
  full_vip <- bind_rows(obj) %>%
    filter(variable != "_baseline_")
  
  perm_vals <- full_vip %>% 
    filter(variable == "_full_model_") %>% 
    group_by(label) %>% 
    summarise(dropout_loss = mean(dropout_loss))
  
  
  
  top_vars <- full_vip %>%
    filter(variable != "_full_model_") %>%
    group_by(variable) %>%
    summarise(mean_dropout_loss = mean(dropout_loss, na.rm = TRUE)) %>%
    arrange(desc(mean_dropout_loss)) %>%
    slice_head(n = top_n) %>%
    pull(variable)
  
  p <- full_vip %>%
    filter(variable %in% top_vars, variable != "_full_model_") %>%
    mutate(variable = fct_reorder(variable, dropout_loss)) %>%
    ggplot(aes(dropout_loss, variable))
  
  if (length(obj) > 1) {
    p <- p + 
      facet_wrap(vars(label)) +
      geom_vline(data = perm_vals, aes(xintercept = dropout_loss, color = label),
                 linewidth = 1.4, lty = 2, alpha = 0.7) +
      geom_boxplot(aes(color = label, fill = label), alpha = 0.2)
  } else {
    p <- p + 
      geom_vline(data = perm_vals, aes(xintercept = dropout_loss),
                 linewidth = 1.4, lty = 2, alpha = 0.7) +
      geom_boxplot(fill = "#91CBD765", alpha = 0.4)
  }
  
  p +
    theme(legend.position = "none") +
    labs(x = metric_lab, 
         y = NULL, fill = NULL, color = NULL) +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 14)
    )
}

ggplot_imp(vip_rf_rmse, top_n = 10, name_xlab = "Root mean squared error\nafter permutations")
ggsave(filename = "6_results/ML/w_cutoff_threshhold/dynamic_model/DALEX/20250806_variable.importance.test_globaltop10.png", height = 8, width = 6, dpi = 500)

ggplot_imp(vip_rf_rmse, top_n = 20, name_xlab = "Root mean squared error\nafter permutations")
ggsave(filename = "6_results/ML/w_cutoff_threshhold/dynamic_model/DALEX/20250806_variable.importance.test_globaltop20.png", height = 8, width = 6, dpi = 500)

##--------------------------------------------------------------------------------------##


# ggplot_pdp <- function(obj, x) {
#   p <- 
#     as_tibble(obj$agr_profiles) %>%
#     mutate(`_label_` = stringr::str_remove(`_label_`, "^[^_]*_")) %>%
#     ggplot(aes(`_x_`, `_yhat_`)) +
#     geom_line(data = as_tibble(obj$cp_profiles), # remove this to leave out the grey lines from the cp profiles
#               aes(x = {{ x }}, group = `_ids_`),
#               linewidth = 0.5, alpha = 0.2, color = "gray50")
#   
#   num_colors <- n_distinct(obj$agr_profiles$`_label_`)
#   
#   if (num_colors > 1) {
#     p <- p + geom_line(aes(color = `_label_`), linewidth = 1.2, alpha = 0.8)
#   } else {
#     p <- p + geom_line(color = "midnightblue", linewidth = 1.2, alpha = 0.8)
#   }
#   
#   p
# }
# 
# set.seed(1998)
# vars <- c("dynamic.progeny.PI3K", "basic.ejc_density", "dynamic.progeny.p53", "dynamic.progeny.Estrogen", "basic.UTR5_length")
# 



