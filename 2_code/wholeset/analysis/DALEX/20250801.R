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
data_train <- ml_data %>%
  filter(!(group %in% c("Mulder_4sU_KC.AG1478", "Rissland_4sU_HEK293"))) # about a 70/30 split


rf_fit <- readRDS(here("6_results/ML/w_cutoff_threshhold/dynamic_model/fit/20250722_dynamic_model_fit.RDS"))
# importance <- extract_fit_parsnip(fitted_run)$fit
print("explaining")
set.seed(1998)
explainer_rf <- 
  explain_tidymodels(
    rf_fit, 
    data = data_train, 
    y = data_train$mean_hl,
    label = "random forest",
    verbose = TRUE
  )

set.seed(1998)
print("model parts")
vip_rf_rmse <- DALEX::model_parts(explainer_rf, loss_function = loss_root_mean_square, N = NULL)
saveRDS(vip_rf_rmse, here("6_results/ML/w_cutoff_threshhold/dynamic_model/fit/20250801_dynamicRF_vip"))
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



