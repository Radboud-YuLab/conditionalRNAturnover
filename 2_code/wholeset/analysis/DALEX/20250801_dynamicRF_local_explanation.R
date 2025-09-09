library(tidymodels)
library(DALEXtra)
library(tidyverse)
library(here)
library(vip)
library(ggpubr)
library(stringr)


# retrieve the training data
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


ggplot_pdp <- function(obj, x) {
  
  p <- 
    as_tibble(obj$agr_profiles) %>%
    mutate(`_label_` = stringr::str_remove(`_label_`, "^[^_]*_")) %>%
    ggplot(aes(`_x_`, `_yhat_`)) 
  #   geom_line(data = as_tibble(obj$cp_profiles), # remove this to leave out the grey lines from the cp profiles
  #             aes(x = {{ x }}, group = `_ids_`),
  #             linewidth = 0.5, alpha = 0.2, color = "gray50")
  # 
  num_colors <- n_distinct(obj$agr_profiles$`_label_`)
  
  if (num_colors > 1) {
    p <- p + geom_line(aes(color = `_label_`), linewidth = 1.2, alpha = 0.8)
  } else {
    p <- p + geom_line(color = "midnightblue", linewidth = 1.2, alpha = 0.8)
  }
  
  p
}

set.seed(1998)
vars <- c("dynamic.progeny.PI3K", "basic.ejc_density", "dynamic.progeny.p53", "dynamic.progeny.Estrogen", "basic.UTR5_length")

pdp_ejc <- model_profile(explainer_rf, 
                         variables = "basic.ejc_density", 
                         groups = "group")


pdp_pi3k <- model_profile(explainer_rf, 
                          variables = "dynamic.progeny.PI3K", 
                          groups = "group", 
                          processes = 1)
pdp_utr5 <- model_profile(explainer_rf, 
                          variables = "basic.UTR5_length", 
                          groups = "group")
pdp_estrogen <- model_profile(explainer_rf, 
                              variables = "dynamic.progeny.Estrogen", 
                              groups = "group")
pdp_p53 <- model_profile(explainer_rf, 
                         variables = "dynamic.progeny.p53", 
                         groups = "group")



plot_ejc <- ggplot_pdp(pdp_ejc, basic.ejc_density) +
  labs(x = "exon junction density",
       y = "predicted half-life",
       colour = NULL) + geom_point(aes(color = `_label_`)) +
  theme(
    legend.position = "none"
  )

plot_pi3k <- ggplot_pdp(pdp_pi3k, dynamic.progeny.PI3K) +
  labs(x = "PI3K",
       y = "predicted half-life",
       colour = NULL) + geom_point(aes(color = `_label_`)) +
  theme(
    legend.position = "none"
  )
plot_utr5 <- ggplot_pdp(pdp_utr5, basic.UTR5_length) + 
  labs(x = "5'UTR length",
       y = "predicted half-life",
       colour = NULL) + geom_point(aes(color = `_label_`))

plot_estrogen <- ggplot_pdp(pdp_estrogen, dynamic.progeny.Estrogen) + 
  labs(x = "Estrogen",
       y = "predicted half-life",
       colour = NULL) + geom_point(aes(color = `_label_`)) +
  theme(
    legend.position = "none"
  )

plot_p53 <- ggplot_pdp(pdp_p53, dynamic.progeny.p53) + 
  labs(x = "p53",
       y = "predicted half-life",
       colour = NULL) + geom_point(aes(color = `_label_`))  +
  theme(
    legend.position = "none"
  )

# save.image(file = here("2_code/wholeset/analysis/DALEX/20250805_localexplanations_dynamicRF.RData"))
# 
# load("data_KC_A/analysis/dynamic_best_config/local_explanation/20240923_pdp.data.RData")
# 
combined_plot <- ggarrange(plot_pi3k,plot_ejc, plot_p53, plot_estrogen, plot_utr5,
                           ncol = 3, nrow = 2, common.legend = TRUE, legend = "right")
print(combined_plot)
ggsave(filename = "6_results/ML/w_cutoff_threshhold/dynamic_model/DALEX/20250805_dynamicRF_localexplanationstop5.png", plot = combined_plot, dpi = 500, width = 12, height = 6)
