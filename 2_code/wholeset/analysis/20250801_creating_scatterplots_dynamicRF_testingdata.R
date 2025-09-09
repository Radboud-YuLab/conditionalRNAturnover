# create scatterplot on the performance of the dynamic model on the testing datasets

library(tidymodels)
library(here)

# retrieve testing datasets
initial_data <- readRDS(here("4_processed_data/ML_input/20250527_dynamicHL_basic.genetic.RDS"))
ml_data <- initial_data[,grep(pattern = "basic|genetic.codon|dynamic.progeny|mean_hl|ensembl_gene_id|group", x = colnames(initial_data))]
mulder_AG1478 <- ml_data %>% filter(group == "Mulder_4sU_KC.AG1478")
rissland_HEK293 <- ml_data %>% filter(group == "Rissland_4sU_HEK293")

# load dynamic model fit
dynamic_rf_fit <- readRDS(here("6_results/ML/w_cutoff_threshhold/dynamic_model/fit/20250722_dynamic_model_fit.RDS"))

# predict

preds_mulder <- mulder_AG1478 %>% select(ensembl_gene_id,mean_hl) %>% 
  mutate(predictions = predict(dynamic_rf_fit, new_data = mulder_AG1478) %>% pull(.pred))
preds_rissland <- rissland_HEK293 %>% select(ensembl_gene_id, mean_hl) %>% 
  mutate(predictions = predict(dynamic_rf_fit, new_data = rissland_HEK293) %>% pull(.pred))

ggplot(preds_mulder, aes(x = mean_hl, y = predictions)) +
  geom_point(alpha = 0.5) + 
  xlim(0,2.5) +
  ylim(0,2.5) +
  ggpmisc::stat_poly_eq(aes(label = after_stat(rr.label),), 
                        show.legend = FALSE) +
  labs(
    y = "Dynamic-RF model RNA half-life predictions\n(log10-transformed)",
    x = "RNA half-life KC_Mulder_AG1478\n(log10-transformed)"
  ) +
  theme(
    axis.title.y.left = element_text(margin = margin(r = 5, l = 5)),
    axis.text = element_text(size = 10)
  )
ggsave(file = here("6_results/ML/w_cutoff_threshhold/dynamic_model/tune_results/graphs/20250801_dynamic.rf_predictions_mulderAG1478.png"),
       dpi = 400, height = 5, width = 6)

ggplot(preds_rissland, aes(x = mean_hl, y = predictions)) +
  geom_point(alpha = 0.5) + 
  xlim(-1,4) +
  ylim(-1,4) +
  ggpmisc::stat_poly_eq(aes(label = after_stat(rr.label),), 
                        show.legend = FALSE) +
  labs(
    y = "Dynamic-RF model RNA half-life predictions\n(log10-transformed)",
    x = "RNA half-life HEK293_Rissland\n(log10-transformed)"
  ) +
  theme(
    axis.title.y.left = element_text(margin = margin(r = 5, l = 5)),
    axis.text = element_text(size = 10)
  )

ggsave(file = here("6_results/ML/w_cutoff_threshhold/dynamic_model/tune_results/graphs/20250801_dynamic.rf_predictions_risslandHEK293.png"),
       dpi = 400, height = 5, width = 6)

