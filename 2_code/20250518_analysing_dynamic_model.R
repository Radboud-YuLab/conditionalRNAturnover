# script for analysing the results of the dynamic RF model


library(tidyverse)
library(tidymodels)
library(ggpubr)
library(here)


## inputs
# tune result
tune_result <- readRDS(here("6_results/ML/w_cutoff_threshhold/dynamic_model/tune_results/20250716_rand.forest_tune.result_run1.RDS"))
predictions <- tune_result %>% collect_predictions()
# used data of that run
used_data <- readRDS("6_results/ML/w_cutoff_threshhold/dynamic_model/tune_results/used_data/20250716_rand.forest_prepped.data_run1.RDS")


# looking at the autoplot
tune_result %>% autoplot()


# check what each resample means
row_resample_ident <- predictions %>% group_by(id) %>% arrange(.row) %>% slice_head(n=1) %>% select(id, .row, mean_hl)

# link the .row to which group it corresponds to in the used_data df
row_resample_ident
for (i in row_resample_ident$.row) {
  print(used_data[i,] %>% select(group, mean_hl))
}

id_mapping <- c(
  "Resample1" = "Cramer_4sU_K562",
  "Resample2" = "Bazzini_4sU_K562",
  "Resample3" = "Dieterich_4sU_MCF7",
  "Resample4" = "Mortazavi_4sU_LCL",
  "Resample5" = "Dieterich_4sU_HEK293",
  "Resample6" = "Mulder_4sU_KC.Vehicle", 
  "Resample7" = "Mulder_4sU_KC.BMP27"
)

# determine how to choose the best configuration based on these results
metrics <- tune_result %>% collect_metrics()
rsq_metrics <- metrics  %>% filter(.metric == "rsq") %>% arrange(desc(mean))
rmse_metrics <- metrics  %>% filter(.metric == "rmse") %>% arrange(mean)




# taking a higher number of mtry will make the random forest look more like decision trees which might cause issues with overfitting, better to use a lower number
# choosing the configuration

best_config_name <- "Preprocessor1_Model37" # similar to this

tune_result %>% collect_metrics() %>% filter(.config == best_config_name, .metric == "rsq")
tune_result %>% collect_metrics() %>% filter(.config == best_config_name, .metric == "rmse")
rf_preds_best <- predictions %>% filter(.config == best_config_name)

# apply the right ids
rf_preds_best <- rf_preds_best %>%
  mutate(id = id_mapping[as.character(id)])

# add ensembl_gene_id
rf_preds_with_identifiers <- rf_preds_best %>% arrange(.row) %>% mutate(ensembl_gene_id = used_data$ensembl_gene_id)
saveRDS(rf_preds_with_identifiers, here("6_results/ML/w_cutoff_threshhold/dynamic_model/fit/20250731_rand.forest_predictions_best_config.RDS"))


p <- ggplot(rf_preds_best, aes(x = mean_hl, y = .pred, colour = id, group = id)) +  # map id to both colour and group
  geom_point(alpha = 0.5) +  # Plot the points
  geom_smooth(method = "lm", se = FALSE) +  # Add regression lines for each group
  xlim(-0.5, 5) +
  ylim(-0.5, 5) +
  ggpmisc::stat_poly_eq(aes(label = after_stat(rr.label), 
                            colour = id, group = id),
                        formula = y ~ x, 
                        parse = TRUE, 
                        label.x.npc = 'right', label.y.npc = 'top', 
                        show.legend = FALSE) +
  labs(x = "RNA half-life\n(log10-transformed)",
       y = "Dynamic-RF model RNA half-life predictions\n(log10-transformed)") +
  guides(colour = guide_legend(title = "Left-out dataset")) +
  theme(
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10))
  )
print(p)

ggsave(filename = here("6_results/ML/w_cutoff_threshhold/dynamic_model/tune_results/graphs/20250717_scatterplot_performance_dynamic.RF.png"),
       plot = p, height = 5, width = 7, dpi = 500)







