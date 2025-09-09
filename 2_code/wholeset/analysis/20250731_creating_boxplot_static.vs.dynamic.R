# goal of this script is to compare the metrics from the static-RF and dynamic-RF and visualize them on boxplot


library(tidymodels)
library(here)

# 1. retrieve predictions from both static_RF and dynamic_RF

# static_RF 
static_preds <- readRDS(here("6_results/ML/w_cutoff_threshhold/static_model/fit/20250731_static_rf.fit.whole.predictions.RDS"))
static_preds <- static_preds %>% select(ensembl_gene_id, static_pred) # no need for static average

static_preds %>% head()

# dynamic_RF
dynamic_preds <- readRDS(here("6_results/ML/w_cutoff_threshhold/dynamic_model/fit/20250731_rand.forest_predictions_best_config.RDS"))
dynamic_preds <- dynamic_preds %>% select(ensembl_gene_id,.pred, id, mean_hl)
colnames(dynamic_preds)[2] <- "dynamic_pred"


# 2. combine predictions and create rsq values from mean_hl based on group
combined_preds <- left_join(dynamic_preds, static_preds, by = "ensembl_gene_id")

rsq <- combined_preds %>%
  group_by(id) %>%
  summarise(Static = rsq_vec(truth = mean_hl, estimate = static_pred),
            Dynamic = rsq_vec(truth = mean_hl, estimate = dynamic_pred))

rsq_long <- rsq %>% pivot_longer(cols = c(Static,Dynamic))

# 3. create box plot
ggplot(rsq_long, aes(x = name, y = value)) +
  geom_boxplot() + 
  geom_point(position = position_jitter(width = 0.2), aes(colour = id),size = 4, show.legend = T) +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 12),
    axis.title.y.left = element_text(size = 13, margin = margin(r= 15)),
    axis.title.x.bottom = element_text(margin = margin(t = 10))
  ) +
  labs(x = "Random Forest",
       y = "R-squared\n(higher is better)",
       color = "Left-out dataset") +
  ylim(0, 1)
ggsave(here("6_results/ML/w_cutoff_threshhold/comparison/20250731_boxplot_static.vs.dynamic.png"), dpi = 500, width = 7, height = 5)



