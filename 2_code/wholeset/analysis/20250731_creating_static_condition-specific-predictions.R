# goal of this script is to predict condition-specific half-life from the static model that has been fit with the whole dataset

library(here)
library(tidyverse)
library(tidymodels)

set.seed(1998)
# 1. retrieve fit, and static dataset 
static.fit_whole <- readRDS(here("6_results/ML/w_cutoff_threshhold/static_model/fit/20250722_static_rf.fit.whole.RDS"))
initial_static_data <- readRDS(here("4_processed_data/ML_input/20250703_b.g.w.cutoff_static.4sU_ML.input.RDS"))
static_data <- initial_static_data[,grep("basic|genetic.codon|mean_hl|ensembl_gene_id", x = colnames(initial_static_data))]

# 2. predict every geneid RNA half-life using static model
static_predictions <- static_data %>% select(ensembl_gene_id, mean_hl) %>% 
  mutate(
    static_pred = predict(static.fit_whole, new_data = static_data) %>% pull(.pred)
  )

saveRDS(static_predictions, file = here("6_results/ML/w_cutoff_threshhold/static_model/fit/20250731_static_rf.fit.whole.predictions.RDS"))
# checking correlation 
ggplot(data = static_predictions, aes(x = mean_hl, y = static_pred)) + 
  geom_point()

# very high correlation as it should be as training data and testing is the same (thus introducing bias)


# 3. load dynamic data
set.seed(1998)
initial_dynamic_data <- readRDS(here("4_processed_data/ML_input/20250527_dynamicHL_basic.genetic.RDS"))

# only required the groupnames, ensembl_gene_ids and condition_specific HL values
dynamic_data <- initial_dynamic_data[,grep("mean_hl|ensembl_gene_id|group", x = colnames(initial_dynamic_data))] %>% 
  filter(!(group %in% c("Mulder_4sU_KC.AG1478", "Rissland_4sU_HEK293")))

# 4. compare the static model's predictions to condition-specific HL values in dynamic data

# check all unique values of dynamic with static
setdiff(dynamic_data$ensembl_gene_id, static_data$ensembl_gene_id) # should be 0

static_on_dynamic <- left_join(dynamic_data, static_predictions %>% select(ensembl_gene_id, static_pred), by = "ensembl_gene_id")
any(is.na(static_on_dynamic)) # should be false


# 5. create correlation plot
p <- ggplot(static_on_dynamic, aes(x = mean_hl, y = static_pred, colour = group, group = group)) +
  geom_point(alpha = 0.5) + # Plot the points
  geom_smooth(method = "lm", se = FALSE) +  # Add regression lines for each group
  xlim(-0.5, 5) +
  ylim(-0.5, 5) +
  ggpmisc::stat_poly_eq(aes(label = after_stat(rr.label), 
                            colour = group, group = group),
                        formula = y ~ x, 
                        parse = TRUE, 
                        show.legend = FALSE) + 
  labs(x = "RNA half-life\n(log10-transformed)",
       y = "Static-RF model RNA half-life predictions\n(log10-transformed)") +
  guides(colour = guide_legend(title = "Left-out dataset")) +
  theme(
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10))
  )

ggsave(filename = here("6_results/ML/w_cutoff_threshhold/static_model/tune_results/graphs/20250731_static.RF_condition-specific-predictions_scatterplot.png"),
       plot = p, height = 5, width = 7, dpi = 500)




