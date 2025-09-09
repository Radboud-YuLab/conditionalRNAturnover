# Creating figure where the 4sU studies half-life are represented in a heatmap and scatterplots showing why replicate 1 and 4
# from the Rissland study should be left out of the dataset
# WHOLESET


library(tidymodels)
library(tidyverse)
library(ComplexHeatmap)
library(ggplot2)
library(gplots)
library(circlize)
library(ggpubr)

source("data_KC_A/scripts/20240723_make_scatterplots.R")

merged_df_a.studies_mulder <- readRDS("4_processed_data/20250516_hl_4sU.RDS")

colnames(merged_df_a.studies_mulder)
# rename Mortazavi
cor.df_selection <- cor(merged_df_a.studies_mulder,use = "pairwise.complete.obs")

col_fun = colorRamp2(c(0, 0.5, 1), c("darkblue", "whitesmoke","indianred"))
png("6_results/20250527_heatmap_selected_4sU.studies.png", res = 300, width = 8, height = 7, units = "in")

ComplexHeatmap::Heatmap(cor.df_selection,col = col_fun,
                        name = "Pearson\ncorrelation",
                        row_names_gp =gpar(fontsize = 9),
                        column_names_gp = gpar(fontsize = 9),
                        width= unit(10, "cm"),
                        height = unit(10, "cm"))
dev.off()


rissland_rep_combinations <- combn(colnames(merged_df_selection)[grep(pattern = "Rissland_4sU_HEK293", x = colnames(merged_df_selection))], m =2)
rissland_rep_combinations <- split(rissland_rep_combinations, col(rissland_rep_combinations))
condition_plots <- list()
for (i in seq_along(rissland_rep_combinations)) {
  comb <- rissland_rep_combinations[[i]]
  print(comb)
  p <- ggplot(merged_df_selection, aes(x = !!sym(comb[1]), y = !!sym(comb[2]))) +
    geom_point(alpha = 0.5) +
    ggtitle("") +
    xlab(comb[1]) +
    ylab(comb[2]) +
    theme_minimal() +
    xlim(-1,5) +
    ylim(-1,5) +
    theme(
      aspect.ratio = 1
    )
  p <- p + stat_cor(method = "pearson",aes(label = ..r.label..),label.x.npc = "left", label.y.npc = "top", digits = 2)
  condition_plots[[i]] <- p
}

combined_plot <- ggarrange(plotlist = condition_plots, ncol = 3, nrow = 2)
print(combined_plot)

ggsave("6_results/20250527_Rissland_HL_scatterplots.png", plot = combined_plot, dpi = 500, bg = "white", width = 10, height = 7)




