library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# make this more intuitive
filter_rows_with_single_non_na <- function(df) {
  df %>%
    rowwise() %>%
    filter(sum(!is.na(c_across(everything()))) > 1) %>%
    ungroup()
}


get_scatterplots <- function(df, combis, title, output_file) {
  # Create an empty list to store all plots
  # Function only works if the replicates are called like this "_1" at the end
  all_plots <- list()
  
  for (condition in combis) {
    first_selection <- colnames(df)[grep(condition, x = colnames(df))]
    
    condition_length_name <- nchar(condition)
    replicate_length <- condition_length_name + 2
    
    second_selection <- c()
    for (repl in first_selection) {
      if (nchar(repl) <= replicate_length) {
        second_selection <- c(second_selection, repl)
      }
    }
    
    print(second_selection)
    subset_df <- df[, second_selection]
    
    valid_df <- filter_rows_with_single_non_na(subset_df)
    
    # Determine common axis limits for the condition
    x_limits <- range(valid_df, na.rm = TRUE)
    y_limits <- range(valid_df, na.rm = TRUE)
    
    combinations_test <- combn(colnames(subset_df), 2)
    combinations_list <- split(combinations_test, col(combinations_test))
    
    # contains the plots of 1 condition
    condition_plots <- list()
    for (i in seq_along(combinations_list)) {
      comb <- combinations_list[[i]]
      if (i == 1){
        p <- ggplot(subset_df, aes(x = !!sym(comb[1]), y = !!sym(comb[2]))) +
          geom_point(alpha = 0.5) +
          ggtitle(condition) +
          xlab(comb[1]) +
          ylab(comb[2]) +
          theme_minimal() +
          xlim(x_limits) +
          ylim(y_limits)
        # add the pearson correlation
        p <- p + stat_cor(method = "pearson",aes(label = ..r.label..),label.x.npc = "left", label.y.npc = "top", digits = 2)
        condition_plots[[i]] <- p
      } else {
        p <- ggplot(subset_df, aes(x = !!sym(comb[1]), y = !!sym(comb[2]))) +
          geom_point(alpha = 0.5) +
          ggtitle("") +
          xlab(comb[1]) +
          ylab(comb[2]) +
          theme_minimal() +
          xlim(x_limits) +
          ylim(y_limits)
        # add the pearson correlation
        p <- p + stat_cor(method = "pearson",aes(label = ..r.label..),label.x.npc = "left", label.y.npc = "top", digits = 2)
        condition_plots[[i]] <- p 
      }
    }
    #####
    print("arranging plots")

    
    # Split condition plots into groups of 3 or less
    condition_plot_groups <- split(condition_plots, ceiling(seq_along(condition_plots) / 3))
    
    for (group in condition_plot_groups) {
      if (length(group) > 1) {
        row_plot <- ggarrange(plotlist = group, ncol = length(group))
      } else {
        row_plot <- group[[1]] + theme(aspect.ratio = 1)
      }
      all_plots[[length(all_plots) + 1]] <- row_plot
    }
  }
  print("all plots made")
  ####
  
  # Combine all rows into one figure
  final_plot <- ggarrange(plotlist = all_plots, ncol = 1, nrow = length(all_plots))
  

  annotated_plot <- annotate_figure(
    final_plot,
    top = text_grob(title, face = "bold", size = 20)
  )
  # Save the combined plot as a single PDF file
  # ggsave(paste0(output_file, ".pdf"), plot = annotated_plot, width = 15, height = 5 * length(all_plots), limitsize = FALSE)
  ggsave(paste0(output_file, ".png"), plot = annotated_plot, width = 15, height = 5 * length(all_plots), limitsize = FALSE, dpi = 300, bg = "white")
  # return(annotated_plot)
}

