# Date: 2024-04-17
# Author: Moelong Yu

library(here)
library(ggplot2)


# function to make histogram of numeric variables in data
make_hist_from_df <- function(df, pdfname, number_of_bins = 100, filter = NA) {
  numeric_variables <- colnames(df[, sapply(df, is.numeric)]) # get list of numeric variables in df

  list_plots = list()
  if (is.na(filter)) {
    print("no filter")
    for (i in numeric_variables) { 
      temp_plot <- ggplot(df) +
        geom_histogram(aes_string(x = i), bins = number_of_bins) + # important to use aes_string function instead of aes 
        expand_limits(x = 0) # always show number of counts of the variable at 0
      list_plots[[i]] = temp_plot
    }
  } else {
     print("filter specified")
     numeric_variables <- numeric_variables[grepl(filter, numeric_variables, fixed = TRUE)]
     for (i in numeric_variables) { 
      temp_plot <- ggplot(df) +
      geom_histogram(aes_string(x = i), bins = number_of_bins) + # important to use aes_string function instead of aes 
      expand_limits(x = 0) # always show number of counts of the variable at 0
    list_plots[[i]] = temp_plot
     }
   }
  
  
  # combine plots in one pdf
  pdf(pdfname) #specify the filename
  for (i in numeric_variables) { # iterate over the variables that contain numbers
    print(list_plots[[i]])
  }
  dev.off()
}



