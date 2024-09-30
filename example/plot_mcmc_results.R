#install("gridExtra")
library(ggplot2)
library(gridExtra)
library(grid)

# Define the plotting function
plot_psd_matrix <- function(p, psd_matrix, lower_quantile, upper_quantile, freq_vector) {
  
  # Create an empty list to store the individual plots
  plot_list <- list()
  
  # Loop through each combination of dimensions
  for (i in 1:p) {
    for (j in 1:p) {
      
      # Get the index for this subplot (row i, column j in the matrix)
      idx <- (i - 1) * p + j
      
      # Extract PSD values for this (i, j) combination
      psd_median <- psd_matrix[i, j, ]
      psd_lower <- lower_quantile[i, j, ]
      psd_upper <- upper_quantile[i, j, ]
      
      # Create the data frame for ggplot
      psd_df <- data.frame(
        Frequency = freq_vector,
        Median = psd_median,
        Lower = psd_lower,
        Upper = psd_upper
      )
      
      # Create the ggplot for this subplot
      p_plot <- ggplot(psd_df, aes(x = Frequency)) +
        geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "lightblue", alpha = 0.5) + # Fill between
        geom_line(aes(y = Median), color = "blue", linewidth = 1) +                          # Median line
        labs(title = paste("PSD for dim", i, "vs", j),
             x = "Frequency", y = "PSD") +
        theme_minimal()
      
      # Store the plot in the list
      plot_list[[idx]] <- p_plot
      
      # Save each plot immediately
      ggsave(paste0("psd_plot_dim", i, "_vs_", j, ".png"), plot = p_plot, width = 10, height = 10)
    }
  }
  
  # Arrange the plots in a grid for display
  grid_plot <- marrangeGrob(plot_list, nrow = p, ncol = p)
  
  # Return the combined plot object (if needed for further operations)
  return(grid_plot)
}

load("mcmc_results.RData")
dims <- dim(mcmc_vnp_avg$data)  
npts <- dims[1]                 
p <- dims[2] 
psd_matrix <- mcmc_vnp_avg$psd.median 
lower_quantile <- mcmc_vnp_avg$psd.u05
upper_quantile <- mcmc_vnp_avg$psd.u95
freq_vector <- seq(0, 0.5, length.out = dim(mcmc_vnp_avg$psd.u05)[3])
myplt <- plot_psd_matrix(p, psd_matrix, lower_quantile, upper_quantile, freq_vector)


ggsave("psd_matrix_plot.png", plot = myplt, width = 10, height = 10)


