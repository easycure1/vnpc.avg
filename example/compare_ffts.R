library(vnpc.avg)
library(reticulate)
library(ggplot2)
library(gridExtra)
library(grid)
use_python("/Users/avaj0001/Documents/projects/deep_gw_pe_followup/venv/bin/python")
py_fft_fuction <- import("sgvb_psd.backend.analysis_data")$compute_chunked_fft
set.seed(0) 


n <- 2048
d <- 2
n_segments <- as.integer(4)

simulate_data <- function(n, d){
  param_ar <- rbind(c(.1, 0, 0, 0), c(0, -.3, 0, -.5))
  param_ma <- matrix(nrow=2, ncol=0)
  param_sigma <- matrix(data=c(1, .9, .9, 1), nrow=2, ncol=2)
  data <- beyondWhittle::sim_varma(
    model=list(ar=param_ar, ma=param_ma, sigma=param_sigma), 
    n=n, d=d
    )
  return(data)
}
  

omegaFreq <- function(n) {
  return(2 / n * (1:(n / 2 + 1) - 1))
}

get_noise <- function(data, theta, ...) {
  # mean centered version
  if (class(data)[1] == "matrix" || "array") {
    apply(data, 2, center)
  } else {
    center(data)
  }
}


get_r_fft <- function(data, n_segments){
  seg_n <- n_segments
  d <- ncol(data)
  n <- nrow(data) / seg_n
  FZ <- n ## Store the original full length of each segment
  omega <- omegaFreq(n)
  N <- length(omega)
  data_new <- array(rep(NA, n*d*seg_n), dim = c(n, d, seg_n))
  for (ii in 1:seg_n) {
    data_new[,,ii] <- data[(n*(ii-1)+1):(n*ii),]
  }
  # theta_dim <- model_params$theta_dim
  theta_dim <- 0
  #Ntotal <- mcmc_params$Ntotal
  Ntotal <- 1000
  theta <- matrix(NA, nrow=theta_dim, ncol=Ntotal)
  noise <- array(rep(NA, n_init*d*seg_n), dim = c(n_init, d, seg_n))
  FZ <- array(rep(NA, N*d*seg_n), dim = c(N, d, seg_n))
  for (ii in 1:seg_n) {
    noise[,,ii] <- get_noise(data_new[,,ii], theta[,i])
    FZ[,,ii] <- vnpc.avg:::mdft(noise[,,ii])[1:N,]
  }
  return(FZ)
}



get_py_fft <- function(data, n_segments){
  x = (data - mean(data))/sd(data)
  py_fft = py_fft_fuction(
    x=r_to_py(x), 
    nchunks=n_segments, 
    fmax_for_analysis = Inf,
    fs=0.5
  )
  py_fft_data = py_fft[[1]]
  reshaped_data = aperm(py_fft_data, c(2, 3, 1))
  return (reshaped_data)
}


data = simulate_data(n, d)
r_fft = get_r_fft(data, n_segments)
py_fft = get_py_fft(data, n_segments)


# Function to create the plots using ggplot2
plot_data_gg <- function(data_array1, data_array2, color1 = "blue", color2 = "red", label1 = NULL, label2 = NULL, title="Title") {
  
  # Create a list to hold all plots
  plot_list <- list()
  
  shape <- dim(data_array1)
  print(dim(data_array1))
  print(dim(data_array2))
  
  for (i in 1:shape[3]) {
    # Get the index for the x-axis
    index1 <- 1:nrow(data_array1[, , i])
    index2 <- 1:nrow(data_array2[, , i])
    
    # Create a data frame for plotting
    df1 <- data.frame(Index = index1, Value1 = data_array1[, 1, i], Value2 = data_array1[, 2, i])
    df2 <- data.frame(Index = index2, Value1 = data_array2[, 1, i], Value2 = data_array2[, 2, i])
    
    # Create the plot for the current segment
    p <- ggplot() +
      geom_line(data = df1, aes(x = Index, y = Value1, color = color1), linetype = 1) +
      geom_line(data = df1, aes(x = Index, y = Value2, color = color1), linetype = 2) +
      geom_line(data = df2, aes(x = Index, y = Value1, color = color2), linetype = 1) +
      geom_line(data = df2, aes(x = Index, y = Value2, color = color2), linetype = 2) +
      labs(title = paste("Plot for Chunk ", i), x = "Index", y = "Value") +
      theme_minimal() +
      theme(legend.position = "top") +
      scale_color_manual(values = c(color1, color2), 
                         labels = c(label1, label2))
    
    # Add the plot to the list
    plot_list[[i]] <- p
  }
  
  # Combine the plots into a grid
  grid_plot <- grid.arrange(grobs = plot_list, ncol = 2)
  
  # Add a supertitle
  grid.text(title, x = 0.5, y = 0.99, gp = gpar(fontsize = 20, fontface = "bold"), just = "center")
  
  # Return the combined plot
  return(grid_plot)
}

# Save Real Components Plot
real_plot <- plot_data_gg(Re(r_fft), Re(py_fft), color1 = 'blue', color2 = 'red', label1 = 'R', label2 = 'Py', title="Real Components of FFT Data")
ggsave("real_fft.png", plot = real_plot, width = 10, height = 10)

# Save Imaginary Components Plot
imag_plot <- plot_data_gg(Im(r_fft), Im(py_fft), color1 = 'blue', color2 = 'red', label1 = 'R', label2 = 'Py', title="Imaginary Components of FFT Data")
ggsave("img_fft.png", plot = imag_plot, width = 10, height = 10)
