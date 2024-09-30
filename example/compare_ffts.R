library(vnpc.avg)
library(reticulate)
use_python("/Users/avaj0001/Documents/projects/deep_gw_pe_followup/venv/bin/python")
set.seed(0) 


n <- 2 ^ 9
d <- 2
n_segments <- as.integer(4)
n_freq_per_seg = n / n_segments / 2
print(paste("Simulating data of", d, "x", n, " (dividing into ", n_segments, ")"))
param_ar <- rbind(c(.1, 0, 0, 0), c(0, -.3, 0, -.5))
param_ma <- matrix(nrow=2, ncol=0)
param_sigma <- matrix(data=c(1, .9, .9, 1), nrow=2, ncol=2)
data <- beyondWhittle::sim_varma(model=list(ar=param_ar, ma=param_ma, sigma=param_sigma), n=n, d=d)


omegaFreq <- function(n) {
  return(2 / n * (1:(n / 2 + 1) - 1))
}

seg_n <- n_segments
d <- ncol(data)
n <- nrow(data) / seg_n
n_init <- n ## Store the original full length of each segment
omega <- omegaFreq(n)
N <- length(omega)
data_new <- array(rep(NA, n*d*seg_n), dim = c(n, d, seg_n))
for (ii in 1:seg_n) {
  data_new[,,ii] <- data[(n*(ii-1)+1):(n*ii),]
}
noise <- array(rep(NA, n_init*d*seg_n), dim = c(n_init, d, seg_n))
FZ <- array(rep(NA, N*d*seg_n), dim = c(N, d, seg_n))
for (ii in 1:seg_n) {
  noise[,,ii] <- get_noise(data_new[,,ii], theta[,i])
  FZ[,,ii] <- mdft(noise[,,ii])[1:N,]
}



compute_chunked_fft <- import("sgvb_psd.backend.analysis_data")$compute_chunked_fft
py_fft = compute_chunked_fft(
  x=r_to_py(data), 
  nchunks=n_segments, 
  fs=0.5
)[1]
