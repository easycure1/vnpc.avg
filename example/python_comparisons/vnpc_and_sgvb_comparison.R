# pip install sgvb_psd (version >= 0.0.5)
# install.packages("reticulate")
library(vnpc.avg)
library(reticulate)
use_python("/Users/avaj0001/Documents/projects/deep_gw_pe_followup/venv/bin/python")
set.seed(0) 


n <- 2 ^ 9
d <- 2
n_segments <- 4
n_freq_per_seg = n / n_segments / 2
print(paste("Simulating data of", d, "x", n, " (dividing into ", n_segments, ")"))
param_ar <- rbind(c(.1, 0, 0, 0), c(0, -.3, 0, -.5))
param_ma <- matrix(nrow=2, ncol=0)
param_sigma <- matrix(data=c(1, .9, .9, 1), nrow=2, ncol=2)
data <- beyondWhittle::sim_varma(model=list(ar=param_ar, ma=param_ma, sigma=param_sigma), n=n, d=d)



### RUN MCMC
Ntotal <- 300 # TODO 80000
burnin <- 100  # TODO 30000
Nadaptive <- burnin
thin <- 1 # TODO 5
print_interval <- 50 # TODO 10000
mcmc_vnp_avg <- gibbs_vnpc_avg(
  data = data,
  seg_n = n_segments,
  Ntotal = Ntotal,
  burnin = burnin,
  print_interval = print_interval
)
save(mcmc_vnp_avg, file = "mcmc_results_test.RData")




## PLOT

# format R data for python
psdq <- array(0, dim = c(3,  d, d, n_freq_per_seg))
psdq[1, , , ] <-  mcmc_vnp_avg$psd.u05[ , , 1:n_freq_per_seg]  # Lower quantile
psdq[2, , , ] <- mcmc_vnp_avg$psd.median[ , , 1: n_freq_per_seg]       # PSD median
psdq[3, , , ] <- mcmc_vnp_avg$psd.u95[ , , 1:n_freq_per_seg]    # Upper quantile
psdq <- aperm(psdq, c(1, 4, 2, 3)) # from (3, d,d, nfreq) -> (3, nfreq, d,d)
vnpc_psdq <- r_to_py(psdq)
py_data <- r_to_py(data)
nchunks <- r_to_py(n_segments)



py_run_string("
from sgvb_psd.psd_estimator import PSDEstimator
from sgvb_psd.postproc import plot_psdq
import numpy as np

sgvb_runner = PSDEstimator(
  x=r.py_data,
  N_theta=40,
  nchunks=int(r.nchunks),
  ntrain_map=1000,
  fs=2 * np.pi,
)
sgvb_runner.run(lr=0.003)


axs = sgvb_runner.plot(
  off_symlog=False,
  xlims=(0, 3.14),
  off_ylims=(-0.2, 0.2),
)
axs = plot_psdq(
  r.vnpc_psdq,
  sgvb_runner.freq,
  color='tab:red',
  axs=axs,
  off_symlog=False,
  off_ylims=(-0.2, 0.2),
  xlims=(0, 3.14),
)
axs[0,0].get_figure().savefig('comparison.png')
")


