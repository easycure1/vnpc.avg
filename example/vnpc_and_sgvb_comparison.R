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




py_psd_estimator <- import("sgvb_psd.psd_estimator")$PSDEstimator
plot_psd <- import("sgvb_psd.postproc.plot_psd")$plot_psdq

py_psd_optim <- py_psd_estimator(
  x = r_to_py(data),
  N_theta = 50,
  nchunks = as.integer(n_segments),
  ntrain_map = 1000,
  max_hyperparm_eval = 1,
  fs = 2 * pi,
)
py_psd_optim$run(lr = 0.003)


## PLOT

# format R data for python
psdq <- array(0, dim = c(3,  d, d, n_freq_per_seg))
psdq[1, , , ] <-  mcmc_vnp_avg$psd.u05[ , , 1:n_freq_per_seg]  # Lower quantile
psdq[2, , , ] <- mcmc_vnp_avg$psd.median[ , , 1: n_freq_per_seg]       # PSD median
psdq[3, , , ] <- mcmc_vnp_avg$psd.u95[ , , 1:n_freq_per_seg]    # Upper quantile
psdq <- aperm(psdq, c(1, 4, 2, 3)) # from (3, d,d, nfreq) -> (3, nfreq, d,d)
py$psdq <- psdq

py_run_string("
axs = sgvb.plot(
  off_symlog=False,
  xlims=(0, 3.14),
  off_ylims=(-0.2, 0.2),
)
axs[0,0].get_figure().savefig('sgvb.png')

axs = plot_psdq(
  psdq,
  sgvb.freq,
  color='tab:red',
  axs=axs,
  off_symlog=False,
  off_ylims=(-0.2, 0.2),
  xlims=(0, 3.14),
)
axs[0,0].get_figure().savefig('comparison.png')
")


