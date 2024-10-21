library(vnpc.avg)
library(reticulate)
use_python("/Users/avaj0001/Documents/projects/deep_gw_pe_followup/venv/bin/python")

set.seed(0)

n_segments = 1

# Load VNPC results
load('/Users/avaj0001/Documents/projects/vnpc_avg/example/python_comparisons/et_results.RData')
dims <- dim(mcmc_vnp_avg$psd.median)
d <- dims[1]
n_freq_per_seg <- dims[3]

# format R data for python
psdq <- array(0, dim = c(3,  d, d, n_freq_per_seg))
psdq[1, , , ] <-  mcmc_vnp_avg$psd.u05[ , , 1:n_freq_per_seg]  # Lower quantile
psdq[2, , , ] <- mcmc_vnp_avg$psd.median[ , , 1: n_freq_per_seg]       # PSD median
psdq[3, , , ] <- mcmc_vnp_avg$psd.u95[ , , 1:n_freq_per_seg]    # Upper quantile
psdq <- aperm(psdq, c(1, 4, 2, 3)) # from (3, d,d, nfreq) -> (3, nfreq, d,d)
# drop freq[0] as it is 0
psdq <- psdq[,2:n_freq_per_seg,,]
py_vnpc_psdq <- r_to_py(psdq)


# Load raw data
load('/Users/avaj0001/Documents/projects/vnpc_avg/example/XYZ_nlized.RData')
data <- XYZ_nlized
data <- data[1:1024,]
py_data <- r_to_py(data)



py_run_string("
from sgvb_psd.psd_estimator import PSDEstimator
from sgvb_psd.postproc import plot_psdq, plot_peridogram, format_axes
import numpy as np

sgvb_runner = PSDEstimator(
  x=r.py_data,
  N_theta=40,
  nchunks=1,
  ntrain_map=1000,
  fs=np.pi*2,
)
sgvb_runner.run(lr=0.003)
")



py_run_string("
from sgvb_psd.postproc import plot_psdq, plot_peridogram, format_axes

# rescale freq from [0, pi] to [0, 512]
freq = sgvb_runner.freq * 1024.0 / np.pi

# rescale PSDs (1e+22 is approx std for data)
py_psd = sgvb_runner.pointwise_ci / 1e+22**2 / (1024.0/ np.pi)
r_psd = r.py_vnpc_psdq / 1e+22**2 / (1024.0/ np.pi)
r_psd = r_psd.astype(np.complex128)
# make lower triangle elements imaginary by multiplying by 1j
for i in range(3):
  for j in range(3):
    if i > j:
      r_psd[:, :, i, j] *= 1j

# rescale pdgrm
pdgrm = sgvb_runner.pdgrm / 1e+22**2 / (1024.0/ np.pi)
pdgrm_freq = sgvb_runner.pdgrm_freq * 1024.0 / np.pi

kwgs = dict(
  off_symlog=True,
  xlims=[5, 128],
  diag_ylims=[1e-52, 1e-46],
  off_ylims=[-3e-47, 3e-47],
  labels='XYZ',
)
axs = plot_peridogram(pdgrm, pdgrm_freq)
axs = plot_psdq(py_psd,freq,axs=axs)
axs = plot_psdq(r_psd, freq, color='tab:red', axs=axs)
format_axes(axs, **kwgs)
axs[0,0].get_figure().savefig('comparison.png')
")









## PLOT


