from sgvb_psd.psd_estimator import PSDEstimator
from sgvb_psd.postproc import plot_psdq, plot_peridogram, format_axes

from rpy2 import robjects
import numpy as np


def load_vnpc_results():
  robjects.r['load']('et_results.RData')
  mcmc_vnp_avg = robjects.r['mcmc_vnp_avg']

  # dims < - dim(mcmc_vnp_avg$psd.median)
  # d < - dims[1]
  # n_freq_per_seg < - dims[3]
  #
  # # format R data for python
  # psdq < - array(0, dim=c(3, d, d, n_freq_per_seg))
  # psdq[1,,, ] < -  mcmc_vnp_avg$psd.u05[, , 1:n_freq_per_seg]  # Lower quantile
  # psdq[2,,, ] < - mcmc_vnp_avg$psd.median[, , 1: n_freq_per_seg]  # PSD median
  # psdq[3,,, ] < - mcmc_vnp_avg$psd.u95[, , 1:n_freq_per_seg]  # Upper quantile
  # psdq < - aperm(psdq, c(1, 4, 2, 3))  # from (3, d,d, nfreq) -> (3, nfreq, d,d)
  # # drop freq[0] as it is 0
  # psdq < - psdq[, 2:n_freq_per_seg, , ]


# Load R data
robjects.r['load']('XYZ_nlized.RData')


data = np.array(robjects.r['XYZ_nlized'])

print("Loaded data shape: ", data.shape)
data = data[:1024]
print("Truncating data to: ", data.shape)



sgvb_runner = PSDEstimator(
  x=data,
  N_theta=40,
  nchunks=1,
  ntrain_map=100,
  fs=np.pi*2,
)
sgvb_runner.run(lr=0.003)



freq = sgvb_runner.freq
# rescale freq from [0, pi] to [0, 512]
freq_scaling = 512.0 / np.pi
freq = freq * freq_scaling
# rescale psd
psd = sgvb_runner.pointwise_ci

pdgrm = sgvb_runner.pdgrm / 1e+22**2 / (1024.0/ np.pi)
pdgrm_freq = sgvb_runner.pdgrm_freq * freq_scaling


py_psd = sgvb_runner.pointwise_ci / 1e+22**2 / (1024.0/ np.pi)


kwgs = dict(
  off_symlog=True,
  xlims=[5, 128],
  diag_ylims=[1e-52, 1e-46],
  off_ylims=[-3e-47, 3e-47],
)
axs = plot_peridogram(pdgrm, pdgrm_freq)
axs = plot_psdq(
  py_psd,
  freq,
  axs=axs,
)
format_axes(axs, **kwgs)
axs[0,0].get_figure().savefig('sgvb.png')










