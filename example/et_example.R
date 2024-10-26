load('XYZ_nlized.RData')
data <- XYZ_nlized

Ntotal <- 3000 # TODO 80000
burnin <- 1000  # TODO 30000
Nadaptive <- burnin
thin <- 1 # TODO 5
print_interval <- 500 # TODO 10000



mcmc_vnp_avg <- vnpc.avg::gibbs_vnpc_avg(data = data[1:32768,],
                                         samp_freq = 16, #sampling frequency
                                         truncation = T,
                                         trunc_freq_lim = c(5, 128), #truncation bounds
                                         seg_n = 1,
                                         Ntotal = Ntotal,
                                         burnin = burnin)
save(mcmc_vnp_avg, file = "et_results.RData")