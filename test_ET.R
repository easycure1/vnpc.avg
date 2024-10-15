load('XYZ_nlized.RData')
data <- XYZ_nlized

Ntotal <- 3000 # TODO 80000
burnin <- 1000  # TODO 30000
Nadaptive <- burnin
thin <- 1 # TODO 5
print_interval <- 500 # TODO 10000

debug(vnpc.avg::gibbs_vnpc_avg)
undebug(vnpc.avg::gibbs_vnpc_avg)

mcmc_vnp_avg <- vnpc.avg::gibbs_vnpc_avg(data = data[1:1024,],
                                         seg_n = 1,
                                         Ntotal = Ntotal,
                                         burnin = burnin,
                                         print_interval = print_interval)
