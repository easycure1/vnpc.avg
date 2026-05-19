#'@title Gibbls sampler for VNP based on average periodogram
#'@description
#'Gibbs sampler for multivariate Bayesian nonparametric inference with corrected Whittle likelihood
#'
#'@importFrom Rcpp evalCpp
#'@useDynLib vnpc.avg, .registration = TRUE
#'@export

gibbs_vnpc_avg <- function(mpg_avg,
                           f_param_avg,
                           Nb,
                           Ntotal,
                           burnin,
                           thin=1,
                           print_interval=100,
                           numerical_thresh=1e-7,
                           adaption.N=burnin,
                           adaption.batchSize=50,
                           adaption.tar=.44,
                           eta=dim(mpg_avg)[2],
                           omega=dim(mpg_avg)[2],
                           Sigma=1e4*diag(dim(mpg_avg)[2]),
                           k.theta=0.01,
                           kmax = 100*coars + 500*(!coars),
                           trunc_l = 0.1,
                           trunc_r = 0.9,
                           coars=F,
                           L = max(20, (dim(mpg_avg)[1] * 2) ^ (1 / 3)),
                           sqrt_d = F) {  
  
  # if (!is.matrix(data) || !is.numeric(data)) {
  #   stop("'data' must be numeric matrix with d columns and n rows")
  # }
  
  d <- dim(mpg_avg)[2]
  # if (d<2) {
  #   stop("This function is not suited for univariate time series. Use gibbs_NP instead")
  # }
  
  # if (max(abs(apply(data,2,mean,na.rm=T))) > 1e-4) {
  #   data <- apply(data,2,center,na.rm=T)
  #   warning("Data has been mean centered")
  # }
  
  if (eta <= d-1) {
    stop("eta must be a number greater than d-1")
  }
  if (omega <= 0) {
    stop("omega must be a positive number")
  }
  if (class(Sigma)[1]!="matrix" || (!is_hpd(Sigma)) || any(dim(Sigma)!=c(d,d))) {
    stop("Sigma must be a Hermitian positive definite d times d matrix")
  }
  
  cl <- match.call
  
  mcmc_params <- list(Ntotal=Ntotal,
                      burnin=burnin,
                      thin=thin,
                      print_interval=print_interval,
                      numerical_thresh=numerical_thresh,
                      verbose=F,
                      Nadaptive=adaption.N,
                      adaption.batchSize=adaption.batchSize,
                      adaption.targetAcceptanceRate=adaption.tar)
  prior_params <- list(eta=eta,
                       omega=omega,
                       Sigma=Sigma,
                       k.theta=k.theta,
                       kmax=kmax,
                       bernstein_l=trunc_l, # note
                       bernstein_r=trunc_r, # note
                       coarsened=coars, # note
                       L=L,
                       f_param_avg=f_param_avg,
                       sqrt_d=sqrt_d)
  model_params <- psd_dummy_model()
  
  # Call internal MCMC algorithm
  mcmc_VNP <- gibbs_multivariate_nuisance(mpg_avg=mpg_avg,
                                          Nb=Nb,
                                          mcmc_params=mcmc_params, 
                                          corrected=T, 
                                          prior_params=prior_params, 
                                          model_params=model_params)
  
  #return(mcmc_VNPC)
  return(structure(list(call=cl,
                        mpg_avg=mpg_avg,
                        psd.median=complexValuedPsd(mcmc_VNP$fpsd.s),
                        psd.p05=complexValuedPsd(mcmc_VNP$fpsd.s05),
                        psd.p95=complexValuedPsd(mcmc_VNP$fpsd.s95),
                        psd.mean=complexValuedPsd(mcmc_VNP$fpsd.mean),
                        psd.u05=complexValuedPsd(mcmc_VNP$fpsd.uuci05),
                        psd.u95=complexValuedPsd(mcmc_VNP$fpsd.uuci95),
                        coherence.median = mcmc_VNP$coherence.s,
                        coherence.p05 = mcmc_VNP$coherence.s05,
                        coherence.p95 = mcmc_VNP$coherence.s95,
                        # missing_values=mcmc_VNP$missingValues_trace,
                        k = mcmc_VNP$k,
                        r = mcmc_VNP$r,
                        Z = mcmc_VNP$Z,
                        U = mcmc_VNP$U,
                        lpost=mcmc_VNP$lpostTrace,
                        algo="gibbs_vnpc_avg"),
                        class="gibbs_psd"))
}



