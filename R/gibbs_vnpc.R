#' Gibbs sampler for multivariate Bayesian inference with the nonparametrically corrected VAR likelihood
#'
#' Obtain samples of the posterior of the multivariate corrected likelihood in conjuction with an Hpd AGamma process prior on the spectral density matrix
#' @param data numerical matrix
#' @param samp_freq sampling frequency (default = 2*pi)
#' @param trunc_freq sampling frequency for truncated data (default = 2*pi)
#' @param seg_n number of segments (integer >= 1)
#' @param truncation flag indicating whether the data needs to be truncated (default = FALSE).
#' @param trunc_freq_lim frequency bounds of the truncated data, an positive integer or a 2-dimensional vector with lower and upper bounds being integers, used for truncation = TRUE only.
#' @param corrected flag indicating whether the corrected likelihood is used (default = FALSE).
#' @param var.order VAR order for the parametric working model
#' @param Ntotal total number of iterations to run the Markov chain
#' @param burnin number of initial iterations to be discarded
#' @param thin thinning number (postprocessing)
#' @param print_interval number of iterations, after which a status is printed to console
#' @param numerical_thresh lower (numerical pointwise) bound for the eigenvalues of the spectral density
#' @param adaption.N total number of iterations, in which the proposal variances (of r, U and VAR coefficients) are adapted
#' @param adaption.batchSize  batch size of proposal adaption
#' @param adaption.tar target acceptance rate for adapted parameters
#' @param eta AGamma process parameter, real number > ncol(data) - 1
#' @param omega_fun AGamma process parameter, positive constant
#' @param Sigma_fun AGamma process parameter, Hpd matrix
#' @param k.theta prior parameter for polynomial degree k (propto exp(-k.theta*k*log(k)))
#' @param kmax upper bound for polynomial degree of Bernstein-Dirichlet mixture (can be set to Inf, algorithm is faster with kmax<Inf due to pre-computation of basis functions, but values 500<kmax<Inf are very memory intensive)
#' @param trunc_l,trunc_r left and right truncation of Bernstein polynomial basis functions, 0<=trunc_l<trunc_r<=1
#' @param coars flag indicating whether coarsened or default bernstein polynomials are used (see Appendix E.1 in Ghosal and van der Vaart 2017)
#' @param L truncation parameter of Gamma process
#' @param mu_beta prior parameter for VAR coefficients, stacked numerical vector
#' @param V_beta prior parameter for VAR coefficients, Hpd matrix
#' @param sqrt_d flag indicating whether the regular square root of a Hermitian matrix is used, otherwise the Cholesky decomposition is used to approximate the square root (see Remark 5.2 in Liu (2023)) 
#'
#' @return list containing the following fields:
#'  \item{data}{(centerd) data}
#'  \item{psd.median, psd.mean}{(pointwise) posterior median and mean}
#'  \item{psd.p05, psd.p95}{90\% pointwise credible interval}
#'  \item{psd.u05, psd.u95}{90\% uniform credible interval}
#'  \item{coherence.median}{(pointwise) posterior median coherence}
#'  \item{coherence.p05, coherence.p95}{(pointwise) 90\% pointwise credible interval for the coherence}
#' @references Y. Liu (2023)
#' \emph{A Nonparametrically corrected likelihood for Bayesian spectral analysis of multivariate time series}
#' PhD thesis, University of Auckland
#' <https://hdl.handle.net/2292/65154> 
#' @importFrom Rcpp evalCpp
#' @useDynLib vnpc.avg, .registration = TRUE
#' @export
#'
gibbs_vnpc_avg <- function(data=NULL,
                           mpg_avg,
                           Nb,
                           samp_freq=2*pi,
                           trunc_freq=2*pi,
                           seg_n=1,
                           truncation=FALSE,
                           trunc_freq_lim=NULL,
                           bspline=FALSE,
                           degree_bs=3,
                           corrected=FALSE,
                           toggle=FALSE,
                           f_param_avg=NULL,
                           var.order=NULL,
                           Ntotal,
                           burnin,
                           thin=1,
                           print_interval=100,
                           numerical_thresh=1e-12,
                           adaption.N=burnin,
                           adaption.batchSize=50,
                           adaption.tar=0.44,
                           # eta=ncol(data),
                           eta=dim(mpg_avg)[2],
                           omega_fun=create_omega_fun_from_beta_density(1,1,1),
                           Sigma_fun=my_Sigma_fun,
                           k.theta=0.01,
                           kmax=100*coars + 500*(!coars),
                           trunc_l=0.1,
                           trunc_r=0.9,
                           coars=F,
                           H0.alpha=NULL,
                           H0.beta=NULL,
                           MH=NULL,
                           L=20,
                           mu_beta=NULL,
                           V_beta=NULL,
                           sqrt_d=F) {
  # if (!is.matrix(data) || !is.numeric(data)) {
  #   stop("'data' must be numeric matrix with d columns and n rows")
  # }
  # 
  # d <- ncol(data)
  # if (d<2) {
  #   stop("This function is not suited for univariate time series. Use gibbs_npc instead")
  # }
  # 
  # if (max(abs(apply(data,2,mean,na.rm=T))) > 1e-4) {
  #   data <- apply(data,2,center,na.rm=T)
  #   warning("Data has been mean centered")
  # }
  d <- dim(mpg_avg)[2]
  if (eta <= d-1) {
    stop("eta must be a number greater than d-1")
  }
  # if (omega_fun <= 0) {
  #  stop("omega must be a positive number")
  # }
  mcmc_params <- list(Ntotal=Ntotal,
                      burnin=burnin,
                      thin=thin,
                      print_interval=print_interval,
                      numerical_thresh=1e-12,
                      Nadaptive=adaption.N,
                      adaption.batchSize=adaption.batchSize,
                      adaption.targetAcceptanceRate=adaption.tar,
                      verbose=F)
  
  # if (truncation) {
  #   if (length(trunc_freq_lim) == 1) {
  #     trunc_N <- trunc_freq * trunc_freq_lim
  #   } else if (length(trunc_freq_lim) == 2) {
  #     trunc_N_upper <- trunc_freq_lim[2] * trunc_freq
  #     trunc_N_lower <- trunc_freq_lim[1] * trunc_freq
  #     trunc_N <- trunc_N_upper - trunc_N_lower
  #   }
  #   L <- max(L, length(trunc_N) ^ (1 / 3))
  # } else {
  #   L <- max(L, (length(data) / seg_n) ^ (1 / 3))
  # }
  # L <- ceiling(L)
  
  # if (corrected) {
  #   mu_beta <- rep(1e-4, ncol(data)*ncol(data)*var.order)
  #   V_beta <- diag(ncol(data)*ncol(data)*var.order)*1e4
  # }
  prior_params <- list(prior.cholesky=F,
                       var.order=var.order,
                       eta=eta,
                       omega_fun=omega_fun,
                       Sigma_fun=Sigma_fun,
                       k.theta=k.theta,
                       kmax=kmax,
                       bernstein_l=trunc_l, # note
                       bernstein_r=trunc_r, # note
                       coarsened=coars, # note
                       L=L,
                       toggle=toggle,
                       prior.q=T,
                       mu_beta=mu_beta,
                       V_beta=V_beta,
                       sqrt_d=sqrt_d)
  model_params <- psd_dummy_model()
  mcmc_VNPC <- gibbs_m_avg_nuisance(data=data,
                                    mpg_avg=mpg_avg,
                                    Nb=Nb,
                                    samp_freq=samp_freq,
                                    trunc_freq=trunc_freq,
                                    mcmc_params=mcmc_params,
                                    seg_n=seg_n,
                                    truncation=truncation,
                                    trunc_freq_lim=trunc_freq_lim,
                                    bspline=bspline,
                                    degree_bs=degree_bs,
                                    H0.alpha=H0.alpha,
                                    H0.beta=H0.beta,
                                    MH=MH,
                                    corrected=corrected,
                                    f_param_avg=f_param_avg,
                                    prior_params=prior_params,
                                    model_params=model_params)
  #mcmc_VNPC <- reduceMemoryStorage_matrixGamma(mcmc)
  #return(mcmc_VNPC)
  return(structure(list(data=data,
                        k=mcmc_VNPC$k,
                        psd.post=mcmc_VNPC$fpsd.sample,
                        psd.median=mcmc_VNPC$fpsd.s,
                        psd.p05=mcmc_VNPC$fpsd.s05,
                        psd.p95=mcmc_VNPC$fpsd.s95,
                        psd.mean=mcmc_VNPC$fpsd.mean,
                        psd.u05=mcmc_VNPC$fpsd.uuci05,
                        psd.u95=mcmc_VNPC$fpsd.uuci95,
                        coherence.median = mcmc_VNPC$coherence.s,
                        coherence.p05 = mcmc_VNPC$coherence.s05,
                        coherence.p95 = mcmc_VNPC$coherence.s95,
                        post=mcmc_VNPC$lpostTrace)))
}





























