#ifndef __GIBBS_VNPC_HELP_INCLUDED__
#define __GIBBS_VNPC_HELP_INCLUDED__


#if !defined(ARMA_WARN_LEVEL)
  #define ARMA_WARN_LEVEL 0
#endif


#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "bernstein_gamma_psd.h"
#include "matrix_cube.h"

double lprior_parametric(Rcpp::List phi);
double llikelihood_corrected_q(const arma::cx_mat& FZ,
                               arma::mat ar,
                               const arma::cx_cube& f_param_half,
                               arma::mat sigma,
                               const bernsteinGammaPsd& f,
                               bool sqrt_d);
double llike_var_partial_cpp(const arma::mat& FCFZ, 
                             arma::mat ar, 
                             arma::mat sigma);
arma::cx_mat get_CFZ_q(const arma::cx_mat FZ,
                       const arma::cx_cube q,
                       const arma::cx_cube f_param_half);
arma::cx_mat get_CFZ_q_sq(const arma::cx_mat FZ,
                          const arma::cx_cube q,
                          const arma::cx_cube f_param_half);					   
arma::cx_cube chol_cube(arma::cx_cube f_, bool excludeBoundary);
arma::cx_cube sqrt_cube(arma::cx_cube f_, bool excludeBoundary);
double llike_corrected_avg(const arma::cx_cube& mpg_avg, 
                           const arma::cx_cube& f, 
                           const arma::cx_cube& f_param_avg_half, 
                           const int& Nb);
double logdet_cube_sum(const arma::cx_cube f);
arma::mat phiFromBeta_normalInverseWishart_cpp(arma::vec& beta,
                                               const int& K,
                                               const int& p);

#endif