#ifndef __VARMA_INCLUDED__
#define __VARMA_INCLUDED__

#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::cx_cube varma_transfer2psd(Rcpp::ComplexVector transfer_ar_,
                                 Rcpp::ComplexVector transfer_ma_,
                                 arma::mat sigma);
arma::cx_cube varma_transfer2psd_cpp(arma::cx_cube transfer_ar,
                                     arma::cx_cube transfer_ma,
                                     arma::mat sigma);
arma::cx_cube transfer_polynomial(Rcpp::NumericVector lambda, arma::mat coef);
arma::mat epsilon_var(arma::mat zt, arma::mat ar);
double sldmvnorm(arma::mat z_t, arma::mat Sigma);
arma::cx_cube psd_varma_cpp(Rcpp::NumericVector freq,
                            arma::mat ar,
                            arma::mat ma,
                            arma::mat sigma);

#endif
