#ifndef __MISC_INCLUDED__
#define __MISC_INCLUDED__

#include <RcppArmadillo.h>
#include <Rcpp.h>

Rcpp::NumericMatrix acvMatrix(Rcpp::NumericVector acv);
arma::mat acvBlockMatrix(arma::mat acv);
double acceptanceRate(Rcpp::NumericVector trace);
arma::cx_cube trans_cube(Rcpp::ComplexVector f_);
arma::cx_cube mult_cube(Rcpp::ComplexVector a_, Rcpp::ComplexVector b_);

#endif
