#ifndef __GIBBS_UTIL_INCLUDED__
#define __GIBBS_UTIL_INCLUDED__

#include <RcppArmadillo.h>
#include <Rcpp.h>

Rcpp::NumericVector pFromV(Rcpp::NumericVector v);
Rcpp::NumericVector vFromP(Rcpp::NumericVector p, const double eps);

#endif