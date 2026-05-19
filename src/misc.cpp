#include "misc.h"
#include "matrix_cube.h"
using namespace Rcpp;
using namespace std;

//' Build an n times n Toeplitz matrix from the 
//' autocovariance values gamma(0),...,gamma(n-1)
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix acvMatrix(NumericVector acv) {
  unsigned p = acv.size();
  NumericMatrix m(p, p);
  for (int i=0; i < p; ++i) {
    for (int j=0; j < p; ++j) {
      unsigned index = abs(i-j);
      m(i,j) = acv[index];
    }
  }
  return(m);
}

//' Build an nd times nd Block Toeplitz matrix from the
//' (d times d) autocovariances gamma(0),...,gamma(n-1)
//' @keywords internal
// [[Rcpp::export]]
arma::mat acvBlockMatrix(arma::mat acv) {
  const unsigned d = acv.n_rows;
  const unsigned p = acv.n_cols / d;
  arma::mat m(p*d, p*d);
  for (int i=0; i < p; ++i) {
    for (int j=0; j < p; ++j) {
      unsigned index = abs(i-j);
      m.submat(i*d,j*d,(i+1)*d-1,(j+1)*d-1) = acv.submat(0,
               index*d, d-1, (index+1)*d-1);
    }
  }
  return(m);
}

//' Computing acceptance rate based on trace
//' Note: Only use for traces from continous distributions!
//' @keywords internal
// [[Rcpp::export]]
double acceptanceRate(NumericVector trace) {
  unsigned rejections = 0;
  for (unsigned i=1; i < trace.length(); ++i) {
    rejections += (trace[i]==trace[i-1]);
  }
  double rejectionRate = (double)rejections / (double)trace.length();
  return 1 - rejectionRate;
}

//' Hermitian conjugate of a matrix array
//' @keywords internal
// [[Rcpp::export]]
arma::cx_cube trans_cube(ComplexVector f_) { // ok
  const arma::cx_cube f = cx_cube_from_ComplexVector(f_);
  arma::cx_cube res(f.n_rows, f.n_cols, f.n_slices);
  for (unsigned j=0; j<f.n_slices; ++j) {
    res.slice(j) = arma::trans(f.slice(j)); // Hermitian conjugate
  }
  return res;
}

//' Matrix products of a matrix array
//' @keywords internal
// [[Rcpp::export]]
arma::cx_cube mult_cube(ComplexVector a_, ComplexVector b_) { // ok
  const arma::cx_cube a = cx_cube_from_ComplexVector(a_);
  const arma::cx_cube b = cx_cube_from_ComplexVector(b_);
  arma::cx_cube c(a.n_rows, a.n_cols, a.n_slices);
  for (unsigned j=0; j<a.n_slices; ++j) {
    c.slice(j) = a.slice(j) * b.slice(j);
  }
  return c;
}

