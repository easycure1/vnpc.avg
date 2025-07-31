#include "gibbs_util.h"

using namespace Rcpp;
using namespace std;


//' Get p from V (DP inverse stick breaking)
//' @keywords internal
// [[Rcpp::export]]
NumericVector pFromV(NumericVector v) {
  unsigned L = v.size();
  NumericVector p(L + 1);
  double currentProduct = 1.0;
  double pSum = 0.0;
  for (unsigned l = 0; l < L; ++l) {
    p[l + 1] = currentProduct * v[l];
    currentProduct *= (1.0 - v[l]);
    pSum += p[l + 1];
  }
  p[0] = std::max(1.0 - pSum, 0.0); // account for numerical instabilities
  return p;
}


//' Get v From P (DP inverse stick breaking)
//' Note: p is assumed to have length L, i.e. it does NOT contain p_0
//' @keywords internal
// [[Rcpp::export]]
NumericVector vFromP(NumericVector p, const double eps=1e-8) {
  unsigned L = p.size();
  NumericVector v(L);
  double currentProduct = 1.0;
  for (unsigned l = 0; l < L; ++l) {
    v[l] = std::min(std::max(p[l] / currentProduct, eps),1.0-eps); // numerical stability
    //v[l] = p[l] / currentProduct;
    currentProduct *= (1.0 - v[l]);
  }
  return v;
}
