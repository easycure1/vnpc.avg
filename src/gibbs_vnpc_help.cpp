#include "gibbs_vnpc_help.h"
#include "gibbs_vnp_avg_help.h"
#include "varma.h"


using namespace std;
using namespace Rcpp;


// Unnormalized prior for the parametric (VAR) part of the semiparametric procedure
// (Normal-Inverse-Wishart for beta)
double lprior_parametric(List phi) {
  arma::vec beta = Rcpp::as<arma::vec>(phi["beta"]);
  arma::vec mu_beta = Rcpp::as<arma::vec>(phi["mu_beta"]);
  arma::mat V_beta_inv = Rcpp::as<arma::mat>(phi["V_beta_inv"]);
  
  // Use as_scalar() to convert the matrix output to a scalar
  double res = -.5 * arma::as_scalar((beta - mu_beta).t() * V_beta_inv * (beta - mu_beta));
  return res;
}

// Log corrected likelihood using the q-parametrisation
// Note: omega=0 and omega=pi excluded!
double llikelihood_corrected_q(const arma::cx_mat& FZ,
                               arma::mat ar,
                               const arma::cx_cube& f_param_half,
                               arma::mat sigma,
                               const bernsteinGammaPsd& f,
                               bool sqrt_d) {
  const int N = FZ.n_rows;
  const int d = FZ.n_cols;
  arma::cx_mat CFZ;
  if (sqrt_d) {
    CFZ = get_CFZ_q_sq(FZ, f.eval(), f_param_half);
  } else {
    CFZ = get_CFZ_q(FZ, f.eval(), f_param_half);
  }
  const arma::mat FCFZ(arma::real(midft_cpp(CFZ)));
  double res = -logdet_cube_sum(f.eval())  + llike_var_partial_cpp(FCFZ, ar, sigma);
  return res;
}

// Log partial likelihood of VAR(p), CPP-internal only
// Note: Functionality equivalent to llike_var_partial()
double llike_var_partial_cpp(const arma::mat& FCFZ, 
                             arma::mat ar, 
                             arma::mat sigma) {
  arma::mat epsilon_t(epsilon_var(FCFZ, ar));
  double res = sldmvnorm(epsilon_t, sigma);
  return res;
}

//' Build the correction matrix in the frequency domain based on the q-parametrization with the Cholesky decomposition
//' @keywords internal
// [[Rcpp::export]]
 arma::cx_mat get_CFZ_q(const arma::cx_mat FZ, const arma::cx_cube q,
                        const arma::cx_cube f_param_half) {
   const int N = FZ.n_rows;
   // const arma::cx_cube q = cx_cube_from_ComplexVector(q_);
   // const arma::cx_cube f_param_half = cx_cube_from_ComplexVector(f_param_half_);
   arma::cx_mat res(FZ.n_rows, FZ.n_cols);
   res.row(0) = arma::cx_rowvec(FZ.n_cols, arma::fill::zeros); // exclude boundaries
   res.row(FZ.n_rows-1) = arma::cx_rowvec(FZ.n_cols, arma::fill::zeros); // exclude boundaries
   for (unsigned j=1; j<N-1; ++j) {
     res.row(j) = (
       f_param_half.slice(j) *
         arma::inv(f_param_half.slice(j) * trans(arma::chol(q.slice(j)))) *
         FZ.row(j).st()).st();
   }
   return res;
}

//' Build the correction matrix in the frequency domain based on the q-parametrization with square root of matrices
//' @keywords internal
// [[Rcpp::export]]
arma::cx_mat get_CFZ_q_sq(const arma::cx_mat FZ, const arma::cx_cube q,
                           const arma::cx_cube f_param_half) {
   const int N = FZ.n_rows;
   // const arma::cx_cube q = cx_cube_from_ComplexVector(q_);
   // const arma::cx_cube f_param_half = cx_cube_from_ComplexVector(f_param_half_);
   arma::cx_mat res(FZ.n_rows, FZ.n_cols);
   res.row(0) = arma::cx_rowvec(FZ.n_cols, arma::fill::zeros); // exclude boundaries
   res.row(FZ.n_rows-1) = arma::cx_rowvec(FZ.n_cols, arma::fill::zeros); // exclude boundaries
   for (unsigned j=1; j<N-1; ++j) {
     res.row(j) = (
       f_param_half.slice(j) * arma::inv(arma::sqrtmat(f_param_half.slice(j) *
         q.slice(j) * f_param_half.slice(j))) * FZ.row(j).st()).st();
   }
   return res;
}

//' Cholesky docomposition representation of a matrix array. See Remark 5.2
//' @keywords internal
// [[Rcpp::export]]
arma::cx_cube chol_cube(arma::cx_cube f_, bool excludeBoundary) { // ok
   // const arma::cx_cube f = cx_cube_from_ComplexVector(f_);
   const arma::cx_cube f = f_;
   arma::cx_cube f_half(f.n_rows, f.n_cols, f.n_slices); // Carful: No fill
   if (excludeBoundary) {
     f_half.slice(0) = arma::cx_mat(f.n_rows, f.n_cols, arma::fill::zeros);
     f_half.slice(f.n_slices-1) = arma::cx_mat(f.n_rows, f.n_cols, arma::fill::zeros);
   }
   for (unsigned j=excludeBoundary; j < f.n_slices-excludeBoundary; ++j) {
     f_half.slice(j) = trans(arma::chol(f.slice(j)));
   }
   return f_half;
}

//' Square root of a matrix array. See Remark 5.2
//' @keywords internal
// [[Rcpp::export]]
arma::cx_cube sqrt_cube(arma::cx_cube f_, bool excludeBoundary) { // ok
   // const arma::cx_cube f = cx_cube_from_ComplexVector(f_);
   const arma::cx_cube f = f_;
   arma::cx_cube f_half(f.n_rows, f.n_cols, f.n_slices); // Carful: No fill
   if (excludeBoundary) {
     f_half.slice(0) = arma::cx_mat(f.n_rows, f.n_cols, arma::fill::zeros);
     f_half.slice(f.n_slices-1) = arma::cx_mat(f.n_rows, f.n_cols, arma::fill::zeros);
   }
   for (unsigned j=excludeBoundary; j < f.n_slices-excludeBoundary; ++j) {
     f_half.slice(j) = trans(arma::sqrtmat(f.slice(j)));
   }
   return f_half;
}

//' Log likelihood of a pseudo parametric VAR model
//' @keywords internal
// [[Rcpp::export]]
double llike_corrected_avg(const arma::cx_cube& mpg_avg, 
                           const arma::cx_cube& f, 
                           const arma::cx_cube& f_param_avg_half, 
                           const int& Nb) {
  const int N = mpg_avg.n_slices;
  double res(0.0);
  for (int j=1; j<N-1; ++j) {
    arma::cx_mat F = f_param_avg_half.slice(j) * f.slice(j) * arma::trans(f_param_avg_half.slice(j));
    // (optional but recommended) enforce Hermitian symmetry to kill tiny asymmetry
    F = 0.5 * (F + F.t());
    
    // Cholesky (HPD)
    arma::cx_mat L;
    const bool ok = arma::chol(L, F, "lower");
    if (!ok) {
      double eps = 1e-10 * arma::trace(F).real() / F.n_rows;  // scale-aware
      F += eps * arma::eye<arma::cx_mat>(F.n_rows, F.n_cols);
    }
    
    // log det(F) = 2 * sum(log(diag(L)))
    const double logdet = 2.0 * arma::accu(arma::log(arma::real(L.diag())));
    
    std::complex<double> tr = arma::trace(arma::inv(F) * mpg_avg.slice(j));
    res += Nb * (logdet+ tr.real());
  }
  
  return -res;
}

//' Sum of log determinants of a matrix array excluding boundaries
//' @keywords internal
// [[Rcpp::export]]
double logdet_cube_sum(const arma::cx_cube f) { // ok
  const unsigned N = f.n_slices;
  double res = 0;
  for (unsigned j=1; j<N-1; ++j) {
    std::complex<double> log_det_val;
    double log_det_sign;
    arma::log_det(log_det_val,log_det_sign,f.slice(j));
    res += log_det_val.real();
  }
  return res;
}

// Convert vector parametrization (beta) to matrix-parametrization (phi)
// CPP-internal only
// Functionality equivalent to phiFromBeta_normalInverseWishart(...)
arma::mat phiFromBeta_normalInverseWishart_cpp(arma::vec& beta,
                                               const int& K,
                                               const int& p) {
  arma::mat res = arma::reshape(beta, K*p, K).t();
  return res;
}