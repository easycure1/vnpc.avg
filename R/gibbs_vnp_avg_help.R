#' Function for updating the proposal covariance
#' @keywords internal
updateCov <- function(X,covObj=NA) {
  X <- as.numeric(X)
  
  if (all(is.na(covObj))) {
    return(list(
      n = 1L,
      mean = X,
      C = matrix(0, length(X), length(X))
    ))
  }
  
  n_old <- covObj$n
  n_new <- n_old + 1L
  
  mean_old <- covObj$mean
  delta <- X - mean_old
  mean_new <- mean_old + delta / n_new
  
  C_new <- covObj$C + tcrossprod(delta, X - mean_new)
  
  list(
    n = n_new,
    mean = mean_new,
    C = C_new
  )
}


#' Get empirical covariance matrix from a updateCov() object
#' @keywords internal
getCov <- function(covObj, eps = 1e-6) {
  d <- length(covObj$mean)
  
  if (covObj$n < 2L) {
    return(diag(eps, d))
  }
  
  cov_mat <- covObj$C / (covObj$n - 1)
  cov_mat + diag(eps, d)
}