updateCov <- function(X,covObj=NA) {
  
  if(all(is.na(covObj))) {
    covObj <- list(mean=X,cov=NA,n=1)
    return(covObj)
  }
  
  covObj$n <- covObj$n + 1 # Update number of observations
  
  if(covObj$n==2) {
    X1 <- covObj$mean
    covObj$mean <- X1/2 + X/2
    dX1 <- X1 - covObj$mean
    dX2 <- X - covObj$mean
    covObj$cov <- tcrossprod(dX1,dX1) + tcrossprod(dX2,dX2)
    return(covObj)
  }
  
  dx <- covObj$mean - X # previous mean minus new X
  covObj$cov <- covObj$cov * (covObj$n-2)/(covObj$n-1) + tcrossprod(dx,dx)/covObj$n
  covObj$mean <- covObj$mean*(covObj$n-1)/covObj$n + X/covObj$n
  return(covObj)
}












if ((i < Nadaptive) && (i > 1)) {
  if (i == 1) {
    covObj_r <- updateCov(r[,1], covObj = NA)
    covObj_U <- updateCov(U__phi[,1,,1], covObj = NA) # U (first component suffices)
  } else {
    covObj_r <- updateCov(r[,i], covObj = covObj_r)
    covObj_U <- updateCov(U__phi[,1,,i], covObj = covObj_U)
  }
  ### Update of the proposal covariance matrix for r
  if ((i <= 2*L) || (runif(1, 0, 1) < 0.05)) {
    cov_r <- 0.01/L*diag(L)
  } else {
    cov_r <- 2.38^2/L*covObj_r$cov
  }
  ### Update of the proposal covariance matrix for U
  ### U (first component suffices)
  if ((i <= 2*L) || (runif(1, 0, 1) < 0.05)) {
    cov_U <- 0.01/L*diag(L)
  } else {
    cov_U <- 2.38^2/L*covObj_L$cov
  }
  if (corrected && toggle) {
    batch.param <- param__phi[1,(1:var.order)*(d-1)+1,batch,drop=F]
    batch.param.acceptanceRate <- apply(batch.param, 1, acceptanceRate)
    lsd_param <- lsd_param + ((batch.param.acceptanceRate > adaption.targetAcceptanceRate)*2-1) * adaption.delta
    eps_param <- exp(2*lsd_param)
  }
}







r.old <- r[,i]
r.star <- mvtnorm::rmvnorm(1, mu = log(r.old), sigma = cov_r)
r.star <- exp(r.star)
f.r.star <- lpost_matrixGamma(samp_freq=samp_freq,
                              omega=omega,
                              mpg_avg=mpg_avg,
                              Nb=Nb,
                              FZ=NULL,
                              data=data,
                              r=r.star,
                              U=U[,,,i],
                              Z=Z[,i+1],
                              k=k[,i+1],
                              C_alpha=C_alpha,
                              omega_fun=omega_fun,
                              k.theta=k.theta,
                              db.list=db.list,
                              bspline=bspline,
                              degree_bs=degree_bs,
                              V=V[,i],
                              W=W[,i],
                              H0.alpha=H0.alpha,
                              H0.beta=H0.beta,
                              MH=MH,
                              eta=eta,
                              Sigma_fun=Sigma_fun,
                              corrected=corrected,
                              toggle=toggle,
                              f_param_avg_half=f_param_avg_half,
                              phi=phi.fit,
                              sigma_ar=sigma.fit,
                              prior.q=prior.q,
                              prior.cholesky=prior.cholesky,
                              excludeBoundary=T, # note
                              verbose=verbose)
f.r <- f.store
# Accept / reject
alpha3 <- min(0, f.r.star -
                f.r +
                mvtnorm::dmvnorm(log(r.old), mean = log(r.star), sigma = cov_r, log = T) -
                mvtnorm::dmvnorm(log(r.star), mean = log(r.old), sigma = cov_r, log = T) -
                sum(log(r.old)) + sum(log(r.star)))
if (log(runif(1, 0, 1)) < alpha3) {
  r[,i+1] <- r.star # accept
  r.old <- r.star
  f.store <- f.r.star
} else {
  r[,i+1] <- r[,i] # reject
}







U.old <- U[,,,i]
for (l in 1:L) {
  
  # MH proposal with boundary inverse-reflections (similar to Z's -- see Choudhuri 2004)
  rejectedU <- F
  U__phi.star <- U__phi[,l,i] + runif(length(U__SCALING), -cov_U[l,l], cov_U[l,l]) * U__SCALING
  U__phi.star[U__phi.star < 0] <- U__phi.star[U__phi.star < 0] + U__SCALING[U__phi.star < 0] # put in interval
  U__phi.star[U__phi.star > U__SCALING] <- U__phi.star[U__phi.star > U__SCALING] - U__SCALING[U__phi.star > U__SCALING] # put in interval
  U.star <- U.old
  U.star[,,l] <- cholesky_UFromPhi(U__phi.star)
  if (hasEigenValueSmallerZero(U.star[,,l], TOL=NUMERICAL_THRESH)) { # stay positive definite
    if (verbose) print_warn(paste0("Discaring U_", l, " prosal",  #"with value ", U.star[,,l],
                                   " because of eigenvalues"))
    U[,,l,i+1] <- U[,,l,i] # reject and use previous
    U__phi[,l,i+1] <- U__phi[,l,i]
    rejectedU <- T
  }
  if (matCond(U.star[,,l]) < NUMERICAL_THRESH) { # stay numerically stable
    if (verbose) print_warn(paste0("Discaring U_", l, " prosal ", # "with value ", U.star[,,l],
                                   " because of numerics"))
    U[,,l,i+1] <- U[,,l,i] # reject and use previous
    U__phi[,l,i+1] <- U__phi[,l,i]
    rejectedU <- T
  }
  if (!rejectedU) {
    f.U.star <- lpost_matrixGamma(samp_freq=samp_freq,
                                  omega=omega,
                                  mpg_avg=mpg_avg,
                                  Nb=Nb,
                                  FZ=NULL,
                                  data=data,
                                  r=r[,i+1],
                                  U=U.star,
                                  Z=Z[,i+1],
                                  k=k[,i+1],
                                  C_alpha=C_alpha,
                                  omega_fun=omega_fun,
                                  k.theta=k.theta,
                                  db.list=db.list,
                                  bspline=bspline,
                                  degree_bs=degree_bs,
                                  V=V[,i],
                                  W=W[,i],
                                  H0.alpha=H0.alpha,
                                  H0.beta=H0.beta,
                                  MH=MH,
                                  eta=eta,
                                  Sigma_fun=Sigma_fun,
                                  corrected=corrected,
                                  toggle=toggle,
                                  f_param_avg_half=f_param_avg_half,
                                  phi=phi.fit,
                                  sigma_ar=sigma.fit,
                                  prior.q=prior.q,
                                  prior.cholesky=prior.cholesky,
                                  excludeBoundary=T, # note
                                  verbose=verbose)
    # accept/reject
    f.U <- f.store
    # Note: symmetric Uniform proposals -- but need jacobian of transformation [phi_U |-> U]
    alpha4 <- min(0, f.U.star +
                    cholesky_jacobianLogDeterminant(U__phi.star) -
                    f.U -
                    cholesky_jacobianLogDeterminant(U__phi[,l,i]))
    if (log(runif(1,0,1)) < alpha4) {
      U[,,l,i+1] <- U.star[,,l] # accept
      U__phi[,l,i+1] <- U__phi.star
      U.old <- U.star
      f.store <- f.U.star
    } else {
      U[,,l,i+1] <- U[,,l,i] # reject
      U__phi[,l,i+1] <- U__phi[,l,i]
    }
  }
}





