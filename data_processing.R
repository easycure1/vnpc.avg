library(Rcpp)
library(astsa)

sourceCpp("C:/Users/Yixuan/Documents/PhD/LISA/cpp_misc.cpp")

load('C:/Users/Yixuan/Documents/PhD/LISA/data/caseA_nlized.RData')


data <- caseA_nlized

data_new <- array(rep(NA, 32768*3*125), dim = c(32768, 3, 125))
for (ii in 1:125) {
  data_new[,,ii] <- data[(32768*(ii-1)+1):(32768*ii),]
}


mpg_avg <- array(0, dim = c(3, 3, 2049))
for (ii in 1:125) {
  mpg <- vnpc.avg:::mdft(data_new[,,ii])
  mpg <- vnpc.avg:::mpdgrm(mpg, freq = 2048)[,,1:2049]
  mpg_avg <- mpg_avg + mpg
}
mpg_avg <- mpg_avg/125



#---------------------------generate_time_series--------------------------------#


freq <- c(seq(0, pi, length=1969))
nf <- length(freq)
d <- dim(mpg_avg)[1]


N <- (nf - 1) * 2
S_full <- array(0+0i, dim = c(d,d,N))

S_full[,,1:nf] <- mpg_avg

for (k in 2:(nf-1)) {
  S_full[,,N-k+2] <- Conj(S_full[,,k])
}


X_time <- array(NA, dim = c(N, d, 125))

for (i in 1:125) {
  X_fft <- matrix(0+0i, nrow = N, ncol = d)
  
  for (k in 1:N) {
    # Cholesky factorization of S(f_k)
    Sk <- (S_full[,,k] + Conj(t(S_full[,,k]))) / 2
    Sk <- Re(Sk)
    Lk <- chol(Sk, pivot = FALSE)
    
    # Complex Gaussian noise ~ CN(0, I_d)
    z <- (rnorm(d) + 1i * rnorm(d)) / sqrt(2)
    
    # Fourier coefficient
    X_fft[k, ] <- Lk %*% z
  }
  
  for (j in 1:d) {
    X_time[,j,i] <- Re(fft(X_fft[,j], inverse = TRUE) / N)
  }
}


X_mpg_avg <- array(0, dim = c(3, 3, 1969))

for (i in 1:125) {
  X_mpg <- vnpc.avg:::mdft(X_time[,,i])
  X_mpg <- vnpc.avg:::mpdgrm(X_mpg, freq=2048)
  X_mpg_avg <- X_mpg_avg + X_mpg
}


save(X_time, file = 'C:/Users/Yixuan/Documents/PhD/LISA/data/segmented_caseA_generate.RData')




#------------------------------------VAR_spectral------------------------------#

f_param_avg <- array(0, dim = c(3,3,1969))
for (j in 1:125) {
  out <- ar(X_time[,,j], aic=F, order.max = 21)
  
  beta_start <- NULL
  for (jj in 1:3) {
    beta_start <- c(beta_start, c(t(out$ar[,jj,])))
  }
  param__beta <- beta_start
  param__phi <- vnpc.avg:::phiFromBeta_normalInverseWishart(param__beta, 3, 21)
  
  psd <- beyondWhittle::psd_varma(freq = index, ar = param__phi, Sigma = out$var.pred)
  
  f_param_avg <- f_param_avg + psd
}
f_param_avg <- f_param_avg/125


save(f_param_avg, file = 'C:/Users/Yixuan/Documents/PhD/LISA/data/f_param_avg_gr_21.RData')





# f_param_avg <- array(0, dim = c(3,3,1969))
# index <- seq(5*2*pi/2048, 128*2*pi/2048, length = 1969)
# 
# for (i in 1:125) {
#   out <- ar(data_new[,,i], aic=F, order.max = 246)
#   
#   beta_start <- NULL
#   for (jj in 1:3) {
#     beta_start <- c(beta_start, c(t(out$ar[,jj,])))
#   }
#   param__beta <- beta_start
#   param__phi <- vnpc.avg:::phiFromBeta_normalInverseWishart(param__beta, 3, 246)
#   
#   psd <- beyondWhittle::psd_varma(freq = index, ar = param__phi, Sigma = out$var.pred)
#   
#   f_param_avg <- f_param_avg + psd
# }
# 
# f_param_avg <- f_param_avg/125
# 
# save(f_param_avg, file = 'f_param_avg_246.RData')


# f_param_all <- array(rep(NA, 3*3*1969*125), dim = c(3, 3, 1969, 125))
# for (i in 1:125) {
#   spec_mat <- mvspec(data_new[,,1], spans = c(100,100), plot = F)
#   f_param_all[,,,i] <- spec_mat$fxx[,,(5*16):(128*16)]
# }
# f_param_avg <- array(rep(0, 3*3*1969), dim = c(3, 3, 1969))
# for (i in 1:125) {
#   f_param_avg <- f_param_avg + f_param_all[,,,i]
# }
# f_param_avg <- f_param_avg/125
# 
# save(f_param_avg, file = "estimated_param_psd_avg_100.RData")





# f_param_avg_vec <- rep(NA, 3*3*2*1969)
# for (i in 1:1969) {
#   f_param_avg_vec[((i-1)*18+1):((i-1)*18+1+2)] <- Re(f_param_avg[1,,i])
#   f_param_avg_vec[((i-1)*18+1+3):((i-1)*18+1+5)] <- Im(f_param_avg[1,,i])
#   f_param_avg_vec[((i-1)*18+1+6):((i-1)*18+1+8)] <- Re(f_param_avg[2,,i])
#   f_param_avg_vec[((i-1)*18+1+9):((i-1)*18+1+11)] <- Im(f_param_avg[2,,i])
#   f_param_avg_vec[((i-1)*18+1+12):((i-1)*18+1+14)] <- Re(f_param_avg[3,,i])
#   f_param_avg_vec[((i-1)*18+1+15):((i-1)*18+1+17)] <- Im(f_param_avg[3,,i])
# }
# 
# f_param_avg_new <- array(rep(NA, 3*3*1969), dim = c(3, 3, 1969))
# for (i in 1:1969) {
#   f_param_avg_new[1,,i] <- complex(real = f_param_avg_vec[((i-1)*18+1):((i-1)*18+1+2)], 
#                                    imaginary = f_param_avg_vec[((i-1)*18+1+3):((i-1)*18+1+5)])
#   f_param_avg_new[2,,i] <- complex(real = f_param_avg_vec[((i-1)*18+1+6):((i-1)*18+1+8)],
#                                    imaginary = f_param_avg_vec[((i-1)*18+1+9):((i-1)*18+1+11)])
#   f_param_avg_new[3,,i] <- complex(real = f_param_avg_vec[((i-1)*18+1+12):((i-1)*18+1+14)],
#                                    imaginary = f_param_avg_vec[((i-1)*18+1+15):((i-1)*18+1+17)])
# }
# 
# 
# llike_log <- function(f) {
#   load('C:/Users/Yixuan/Documents/PhD/LISA/data/averaged_periodogram.RData')
#   N <- 1969
#   f_new <- array(rep(NA, 3*3*1969), dim = c(3, 3, 1969))
#   for (i in 1:1969) {
#     f_new[1,,i] <- complex(real = f[((i-1)*18+1):((i-1)*18+1+2)], 
#                            imaginary = f[((i-1)*18+1+3):((i-1)*18+1+5)])
#     f_new[2,,i] <- complex(real = f[((i-1)*18+1+6):((i-1)*18+1+8)],
#                            imaginary = f[((i-1)*18+1+9):((i-1)*18+1+11)])
#     f_new[3,,i] <- complex(real = f[((i-1)*18+1+12):((i-1)*18+1+14)],
#                            imaginary = f[((i-1)*18+1+15):((i-1)*18+1+17)])
#   }
#   res <- llike_log_cpp(mpg_avg, f_new)
#   
#   return(res)
# }
# 
# out <- llike_log(f_param_avg_vec)
# 
# out <- optim(f_param_avg_vec, llike_log)
# 
# save(out, file = "param_optim.RData")
# 
# 
# f_param_optim <- array(rep(NA, 3*3*1969), dim = c(3, 3, 1969))
# for (i in 1:1969) {
#   f_param_optim[1,,i] <- complex(real = out$par[((i-1)*18+1):((i-1)*18+1+2)], 
#                                 imaginary = out$par[((i-1)*18+1+3):((i-1)*18+1+5)])
#   f_param_optim[2,,i] <- complex(real = out$par[((i-1)*18+1+6):((i-1)*18+1+8)],
#                                  imaginary = out$par[((i-1)*18+1+9):((i-1)*18+1+11)])
#   f_param_optim[3,,i] <- complex(real = out$par[((i-1)*18+1+12):((i-1)*18+1+14)],
#                          imaginary = out$par[((i-1)*18+1+15):((i-1)*18+1+17)])
# }
# 
# 


