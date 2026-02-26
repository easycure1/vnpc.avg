#-----------------plots for 1 yr simulated data---------------------------#

library(reticulate)
library(vnpc.avg)
library(ggplot2)
library(patchwork)
library(scales)
library(signal)

np <- import("numpy")
data_npy <- np$load("C:/Users/Yixuan/Documents/PhD/LISA/data/LISA_data/lisa_simulated_noise_dataset_1yr/new/data.npy")
freq_true <- np$load("C:/Users/Yixuan/Documents/PhD/LISA/data/LISA_data/lisa_simulated_noise_dataset_1yr/new/freq_true.npy")
# S_true <- np$load("C:/Users/Yixuan/Documents/PhD/LISA/data/LISA_data/lisa_simulated_noise_dataset_1yr/new/S_true_est.npy")
# S_welch <- np$load("C:/Users/Yixuan/Documents/PhD/LISA/data/LISA_data/lisa_simulated_noise_dataset_1yr/dt_5/S_emp.npy", allow_pickle=T)
# true_matrix_freq <- np$load("C:/Users/Yixuan/Documents/PhD/LISA/data/LISA_data/lisa_simulated_noise_dataset_1yr/dt_5/true_matrix_freq.npy")
true_matrix <- np$load("C:/Users/Yixuan/Documents/PhD/LISA/data/LISA_data/lisa_simulated_noise_dataset_1yr/new/true_matrix.npy")
dt <- np$load("C:/Users/Yixuan/Documents/PhD/LISA/data/LISA_data/lisa_simulated_noise_dataset_1yr/new/delta_t.npy")



# pick points evenly spaced on log scale of index
m <- 10000  # target number of points to plot
idx <- unique(round(exp(seq(log(1), log(length(freq_true)), length.out = m)))) # create indices that are evenly spaced on a log scale

# true_matrix_freq_sub <- true_matrix_freq[idx,,]
# freq_true_sub <- freq_true[idx]

true_matrix_sub <- true_matrix[idx,,]
freq_true_sub <- freq_true[idx]

# freq_S <- seq(freq_true_sub[1], 0.1, len = 1001)

freq_vnp <- seq(0, .1, len = 2049)
# scales <- freq_vnp^2*2/(2000*mean(hanning(2000)^2))/2/pi
# scales <- freq_vnp^2/(4097*mean(hanning(4097)^2))/2/pi


# plot(freq_true_sub[1083:3391], true_matrix_freq_sub[1083:3391,1,1], type = "l", log = "xy")
# lines(freq_S[-1], S_true[-1,1,1], col='red')
# lines(freq_S[-1], S_welch[[1]]$Sxx[-1], col='blue')



coh <- function(Sii, Sjj, Sij) abs(Sij) / sqrt(Sii * Sjj)

coh_mcmc <- function(Sii, Sjj, Sij, Sji) sqrt(Sij^2+Sji^2) / sqrt(Sii * Sjj)

true_psd_gg <- data.frame(
  fx = true_matrix_sub[1634:6321,1,1]/sd(true_matrix_sub[1634:6321,1,1]),
  fy = true_matrix_sub[1634:6321,2,2]/sd(true_matrix_sub[1634:6321,2,2]),
  fz = true_matrix_sub[1634:6321,3,3]/sd(true_matrix_sub[1634:6321,3,3]),
  fxy = true_matrix_sub[1634:6321,1,2]/sd(true_matrix_sub[1634:6321,1,2]),
  fxz = true_matrix_sub[1634:6321,1,3]/sd(true_matrix_sub[1634:6321,1,3]),
  fyz = true_matrix_sub[1634:6321,2,3]/sd(true_matrix_sub[1634:6321,2,3]),
  fyx = true_matrix_sub[1634:6321,2,1]/sd(true_matrix_sub[1634:6321,2,1]),
  fzx = true_matrix_sub[1634:6321,3,1]/sd(true_matrix_sub[1634:6321,3,1]),
  fzy = true_matrix_sub[1634:6321,3,2]/sd(true_matrix_sub[1634:6321,3,2]),
  cohxy = coh(true_matrix_sub[1634:6321,1,1],
              true_matrix_sub[1634:6321,2,2],
              true_matrix_sub[1634:6321,1,2]),
  cohxz = coh(true_matrix_sub[1634:6321,1,1],
              true_matrix_sub[1634:6321,3,3],
              true_matrix_sub[1634:6321,1,3]),
  cohyz = coh(true_matrix_sub[1634:6321,2,2],
              true_matrix_sub[1634:6321,3,3],
              true_matrix_sub[1634:6321,2,3]),
  freq = freq_true_sub[1634:6321]
)

# S_true_gg <- data.frame(
#   fx = S_true[-1,1,1],
#   fy = S_true[-1,2,2],
#   fz = S_true[-1,3,3],
#   fxy = S_true[-1,1,2],
#   fxz = S_true[-1,1,3],
#   fyz = S_true[-1,2,3],
#   fyx = S_true[-1,2,1],
#   fzx = S_true[-1,3,1],
#   fzy = S_true[-1,3,2],
#   cohxy = coh(S_true[-1,1,1],
#               S_true[-1,2,2],
#               S_true[-1,1,2]),
#   cohxz = coh(S_true[-1,1,1],
#               S_true[-1,3,3],
#               S_true[-1,1,3]),
#   cohyz = coh(S_true[-1,2,2],
#               S_true[-1,3,3],
#               S_true[-1,2,3]),
#   freq = freq_S[-1]
# )
# 
# welch_gg <- data.frame(
#   fx = S_welch[[1]]$Sxx[-1],
#   fy = S_welch[[1]]$Syy[-1],
#   fz = S_welch[[1]]$Szz[-1],
#   fxy = Re(S_welch[[1]]$Sxy[-1]),
#   fxz = Re(S_welch[[1]]$Szx[-1]),
#   fyz = Re(S_welch[[1]]$Szz[-1]),
#   fyx = Im(S_welch[[1]]$Sxy[-1]),
#   fzx = Im(S_welch[[1]]$Szx[-1]),
#   fzy = Im(S_welch[[1]]$Szz[-1]),
#   cohxy = coh(S_welch[[1]]$Sxx[-1],
#               S_welch[[1]]$Syy[-1],
#               S_welch[[1]]$Sxy[-1]),
#   cohxz = coh(S_welch[[1]]$Sxx[-1],
#               S_welch[[1]]$Szz[-1],
#               S_welch[[1]]$Szx[-1]),
#   cohyz = coh(S_welch[[1]]$Syy[-1],
#               S_welch[[1]]$Szz[-1],
#               S_welch[[1]]$Syz[-1]),
#   freq = freq_S[-1]
# )


hann_window <- hanning(4096)

load("C:/Users/Yixuan/Documents/PhD/LISA/results/original/simulated_data/1yr/mcmc_vnp_avg_lisa_simu_1yr_new_4096.RData")

vnp_gg <- data.frame(
  fx = mcmc_vnp_avg$psd.median[1,1,3:2049] * mean(hann_window^2)*2,
  fy = mcmc_vnp_avg$psd.median[2,2,3:2049] * mean(hann_window^2)*2,
  fz = mcmc_vnp_avg$psd.median[3,3,3:2049] * mean(hann_window^2)*2,
  fxy = mcmc_vnp_avg$psd.median[1,2,3:2049] * mean(hann_window^2)*2,
  fxz = mcmc_vnp_avg$psd.median[1,3,3:2049] * mean(hann_window^2)*2,
  fyz = mcmc_vnp_avg$psd.median[2,3,3:2049] * mean(hann_window^2)*2,
  fyx = mcmc_vnp_avg$psd.median[2,1,3:2049] * mean(hann_window^2)*2,
  fzx = mcmc_vnp_avg$psd.median[3,1,3:2049] * mean(hann_window^2)*2,
  fzy = mcmc_vnp_avg$psd.median[3,2,3:2049] * mean(hann_window^2)*2,
  cohxy = coh(mcmc_vnp_avg$psd.median[1,1,3:2049],
              mcmc_vnp_avg$psd.median[2,2,3:2049],
              mcmc_vnp_avg$psd.median[1,2,3:2049]),
  cohxz = coh(mcmc_vnp_avg$psd.median[1,1,3:2049],
              mcmc_vnp_avg$psd.median[3,3,3:2049],
              mcmc_vnp_avg$psd.median[1,3,3:2049]),
  cohyz = coh(mcmc_vnp_avg$psd.median[2,2,3:2049],
              mcmc_vnp_avg$psd.median[3,3,3:2049],
              mcmc_vnp_avg$psd.median[2,3,3:2049]),
  freq = freq_vnp[3:2049]
)


load("C:/Users/Yixuan/Documents/PhD/LISA/results/original/simulated_data/1yr/mcmc_vnpc_avg_order_200_lisa_simu_1yr_4096_trunc_new.RData")

vnpc_200_gg <- data.frame(
  fx = mcmc_vnpc_avg$psd.median[1,1,3:2049] * mean(hann_window^2)*2,
  fy = mcmc_vnpc_avg$psd.median[2,2,3:2049] * mean(hann_window^2)*2,
  fz = mcmc_vnpc_avg$psd.median[3,3,3:2049] * mean(hann_window^2)*2,
  fxy = mcmc_vnpc_avg$psd.median[1,2,3:2049] * mean(hann_window^2)*2,
  fxz = mcmc_vnpc_avg$psd.median[1,3,3:2049] * mean(hann_window^2)*2,
  fyz = mcmc_vnpc_avg$psd.median[2,3,3:2049] * mean(hann_window^2)*2,
  fyx = mcmc_vnpc_avg$psd.median[2,1,3:2049] * mean(hann_window^2)*2,
  fzx = mcmc_vnpc_avg$psd.median[3,1,3:2049] * mean(hann_window^2)*2,
  fzy = mcmc_vnpc_avg$psd.median[3,2,3:2049] * mean(hann_window^2)*2,
  cohxy = coh(mcmc_vnpc_avg$psd.median[1,1,3:2049],
              mcmc_vnpc_avg$psd.median[2,2,3:2049],
              mcmc_vnpc_avg$psd.median[1,2,3:2049]),
  cohxz = coh(mcmc_vnpc_avg$psd.median[1,1,3:2049],
              mcmc_vnpc_avg$psd.median[3,3,3:2049],
              mcmc_vnpc_avg$psd.median[1,3,3:2049]),
  cohyz = coh(mcmc_vnpc_avg$psd.median[2,2,3:2049],
              mcmc_vnpc_avg$psd.median[3,3,3:2049],
              mcmc_vnpc_avg$psd.median[2,3,3:2049]),
  freq = freq_vnp[3:2049]
)


load("C:/Users/Yixuan/Documents/PhD/LISA/results/original/simulated_data/1yr/mcmc_vnpc_avg_order_300_lisa_simu_1yr_4096_trunc_new.RData")

vnpc_300_gg <- data.frame(
  fx = mcmc_vnpc_avg$psd.median[1,1,3:2049] * mean(hann_window^2)*2,
  fy = mcmc_vnpc_avg$psd.median[2,2,3:2049] * mean(hann_window^2)*2,
  fz = mcmc_vnpc_avg$psd.median[3,3,3:2049] * mean(hann_window^2)*2,
  fxy = mcmc_vnpc_avg$psd.median[1,2,3:2049] * mean(hann_window^2)*2,
  fxz = mcmc_vnpc_avg$psd.median[1,3,3:2049] * mean(hann_window^2)*2,
  fyz = mcmc_vnpc_avg$psd.median[2,3,3:2049] * mean(hann_window^2)*2,
  fyx = mcmc_vnpc_avg$psd.median[2,1,3:2049] * mean(hann_window^2)*2,
  fzx = mcmc_vnpc_avg$psd.median[3,1,3:2049] * mean(hann_window^2)*2,
  fzy = mcmc_vnpc_avg$psd.median[3,2,3:2049] * mean(hann_window^2)*2,
  cohxy = coh(mcmc_vnpc_avg$psd.median[1,1,3:2049],
              mcmc_vnpc_avg$psd.median[2,2,3:2049],
              mcmc_vnpc_avg$psd.median[1,2,3:2049]),
  cohxz = coh(mcmc_vnpc_avg$psd.median[1,1,3:2049],
              mcmc_vnpc_avg$psd.median[3,3,3:2049],
              mcmc_vnpc_avg$psd.median[1,3,3:2049]),
  cohyz = coh(mcmc_vnpc_avg$psd.median[2,2,3:2049],
              mcmc_vnpc_avg$psd.median[3,3,3:2049],
              mcmc_vnpc_avg$psd.median[2,3,3:2049]),
  freq = freq_vnp[3:2049]
)


load("C:/Users/Yixuan/Documents/PhD/LISA/results/original/simulated_data/1yr/mcmc_vnpc_avg_order_49_lisa_simu_1yr_8193.RData")

vnpc_49_gg <- data.frame(
  fx = mcmc_vnpc_avg$psd.median[1,1,3:2049]*scales,
  fy = mcmc_vnpc_avg$psd.median[2,2,3:2049]*scales,
  fz = mcmc_vnpc_avg$psd.median[3,3,3:2049]*scales,
  fxy = mcmc_vnpc_avg$psd.median[1,2,3:2049]*scales,
  fxz = mcmc_vnpc_avg$psd.median[1,3,3:2049]*scales,
  fyz = mcmc_vnpc_avg$psd.median[2,3,3:2049]*scales,
  fyx = mcmc_vnpc_avg$psd.median[2,1,3:2049]*scales,
  fzx = mcmc_vnpc_avg$psd.median[3,1,3:2049]*scales,
  fzy = mcmc_vnpc_avg$psd.median[3,2,3:2049]*scales,
  cohxy = coh(mcmc_vnpc_avg$psd.median[1,1,3:2049],
              mcmc_vnpc_avg$psd.median[2,2,3:2049],
              mcmc_vnpc_avg$psd.median[1,2,3:2049]),
  cohxz = coh(mcmc_vnpc_avg$psd.median[1,1,3:2049],
              mcmc_vnpc_avg$psd.median[3,3,3:2049],
              mcmc_vnpc_avg$psd.median[1,3,3:2049]),
  cohyz = coh(mcmc_vnpc_avg$psd.median[2,2,3:2049],
              mcmc_vnpc_avg$psd.median[3,3,3:2049],
              mcmc_vnpc_avg$psd.median[2,3,3:2049]),
  freq = freq_vnp
)


load("C:/Users/Yixuan/Documents/PhD/LISA/results/original/simulated_data/1yr/mcmc_vnpc_avg_order_99_lisa_simu_1yr_8193.RData")

vnpc_99_gg <- data.frame(
  fx = mcmc_vnpc_avg$psd.median[1,1,3:2049]*scales,
  fy = mcmc_vnpc_avg$psd.median[2,2,3:2049]*scales,
  fz = mcmc_vnpc_avg$psd.median[3,3,3:2049]*scales,
  fxy = mcmc_vnpc_avg$psd.median[1,2,3:2049]*scales,
  fxz = mcmc_vnpc_avg$psd.median[1,3,3:2049]*scales,
  fyz = mcmc_vnpc_avg$psd.median[2,3,3:2049]*scales,
  fyx = mcmc_vnpc_avg$psd.median[2,1,3:2049]*scales,
  fzx = mcmc_vnpc_avg$psd.median[3,1,3:2049]*scales,
  fzy = mcmc_vnpc_avg$psd.median[3,2,3:2049]*scales,
  cohxy = coh(mcmc_vnpc_avg$psd.median[1,1,3:2049],
              mcmc_vnpc_avg$psd.median[2,2,3:2049],
              mcmc_vnpc_avg$psd.median[1,2,3:2049]),
  cohxz = coh(mcmc_vnpc_avg$psd.median[1,1,3:2049],
              mcmc_vnpc_avg$psd.median[3,3,3:2049],
              mcmc_vnpc_avg$psd.median[1,3,3:2049]),
  cohyz = coh(mcmc_vnpc_avg$psd.median[2,2,3:2049],
              mcmc_vnpc_avg$psd.median[3,3,3:2049],
              mcmc_vnpc_avg$psd.median[2,3,3:2049]),
  freq = freq_vnp
)


load("C:/Users/Yixuan/Documents/PhD/LISA/data/LISA_data/lisa_simulated_noise_dataset_1yr/new/mpgm_avg_4096.RData")

avg_mpg_gg <- data.frame(
  fx = Re(mpgm_avg[1,1,3:2049]) * mean(hann_window^2)*2,
  fy = Re(mpgm_avg[2,2,3:2049]) * mean(hann_window^2)*2,
  fz = Re(mpgm_avg[3,3,3:2049]) * mean(hann_window^2)*2,
  fxy = Re(mpgm_avg[1,2,3:2049]) * mean(hann_window^2)*2,
  fxz = Re(mpgm_avg[1,3,3:2049]) * mean(hann_window^2)*2,
  fyz = Re(mpgm_avg[2,3,3:2049]) * mean(hann_window^2)*2,
  fyx = Im(mpgm_avg[2,1,3:2049]) * mean(hann_window^2)*2,
  fzx = Im(mpgm_avg[3,1,3:2049]) * mean(hann_window^2)*2,
  fzy = Im(mpgm_avg[3,2,3:2049]) * mean(hann_window^2)*2,
  freq = freq_vnp[3:2049]
)


#---------------------elbow criterion----------------------------#

load("C:/Users/Yixuan/Documents/PhD/LISA/data/LISA_data/lisa_simulated_noise_dataset_1yr/var_order_8193.RData")

plot_order <- ggplot(data.frame(order_fit[1:101]), aes(3:103, order_fit[1:101])) +
  geom_line(linewidth = 1) +
  geom_point(colour = 'red3', size = 2) +
  geom_vline(xintercept = seq(3,103,by=2), color = 'blue', alpha = .3) +
  scale_x_continuous(breaks = seq(3, 103, by = 4)) +
  scale_y_continuous(labels = NULL) +
  xlab("lag") +
  ylab("-log-likelihood") +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15))


pdf(file = "C:/Users/Yixuan/Documents/PhD/LISA/results/original/simulated_data/1yr/fit_order_1yr_8193_plot.pdf",
    width = 10, height = 6)

plot_order

dev.off()






plot_xx <- ggplot(true_psd_gg) +
  geom_line(aes(x = freq, y = fx, colour = "true psd"), linewidth = 1) +
  # geom_line(data = S_true_gg, aes(x = freq, y = fx, colour = "S_true"),
  #           linewidth = .8) +
  # geom_line(data = welch_gg, aes(x = freq, y = fx, colour = "Welch psd"),
  #           linewidth = .8) +
  geom_line(data = vnp_gg, aes(x = freq, y = fx, colour = "vnp"),
            linewidth = .8) +
  geom_line(data = vnpc_200_gg, aes(x = freq, y = fx, colour = "vnpc(200)"),
            linewidth = .8) +
  geom_line(data = avg_mpg_gg, aes(x = freq, y = fx, colour = "average periodogram"),
            linewidth = .8) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_colour_manual(values = c('true psd' = 'black', #'S_true' = 'red3', 'Welch psd' = 'cyan3', 
                                 'vnp' = "red3",
                                 'vnpc(200)' = "cyan3",  'average periodogram' = "slategrey")) +
  labs(x = NULL, y = "PSD", title = expression(S[XX])) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none')


plot_yy <- ggplot(true_psd_gg) +
  geom_line(aes(x = freq, y = fy, colour = "true psd"), linewidth = 1) +
  # geom_line(data = S_true_gg, aes(x = freq, y = fy, colour = "S_true"),
  #           linewidth = .8) +
  # geom_line(data = welch_gg, aes(x = freq, y = fy, colour = "Welch psd"),
  #           linewidth = .8) +
  geom_line(data = vnp_gg, aes(x = freq, y = fy, colour = "vnp"),
            linewidth = .8) +
  geom_line(data = vnpc_200_gg, aes(x = freq, y = fy, colour = "vnpc(200)"),
            linewidth = .8) +
  geom_line(data = avg_mpg_gg, aes(x = freq, y = fy, colour = "average periodogram"),
            linewidth = .8) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_colour_manual(values = c('true psd' = 'black', #'S_true' = 'red3', 'Welch psd' = 'cyan3', 
                                 'vnp' = "red3", 'vnpc(200)' = "cyan3",  
                                 'average periodogram' = "slategrey")) +
  labs(x = NULL, y = NULL, title = expression(S[YY])) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none')


plot_zz <- ggplot(true_psd_gg) +
  geom_line(aes(x = freq, y = fz, colour = "true psd"), linewidth = 1) +
  # geom_line(data = S_true_gg, aes(x = freq, y = fz, colour = "S_true"),
  #           linewidth = .8) +
  # geom_line(data = welch_gg, aes(x = freq, y = fz, colour = "Welch psd"),
  #           linewidth = .8) +
  geom_line(data = vnp_gg, aes(x = freq, y = fz, colour = "vnp"),
            linewidth = .8) +
  geom_line(data = vnpc_200_gg, aes(x = freq, y = fz, colour = "vnpc(200)"),
            linewidth = .8) +
  geom_line(data = avg_mpg_gg, aes(x = freq, y = fz, colour = "average periodogram"),
            linewidth = .8) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_colour_manual(values = c('true psd' = 'black', # 'S_true' = 'red3', 'Welch psd' = 'cyan3',
                                 'vnp' = "red3", 'vnpc(200)' = "cyan3",  
                                 'average periodogram' = "slategrey")) +
  labs(x = NULL, y = NULL, title = expression(S[ZZ]), colour = "Method") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20))


plot_xy <- ggplot(true_psd_gg) +
  geom_line(aes(x = freq, y = fxy, colour = "true psd"), linewidth = 1) +
  # geom_line(data = S_true_gg, aes(x = freq, y = fxy, colour = "S_true"),
  #           linewidth = 0.8) +
  # geom_line(data = welch_gg, aes(x = freq, y = fxy, colour = "Welch"),
  #           linewidth = 0.8) +
  geom_line(data = vnp_gg, aes(x = freq, y = fxy, colour = "vnp"),
            linewidth = 0.8) +
  geom_line(data = vnpc_200_gg, aes(x = freq, y = fxy, colour = "vnpc(200)"),
            linewidth = .8) +
  geom_line(data = avg_mpg_gg, aes(x = freq, y = fxy, colour = "average periodogram"),
            linewidth = .8) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_colour_manual(values = c('true psd' = 'black', #'S_true' = 'red3', 'Welch psd' = 'cyan3',
                                 'vnp' = "red3", 'vnpc(200)' = "cyan3",  
                                 'average periodogram' = "slategrey")) +
  labs(x = NULL, y = NULL, title = expression(Re(S[XY]))) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none')


plot_xz <- ggplot(true_psd_gg) +
  geom_line(aes(x = freq, y = fxz, colour = "true psd"), linewidth = 1) +
  # geom_line(data = S_true_gg, aes(x = freq, y = fxz, colour = "S_true"),
  #           linewidth = 0.8) +
  # geom_line(data = welch_gg, aes(x = freq, y = fxz, colour = "Welch"),
  #           linewidth = 0.8) +
  geom_line(data = vnp_gg, aes(x = freq, y = fxz, colour = "vnp"),
            linewidth = 0.8) +
  geom_line(data = vnpc_200_gg, aes(x = freq, y = fxz, colour = "vnpc(200)"),
            linewidth = .8) +
  geom_line(data = avg_mpg_gg, aes(x = freq, y = fxz, colour = "average periodogram"),
            linewidth = .8) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_colour_manual(values = c('true psd' = 'black', #'S_true' = 'red3', 'Welch psd' = 'cyan3', 
                                 'vnp' = "red3", 'vnpc(200)' = "cyan3",
                                 'average periodogram' = "slategrey")) +
  labs(x = NULL, y = NULL, title = expression(Re(S[XZ]))) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none')


plot_yz <- ggplot(true_psd_gg) +
  geom_line(aes(x = freq, y = fyz, colour = "true psd"), linewidth = 1) +
  # geom_line(data = S_true_gg, aes(x = freq, y = fyz, colour = "S_true"),
  #           linewidth = 0.8) +
  # geom_line(data = welch_gg, aes(x = freq, y = fyz, colour = "Welch"),
  #           linewidth = 0.8) +
  geom_line(data = vnp_gg, aes(x = freq, y = fyz, colour = "vnp"),
            linewidth = 0.8) +
  geom_line(data = vnpc_200_gg, aes(x = freq, y = fyz, colour = "vnpc(200)"),
            linewidth = .8) +
  geom_line(data = avg_mpg_gg, aes(x = freq, y = fyz, colour = "average periodogram"),
            linewidth = .8) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_colour_manual(values = c('true psd' = 'black', #'S_true' = 'red3', 'Welch psd' = 'cyan3', 
                                 'vnp' = "red3", 'vnpc(200)' = "cyan3",
                                 'average periodogram' = "slategrey")) +
  labs(x = NULL, y = NULL, title = expression(Re(S[YZ]))) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none')


plot_yx <- ggplot(avg_mpg_gg) +
  geom_line(aes(x = freq, y = fyx, colour = "average periodogram")) +
  # geom_line(data = welch_gg, aes(x = freq, y = fyz, colour = "Welch"),
  #           linewidth = 0.8) +
  geom_line(data = vnp_gg, aes(x = freq, y = fyx, colour = "vnp"),
            linewidth = 0.8) +
  geom_line(data = vnpc_200_gg, aes(x = freq, y = fyx, colour = "vnpc(200)"),
            linewidth = .8) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_colour_manual(values = c('vnp' = "red3", #'Welch psd' = 'cyan3',
                                 'vnpc(200)' = "cyan3", 
                                 'average periodogram' = "slategrey")) +
  labs(x = NULL, y = "PSD", title = expression(Im(S[YX]))) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none')


plot_zx <- ggplot(avg_mpg_gg) +
  geom_line(aes(x = freq, y = fzx, colour = "average periodogram")) +
  geom_line(data = vnp_gg, aes(x = freq, y = fzx, colour = "vnp"),
            linewidth = 0.8) +
  geom_line(data = vnpc_200_gg, aes(x = freq, y = fzx, colour = "vnpc(200)"),
            linewidth = .8) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_colour_manual(values = c('vnp' = "red3", 'vnpc(200)' = "cyan3", 
                                 'average periodogram' = "slategrey")) +
  labs(x = NULL, y = "PSD", title = expression(Im(S[ZX]))) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none')


plot_zy <- ggplot(avg_mpg_gg) +
  geom_line(aes(x = freq, y = fzy, colour = "average periodogram")) +
  geom_line(data = vnp_gg, aes(x = freq, y = fzy, colour = "vnp"),
            linewidth = 0.8) +
  geom_line(data = vnpc_200_gg, aes(x = freq, y = fzy, colour = "vnpc(200)"),
            linewidth = .8) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_colour_manual(values = c('vnp' = "red3", 'vnpc(200)' = "cyan3", 
                                 'average periodogram' = "slategrey")) +
  labs(x = NULL, y = NULL, title = expression(Im(S[ZY]))) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none')




pdf(file = "C:/Users/Yixuan/Documents/PhD/LISA/results/original/simulated_data/1yr/PSD_1yr_4096_new_plot.pdf",
    width = 18, height = 10)

(plot_xx | plot_xy | plot_xz) /
  (plot_yx | plot_yy | plot_yz) /
  (plot_zx | plot_zy | plot_zz)

dev.off()






plot_coh_xy <- ggplot(true_psd_gg) +
  geom_line(aes(freq, y = cohxy, colour = "true psd"), linewidth = 1) +
  # geom_line(data = S_true_gg, aes(freq, y = cohxy, colour = "S_true")) +
  # geom_line(data = welch_gg, aes(freq, y = cohxy, colour = "Welch psd")) +
  geom_line(data = vnp_gg, aes(freq, y = cohxy, colour = "vnp")) +
  geom_line(data = vnpc_200_gg, aes(freq, y = cohxy, colour = "vnpc(200)")) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_colour_manual(values = c("true psd" = "black", #"S_true" = "red3", "Welch psd" = "cyan3", 
                                 "vnp" = "red3", "vnpc(200)" = "cyan3")) +
  labs(x = NULL, y = "Coherence", title = expression(Coherence[XY])) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none')


plot_coh_xz <- ggplot(true_psd_gg) +
  geom_line(aes(freq, y = cohxz, colour = "true psd"), linewidth = 1) +
  # geom_line(data = S_true_gg, aes(freq, y = cohxz, colour = "S_true")) +
  # geom_line(data = welch_gg, aes(freq, y = cohxz, colour = "Welch psd")) +
  geom_line(data = vnp_gg, aes(freq, y = cohxz, colour = "vnp")) +
  geom_line(data = vnpc_200_gg, aes(freq, y = cohxz, colour = "vnpc(200)")) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_colour_manual(values = c("true psd" = "black", #"S_true" = "red3", "Welch psd" = "cyan3",
                                 "vnp" = "red3", "vnpc(200)" = "cyan3")) +
  labs(x = NULL, y = NULL, title = expression(Coherence[XZ])) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none')


plot_coh_yz <- ggplot(true_psd_gg) +
  geom_line(aes(freq, y = cohyz, colour = "true psd"), linewidth = 1) +
  # geom_line(data = S_true_gg, aes(freq, y = cohyz, colour = "S_true")) +
  # geom_line(data = welch_gg, aes(freq, y = cohyz, colour = "Welch psd")) +
  geom_line(data = vnp_gg, aes(freq, y = cohyz, colour = "vnp")) +
  geom_line(data = vnpc_200_gg, aes(freq, y = cohyz, colour = "vnpc(200)")) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_colour_manual(values = c("true psd" = "black", #"S_true" = "red3", "Welch psd" = "cyan3", 
                                 "vnp" = "red3", "vnpc(200)" = "cyan3")) +
  labs(x = NULL, y = NULL, title = expression(Coherence[YZ]), colour = "Method") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20))


pdf(file = "C:/Users/Yixuan/Documents/PhD/LISA/results/original/simulated_data/1yr/Coherence_1yr_4096_new_plot.pdf",
    width = 18, height = 10)

plot_coh_xy | plot_coh_xz | plot_coh_yz

dev.off()









#----------------------VNPC vs VNP------------------------------#

plot_xx <- ggplot(avg_mpg_gg) +
  geom_line(aes(x = freq, y = fx, colour = "average periodogram")) +
  geom_line(data = vnp_gg, aes(x = freq, y = fx, colour = "vnp"),
            linewidth = .8) +
  geom_line(data = vnpc_49_gg, aes(x = freq, y = fx, colour = "vnpc(49)"),
            linewidth = .8) +
  geom_line(data = vnpc_99_gg, aes(x = freq, y = fx, colour = "vnpc(99)"),
            linewidth = .8) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_colour_manual(values = c('average periodogram' = 'black', 'vnp' = "red3",
                                 'vnpc(49)' = "cyan3",  'vnpc(99)' = "olivedrab3")) +
  labs(x = NULL, y = "PSD", title = expression(S[XX])) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none')


plot_yy <- ggplot(avg_mpg_gg) +
  geom_line(aes(x = freq, y = fy, colour = "average periodogram")) +
  geom_line(data = vnp_gg, aes(x = freq, y = fy, colour = "vnp"),
            linewidth = .8) +
  geom_line(data = vnpc_49_gg, aes(x = freq, y = fy, colour = "vnpc(49)"),
            linewidth = .8) +
  geom_line(data = vnpc_99_gg, aes(x = freq, y = fy, colour = "vnpc(99)"),
            linewidth = .8) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_colour_manual(values = c('average periodogram' = 'black', 'vnp' = "red3",
                                 'vnpc(49)' = "cyan3",  'vnpc(99)' = "olivedrab3")) +
  labs(x = NULL, y = NULL, title = expression(S[YY])) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none')


plot_zz <- ggplot(avg_mpg_gg) +
  geom_line(aes(x = freq, y = fz, colour = "average periodogram")) +
  geom_line(data = vnp_gg, aes(x = freq, y = fz, colour = "vnp"),
            linewidth = .8) +
  geom_line(data = vnpc_49_gg, aes(x = freq, y = fz, colour = "vnpc(49)"),
            linewidth = .8) +
  geom_line(data = vnpc_99_gg, aes(x = freq, y = fz, colour = "vnpc(99)"),
            linewidth = .8) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_colour_manual(values = c('average periodogram' = 'black', 'vnp' = "red3",
                                 'vnpc(49)' = "cyan3",  'vnpc(99)' = "olivedrab3")) +
  labs(x = NULL, y = NULL, title = expression(S[ZZ])) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20))


plot_xy <- ggplot(avg_mpg_gg) +
  geom_line(aes(x = freq, y = fxy, colour = "average periodogram")) +
  geom_line(data = vnp_gg, aes(x = freq, y = fxy, colour = "vnp"),
            linewidth = .8) +
  geom_line(data = vnpc_49_gg, aes(x = freq, y = fxy, colour = "vnpc(49)"),
            linewidth = .8) +
  geom_line(data = vnpc_99_gg, aes(x = freq, y = fxy, colour = "vnpc(99)"),
            linewidth = .8) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_colour_manual(values = c('average periodogram' = 'black', 'vnp' = "red3",
                                 'vnpc(49)' = "cyan3",  'vnpc(99)' = "olivedrab3")) +
  labs(x = NULL, y = NULL, title = expression(Re(S[XY]))) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none')


plot_xz <- ggplot(avg_mpg_gg) +
  geom_line(aes(x = freq, y = fxz, colour = "average periodogram")) +
  geom_line(data = vnp_gg, aes(x = freq, y = fxz, colour = "vnp"),
            linewidth = .8) +
  geom_line(data = vnpc_49_gg, aes(x = freq, y = fxz, colour = "vnpc(49)"),
            linewidth = .8) +
  geom_line(data = vnpc_99_gg, aes(x = freq, y = fxz, colour = "vnpc(99)"),
            linewidth = .8) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_colour_manual(values = c('average periodogram' = 'black', 'vnp' = "red3",
                                 'vnpc(49)' = "cyan3",  'vnpc(99)' = "olivedrab3")) +
  labs(x = NULL, y = NULL, title = expression(Re(S[XZ]))) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20))


plot_yz <- ggplot(avg_mpg_gg) +
  geom_line(aes(x = freq, y = fyz, colour = "average periodogram")) +
  geom_line(data = vnp_gg, aes(x = freq, y = fyz, colour = "vnp"),
            linewidth = .8) +
  geom_line(data = vnpc_49_gg, aes(x = freq, y = fyz, colour = "vnpc(49)"),
            linewidth = .8) +
  geom_line(data = vnpc_99_gg, aes(x = freq, y = fyz, colour = "vnpc(99)"),
            linewidth = .8) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_colour_manual(values = c('average periodogram' = 'black', 'vnp' = "red3",
                                 'vnpc(49)' = "cyan3",  'vnpc(99)' = "olivedrab3")) +
  labs(x = NULL, y = NULL, title = expression(Re(S[YZ]))) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20))


plot_yx <- ggplot(avg_mpg_gg) +
  geom_line(aes(x = freq, y = fyx, colour = "average periodogram")) +
  geom_line(data = vnp_gg, aes(x = freq, y = fyx, colour = "vnp"),
            linewidth = .8) +
  geom_line(data = vnpc_49_gg, aes(x = freq, y = fyx, colour = "vnpc(49)"),
            linewidth = .8) +
  geom_line(data = vnpc_99_gg, aes(x = freq, y = fyx, colour = "vnpc(99)"),
            linewidth = .8) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_colour_manual(values = c('average periodogram' = 'black', 'vnp' = "red3",
                                 'vnpc(49)' = "cyan3",  'vnpc(99)' = "olivedrab3")) +
  labs(x = NULL, y = "PSD", title = expression(Im(S[YX]))) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none')


plot_zx <- ggplot(avg_mpg_gg) +
  geom_line(aes(x = freq, y = fzx, colour = "average periodogram")) +
  geom_line(data = vnp_gg, aes(x = freq, y = fzx, colour = "vnp"),
            linewidth = .8) +
  geom_line(data = vnpc_49_gg, aes(x = freq, y = fzx, colour = "vnpc(49)"),
            linewidth = .8) +
  geom_line(data = vnpc_99_gg, aes(x = freq, y = fzx, colour = "vnpc(99)"),
            linewidth = .8) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_colour_manual(values = c('average periodogram' = 'black', 'vnp' = "red3",
                                 'vnpc(49)' = "cyan3",  'vnpc(99)' = "olivedrab3")) +
  labs(x = NULL, y = "PSD", title = expression(Im(S[ZX]))) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none')


plot_zy <- ggplot(avg_mpg_gg) +
  geom_line(aes(x = freq, y = fzy, colour = "average periodogram")) +
  geom_line(data = vnp_gg, aes(x = freq, y = fzy, colour = "vnp"),
            linewidth = .8) +
  geom_line(data = vnpc_49_gg, aes(x = freq, y = fzy, colour = "vnpc(49)"),
            linewidth = .8) +
  geom_line(data = vnpc_99_gg, aes(x = freq, y = fzy, colour = "vnpc(99)"),
            linewidth = .8) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_colour_manual(values = c('average periodogram' = 'black', 'vnp' = "red3",
                                 'vnpc(49)' = "cyan3",  'vnpc(99)' = "olivedrab3")) +
  labs(x = NULL, y = NULL, title = expression(Im(S[ZY]))) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none')



pdf(file = "C:/Users/Yixuan/Documents/PhD/LISA/results/original/simulated_data/1yr/PSD_1yr_8193_vnpc_vs_vnp_plot.pdf",
    width = 18, height = 10)

(plot_xx | plot_xy | plot_xz) /
  (plot_yx | plot_yy | plot_yz) /
  (plot_zx | plot_zy | plot_zz)

dev.off()







plot_coh_xy <- ggplot(true_psd_gg) +
  geom_line(aes(freq, y = cohxy, colour = "true psd")) +
  geom_line(data = vnp_gg, aes(freq, y = cohxy, colour = "vnp")) +
  geom_line(data = vnpc_49_gg, aes(freq, y = cohxy, colour = "vnpc(49)")) +
  geom_line(data = vnpc_99_gg, aes(freq, y = cohxy, colour = "vnpc(99)")) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_colour_manual(values = c("true psd" = "black", "vnp" = "red3",
                                 "vnpc(49)" = "cyan3", "vnpc(99)" = "olivedrab3")) +
  labs(x = NULL, y = "Coherence", title = expression(Coherence[XY])) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none')


plot_coh_xz <- ggplot(true_psd_gg) +
  geom_line(aes(freq, y = cohxz, colour = "true psd")) +
  geom_line(data = vnp_gg, aes(freq, y = cohxz, colour = "vnp")) +
  geom_line(data = vnpc_49_gg, aes(freq, y = cohxz, colour = "vnpc(49)")) +
  geom_line(data = vnpc_99_gg, aes(freq, y = cohxz, colour = "vnpc(99)")) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_colour_manual(values = c("true psd" = "black", "vnp" = "red3",
                                 "vnpc(49)" = "cyan3", "vnpc(99)" = "olivedrab3")) +
  labs(x = NULL, y = "Coherence", title = expression(Coherence[XY])) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none')


plot_coh_yz <- ggplot(true_psd_gg) +
  geom_line(aes(freq, y = cohyz, colour = "true psd")) +
  geom_line(data = vnp_gg, aes(freq, y = cohyz, colour = "vnp")) +
  geom_line(data = vnpc_49_gg, aes(freq, y = cohyz, colour = "vnpc(49)")) +
  geom_line(data = vnpc_99_gg, aes(freq, y = cohyz, colour = "vnpc(99)")) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_colour_manual(values = c("true psd" = "black", "vnp" = "red3",
                                 "vnpc(49)" = "cyan3", "vnpc(99)" = "olivedrab3")) +
  labs(x = NULL, y = "Coherence", title = expression(Coherence[XY])) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20))


pdf(file = "C:/Users/Yixuan/Documents/PhD/LISA/results/original/simulated_data/1yr/Coherence_1yr_8193_vnpc_vs_vnp_plot.pdf",
    width = 18, height = 10)

plot_coh_xy | plot_coh_xz | plot_coh_yz

dev.off()




