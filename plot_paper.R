library(ggplot2)
library(patchwork)
library(scales)
library(reticulate)



# Simulation study --------------------------------------------------------

freq <- seq(0, pi, len = 513)

true_psd <- beyondWhittle::psd_varma(freq,
                                     ar = rbind(c(.5, 0, 0, 0), c(0, -.3, 0, -.5)),
                                     ma = matrix(nrow=2, ncol=0),
                                     Sigma = matrix(data=c(1, .9, .9, 1), nrow=2, ncol=2))


load("RES_MULTI__n1024_imodel1_repN100_algo2.RData")
mcmc_vnpc <- mcmc
load("RES_MULTI__n1024_imodel1_repN100_algo1.RData")
mcmc_vnp <- mcmc


true_psd_gg <- data.frame(
  f11 = Re(true_psd[1,1,-1]),
  f22 = Re(true_psd[2,2,-1]),
  f12 = Re(true_psd[1,2,-1]),
  f21 = Im(true_psd[2,1,-1]),
  freq = freq[-1]
)

vnp_gg <- data.frame(
  f11 = Re(mcmc_vnp$psd.median[1,1,-1]*.2/2/pi),
  f22 = Re(mcmc_vnp$psd.median[2,2,-1]*.2/2/pi),
  f12 = Re(mcmc_vnp$psd.median[1,2,-1]*.2/2/pi),
  f21 = Im(mcmc_vnp$psd.median[2,1,-1]*.2/2/pi),
  freq = freq[-1]
)

vnpc_gg <- data.frame(
  f11 = Re(mcmc_vnpc$psd.median[1,1,-1]*.2/2/pi),
  f22 = Re(mcmc_vnpc$psd.median[2,2,-1]*.2/2/pi),
  f12 = Re(mcmc_vnpc$psd.median[1,2,-1]*.2/2/pi),
  f21 = Im(mcmc_vnpc$psd.median[2,1,-1]*.2/2/pi),
  freq = freq[-1]
)



plot_11 <- ggplot(true_psd_gg) +
  geom_line(aes(x = freq, y = log(f11), 
                linetype = "True PSD", colour = "True PSD"), linewidth = .8) +
  geom_line(data = vnp_gg,
            aes(x = freq, y = log(f11), 
                linetype = "VNP", colour = "VNP"), linewidth = 1) +
  geom_line(data = vnpc_gg,
            aes(x = freq, y = log(f11), 
                linetype = "VNPC", colour = "VNPC"), linewidth = 1) +
  scale_linetype_manual(values = c("True PSD" = "solid",
                                   "VNP" = "dashed",
                                   "VNPC" = "dashed")) +
  scale_colour_manual(values = c("True PSD" = "black",
                                 "VNP" = "red3",
                                 "VNPC" = "blue3")) +
  labs(x = NULL, y = "PSD", title = expression(S[11])) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none')


plot_12 <- ggplot(true_psd_gg) +
  geom_line(aes(x = freq, y = f12, 
                linetype = "True PSD", colour = "True PSD"), linewidth = .8) +
  geom_line(data = vnp_gg,
            aes(x = freq, y = f12, 
                linetype = "VNP", colour = "VNP"), linewidth = .8) +
  geom_line(data = vnpc_gg,
            aes(x = freq, y = f12, 
                linetype = "VNPC", colour = "VNPC"), linewidth = .8) +
  scale_linetype_manual(values = c("True PSD" = "solid",
                                   "VNP" = "dashed",
                                   "VNPC" = "dashed")) +
  scale_colour_manual(values = c("True PSD" = "black",
                                 "VNP" = "red3",
                                 "VNPC" = "blue3")) +
  labs(x = NULL, y = NULL, title = expression(Re(S[12]))) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none')


plot_21 <- ggplot(true_psd_gg) +
  geom_line(aes(x = freq, y = f21, 
                linetype = "True PSD", colour = "True PSD"), linewidth = .8) +
  geom_line(data = vnp_gg,
            aes(x = freq, y = f21, 
                linetype = "VNP", colour = "VNP"), linewidth = .8) +
  geom_line(data = vnpc_gg,
            aes(x = freq, y = f21, 
                linetype = "VNPC", colour = "VNPC"), linewidth = .8) +
  scale_linetype_manual(values = c("True PSD" = "solid",
                                   "VNP" = "dashed",
                                   "VNPC" = "dashed")) +
  scale_colour_manual(values = c("True PSD" = "black",
                                 "VNP" = "red3",
                                 "VNPC" = "blue3")) +
  labs(x = "Frequency", y = "PSD", title = expression(Im(S[21]))) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none')


plot_22 <- ggplot(true_psd_gg) +
  geom_line(aes(x = freq, y = log(f22), 
                linetype = "True PSD", colour = "True PSD"), linewidth = .8) +
  geom_line(data = vnp_gg,
            aes(x = freq, y = log(f22), 
                linetype = "VNP", colour = "VNP"), linewidth = 1) +
  geom_line(data = vnpc_gg,
            aes(x = freq, y = log(f22), 
                linetype = "VNPC", colour = "VNPC"), linewidth = 1) +
  scale_linetype_manual(values = c("True PSD" = "solid",
                                   "VNP" = "dashed",
                                   "VNPC" = "dashed")) +
  scale_colour_manual(values = c("True PSD" = "black",
                                 "VNP" = "red3",
                                 "VNPC" = "blue3")) +
  labs(x = "Frequency", y = NULL, title = expression(S[22]), colour = "", linetype = "") +
  guides(colour = guide_legend(override.aes = list(fill = NA))) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22))


pdf(file = "C:/Users/Yixuan/Documents/PhD/LISA/simulation/simulation_paper/results/VAR(2)_1024_plot.pdf",
    width = 12, height = 9)

(plot_11 + plot_12) /
  (plot_21 + plot_22)
  

dev.off()






coh_mcmc <- function(Sii, Sjj, Sij, Sji) abs(Sij^2+Sji^2) / abs(Sii * Sjj) # coherence formula




# ET data -----------------------------------------------------------------

load("C:/Users/Yixuan/Documents/PhD/LISA/data/ET_data/caseA_true_psd.RData") # true psd data
load("C:/Users/Yixuan/Documents/PhD/LISA/data/ET_data/caseA_freq_psd.RData") # true frequency data

psd_true_gg <- data.frame(
  fx = psd_a[,1],
  fy = psd_a[,2],
  fz = psd_a[,3],
  freq = freq_psd_a
)


load("C:/Users/Yixuan/Documents/PhD/LISA/data/ET_data/caseA_mpg_orig.RData") # periodogram based on the full-length data
load("C:/Users/Yixuan/Documents/PhD/LISA/data/ET_data/caseA_freq_mpg.RData") # frequencies for the periodogram

mpg_orig_gg <- data.frame(
  fx = Re(mpg_a[1,1,-1]),
  fy = Re(mpg_a[2,2,-1]),
  fz = Re(mpg_a[3,3,-1]),
  fxy = Re(mpg_a[1,2,-1]),
  fxz = Re(mpg_a[1,3,-1]),
  fyz = Re(mpg_a[2,3,-1]),
  fyx = Im(mpg_a[2,1,-1]),
  fzx = Im(mpg_a[3,1,-1]),
  fzy = Im(mpg_a[3,2,-1]),
  coh_xy = coh_mcmc(Re(mpg_a[1,1,-1]),
                    Re(mpg_a[2,2,-1]),
                    mpg_a[1,2,-1],
                    mpg_a[2,1,-1]),
  coh_xz = coh_mcmc(Re(mpg_a[1,1,-1]),
                    Re(mpg_a[3,3,-1]),
                    mpg_a[1,3,-1],
                    mpg_a[3,1,-1]),
  coh_yz = coh_mcmc(Re(mpg_a[2,2,-1]),
                    Re(mpg_a[3,3,-1]),
                    mpg_a[2,3,-1],
                    mpg_a[3,2,-1]),
  freq = freq_a[-1]
)


load("C:/Users/Yixuan/Documents/PhD/LISA/data/ET_data/caseA_original.RData") # case A data

sig <- apply(caseA_data, 2, sd) # used for rescaled the psd estimates that are based on the standardized data


freq_mcmc <- seq(5, 128, len = 1969) # frequencies for the estimates

load("C:/Users/Yixuan/Documents/PhD/LISA/results/ET/mcmc_vnp_avg_ET_N20000.RData") # estimates from VNP

vnp_rescaled <- array(dim = c(3, 3, 1969)) # rescaled estimates from VNP
for(i in 1:1969) {                                                  
  vnp_rescaled[,,i] <- mcmc_vnp_avg$psd.median[,,i]*outer(sig, sig)
}

vnp_gg <- data.frame(
  fx = Re(vnp_rescaled[1,1,-1]),
  fy = Re(vnp_rescaled[2,2,-1]),
  fz = Re(vnp_rescaled[3,3,-1]),
  fxy = Re(vnp_rescaled[1,2,-1]),
  fxz = Re(vnp_rescaled[1,3,-1]),
  fyz = Re(vnp_rescaled[2,3,-1]),
  fyx = Im(vnp_rescaled[2,1,-1]),
  fzx = Im(vnp_rescaled[3,1,-1]),
  fzy = Im(vnp_rescaled[3,2,-1]),
  coh_xy = coh_mcmc(Re(vnp_rescaled[1,1,-1]),
                    Re(vnp_rescaled[2,2,-1]),
                    vnp_rescaled[1,2,-1],
                    vnp_rescaled[2,1,-1]),
  coh_xz = coh_mcmc(Re(vnp_rescaled[1,1,-1]),
                    Re(vnp_rescaled[3,3,-1]),
                    vnp_rescaled[1,3,-1],
                    vnp_rescaled[3,1,-1]),
  coh_yz = coh_mcmc(Re(vnp_rescaled[2,2,-1]),
                    Re(vnp_rescaled[3,3,-1]),
                    vnp_rescaled[2,3,-1],
                    vnp_rescaled[3,2,-1]),
  freq = freq_mcmc[-1]
)




load("C:/Users/Yixuan/Documents/PhD/LISA/results/ET/mcmc_vnpc_avg_n_seg_125_order_255_ET_N14000.RData") # estimates from VNPC

vnpc_rescaled <- array(dim = c(3, 3, 1969)) # rescaled estimates from VNPC
for(i in 1:1969) {
  vnpc_rescaled[,,i] <- mcmc_vnpc_avg$psd.median[,,i]*outer(sig, sig)
}

vnpc_gg <- data.frame(
  fx = Re(vnpc_rescaled[1,1,-1]),
  fy = Re(vnpc_rescaled[2,2,-1]),
  fz = Re(vnpc_rescaled[3,3,-1]),
  fxy = Re(vnpc_rescaled[1,2,-1]),
  fxz = Re(vnpc_rescaled[1,3,-1]),
  fyz = Re(vnpc_rescaled[2,3,-1]),
  fyx = Im(vnpc_rescaled[2,1,-1]),
  fzx = Im(vnpc_rescaled[3,1,-1]),
  fzy = Im(vnpc_rescaled[3,2,-1]),
  coh_xy = coh_mcmc(Re(vnpc_rescaled[1,1,-1]),
                    Re(vnpc_rescaled[2,2,-1]),
                    vnpc_rescaled[1,2,-1],
                    vnpc_rescaled[2,1,-1]),
  coh_xz = coh_mcmc(Re(vnpc_rescaled[1,1,-1]),
                    Re(vnpc_rescaled[3,3,-1]),
                    vnpc_rescaled[1,3,-1],
                    vnpc_rescaled[3,1,-1]),
  coh_yz = coh_mcmc(Re(vnpc_rescaled[2,2,-1]),
                    Re(vnpc_rescaled[3,3,-1]),
                    vnpc_rescaled[2,3,-1],
                    vnpc_rescaled[3,2,-1]),
  freq = freq_mcmc[-1]
)



#----------------------------------psd-------------------------------------#

plot_xx <- ggplot(mpg_orig_gg) +
  geom_line(aes(x = freq, y = fx, colour = "Periodogram")) +
  geom_line(data = psd_true_gg,
            aes(x = freq, y = fx, colour = "True PSD")) +
  geom_line(data = vnp_gg,
            aes(x = freq, y = fx, colour = "VNP"), linewidth = .8) +
  geom_line(data = vnpc_gg,
            aes(x = freq, y = fx, colour = "VNPC(255)"), linewidth = .8) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)),
                limits = c(5e-52, NA),              # lower bound slightly lower
                breaks = c(10^-47, 10^-49, 10^-51)) +
  scale_x_continuous(limits = c(5, 128), breaks = seq(25, 125, by = 25)) +
  scale_colour_manual(values = c("Periodogram" = "lightgrey",
                                 "True PSD" = "black",
                                 "VNP" = "red3",
                                 "VNPC(255)" = "cyan3")) +
  labs(x = NULL, y = NULL, title = expression(S[XX])) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none',
        panel.grid = element_blank())


plot_yy <- ggplot(mpg_orig_gg) +
  geom_line(aes(x = freq, y = fy, colour = "Periodogram")) +
  geom_line(data = psd_true_gg,
            aes(x = freq, y = fy, colour = "True PSD")) +
  geom_line(data = vnp_gg,
            aes(x = freq, y = fy, colour = "VNP"), linewidth = .8) +
  geom_line(data = vnpc_gg,
            aes(x = freq, y = fy, colour = "VNPC(255)"), linewidth = .8) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)),
                limits = c(5e-52, NA),              # lower bound slightly lower
                breaks = c(10^-47, 10^-49, 10^-51)) +
  scale_x_continuous(limits = c(5, 128), breaks = seq(25, 125, by = 25)) +
  scale_colour_manual(values = c("Periodogram" = "lightgrey",
                                 "True PSD" = "black",
                                 "VNP" = "red3",
                                 "VNPC(255)" = "cyan3")) +
  labs(x = NULL, y = NULL, title = expression(S[YY])) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none',
        panel.grid = element_blank())


plot_zz <- ggplot(mpg_orig_gg) +
  geom_line(aes(x = freq, y = fz, colour = "Periodogram")) +
  geom_line(data = psd_true_gg,
            aes(x = freq, y = fz, colour = "True PSD")) +
  geom_line(data = vnp_gg,
            aes(x = freq, y = fz, colour = "VNP"), linewidth = .8) +
  geom_line(data = vnpc_gg,
            aes(x = freq, y = fz, colour = "VNPC(255)"), linewidth = .8) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)),
                limits = c(5e-52, NA),              # lower bound slightly lower
                breaks = c(10^-47, 10^-49, 10^-51)) +
  scale_x_continuous(limits = c(5, 128), breaks = seq(25, 125, by = 25)) +
  scale_colour_manual(values = c("Periodogram" = "lightgrey",
                                 "True PSD" = "black",
                                 "VNP" = "red3",
                                 "VNPC(255)" = "cyan3")) +
  labs(x = NULL, y = NULL, title = expression(S[ZZ]), colour = "") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22),
        plot.title = element_text(size = 20),
        panel.grid = element_blank())



plot_xy <- ggplot(mpg_orig_gg) +
  geom_line(aes(x = freq, y = fxy, colour = "Periodogram")) +
  geom_line(data = vnp_gg,
            aes(x = freq, y = fxy, colour = "VNP"), linewidth = .8) +
  geom_line(data = vnpc_gg,
            aes(x = freq, y = fxy, colour = "VNPC(255)"), linewidth = .8) +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1e-50),
                     breaks = c(-10^-49, 0, 10^-49),
                     labels = function(x) {
                       sapply(x, function(v) {
                         if (v == 0) return("0")
                         sign <- ifelse(v < 0, "-", "")
                         paste0(sign, "10^", round(log10(abs(v))))
                       }) |> parse(text = _)
                     }) +
  scale_x_continuous(limits = c(5, 128), breaks = seq(25, 125, by = 25)) +
  scale_colour_manual(values = c("Periodogram" = "lightgrey",
                                 "VNP" = "red3",
                                 "VNPC(255)" = "cyan3")) +
  labs(x = NULL, y = NULL, title = expression(Re(S[xy]))) +
  theme_grey() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none',
        panel.grid = element_blank())



plot_xz <- ggplot(mpg_orig_gg) +
  geom_line(aes(x = freq, y = fxz, colour = "Periodogram")) +
  geom_line(data = vnp_gg,
            aes(x = freq, y = fxz, colour = "VNP"), linewidth = .8) +
  geom_line(data = vnpc_gg,
            aes(x = freq, y = fxz, colour = "VNPC(255)"), linewidth = .8) +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1e-50),
                     breaks = c(-10^-49, 0, 10^-49),
                     labels = function(x) {
                       sapply(x, function(v) {
                         if (v == 0) return("0")
                         sign <- ifelse(v < 0, "-", "")
                         paste0(sign, "10^", round(log10(abs(v))))
                       }) |> parse(text = _)
                     }) +
  scale_x_continuous(limits = c(5, 128), breaks = seq(25, 125, by = 25)) +
  scale_colour_manual(values = c("Periodogram" = "lightgrey",
                                 "VNP" = "red3",
                                 "VNPC(255)" = "cyan3")) +
  labs(x = NULL, y = NULL, title = expression(Re(S[xz]))) +
  theme_grey() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none',
        panel.grid = element_blank())



plot_yz <- ggplot(mpg_orig_gg) +
  geom_line(aes(x = freq, y = fyz, colour = "Periodogram")) +
  geom_line(data = vnp_gg,
            aes(x = freq, y = fyz, colour = "VNP"), linewidth = .8) +
  geom_line(data = vnpc_gg,
            aes(x = freq, y = fyz, colour = "VNPC(255)"), linewidth = .8) +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1e-50),
                     breaks = c(-10^-49, 0, 10^-49),
                     labels = function(x) {
                       sapply(x, function(v) {
                         if (v == 0) return("0")
                         sign <- ifelse(v < 0, "-", "")
                         paste0(sign, "10^", round(log10(abs(v))))
                       }) |> parse(text = _)
                     }) +
  scale_x_continuous(limits = c(5, 128), breaks = seq(25, 125, by = 25)) +
  scale_colour_manual(values = c("Periodogram" = "lightgrey",
                                 "VNP" = "red3",
                                 "VNPC(255)" = "cyan3")) +
  labs(x = NULL, y = "PSD", title = expression(Re(S[yz]))) +
  theme_grey() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none',
        panel.grid = element_blank())



plot_yx <- ggplot(mpg_orig_gg) +
  geom_line(aes(x = freq, y = fyx, colour = "Periodogram")) +
  geom_line(data = vnp_gg,
            aes(x = freq, y = fyx, colour = "VNP"), linewidth = .8) +
  geom_line(data = vnpc_gg,
            aes(x = freq, y = fyx, colour = "VNPC(255)"), linewidth = .8) +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1e-50),
                     breaks = c(-10^-49, 0, 10^-49),
                     labels = function(x) {
                       sapply(x, function(v) {
                         if (v == 0) return("0")
                         sign <- ifelse(v < 0, "-", "")
                         paste0(sign, "10^", round(log10(abs(v))))
                       }) |> parse(text = _)
                     }) +
  scale_x_continuous(limits = c(5, 128), breaks = seq(25, 125, by = 25)) +
  scale_colour_manual(values = c("Periodogram" = "lightgrey",
                                 "VNP" = "red3",
                                 "VNPC(255)" = "cyan3")) +
  labs(x = NULL, y = "Strain PSD [1/Hz]", title = expression(Im(S[yx]))) +
  theme_grey() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none',
        panel.grid = element_blank())



plot_zx <- ggplot(mpg_orig_gg) +
  geom_line(aes(x = freq, y = fzx, colour = "Periodogram")) +
  geom_line(data = vnp_gg,
            aes(x = freq, y = fzx, colour = "VNP"), linewidth = .8) +
  geom_line(data = vnpc_gg,
            aes(x = freq, y = fzx, colour = "VNPC(255)"), linewidth = .8) +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1e-50),
                     breaks = c(-10^-49, 0, 10^-49),
                     labels = function(x) {
                       sapply(x, function(v) {
                         if (v == 0) return("0")
                         sign <- ifelse(v < 0, "-", "")
                         paste0(sign, "10^", round(log10(abs(v))))
                       }) |> parse(text = _)
                     }) +
  scale_x_continuous(limits = c(5, 128), breaks = seq(25, 125, by = 25)) +
  scale_colour_manual(values = c("Periodogram" = "lightgrey",
                                 "VNP" = "red3",
                                 "VNPC(255)" = "cyan3")) +
  labs(x = NULL, y = NULL, title = expression(Im(S[zx]))) +
  theme_grey() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none',
        panel.grid = element_blank())



plot_zy <- ggplot(mpg_orig_gg) +
  geom_line(aes(x = freq, y = fzy, colour = "Periodogram")) +
  geom_line(data = vnp_gg,
            aes(x = freq, y = fzy, colour = "VNP"), linewidth = .8) +
  geom_line(data = vnpc_gg,
            aes(x = freq, y = fzy, colour = "VNPC(255)"), linewidth = .8) +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1e-50),
                     breaks = c(-10^-49, 0, 10^-49),
                     labels = function(x) {
                       sapply(x, function(v) {
                         if (v == 0) return("0")
                         sign <- ifelse(v < 0, "-", "")
                         paste0(sign, "10^", round(log10(abs(v))))
                       }) |> parse(text = _)
                     }) +
  scale_x_continuous(limits = c(5, 128), breaks = seq(25, 125, by = 25)) +
  scale_colour_manual(values = c("Periodogram" = "lightgrey",
                                 "VNP" = "red3",
                                 "VNPC(255)" = "cyan3")) +
  labs(x = "Frequency [Hz]", y = NULL, title = expression(Im(S[zy]))) +
  theme_grey() +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none',
        panel.grid = element_blank())



# pdf(file = "C:/Users/Yixuan/Documents/PhD/LISA/results/ET/psd_ET_order_255_plot.pdf", 
#     width = 18, height = 10, compress = TRUE)
# 
# (plot_xx + plot_xy + plot_xz) /
#   (plot_yx + plot_yy + plot_yz) /
#   (plot_zx + plot_zy + plot_zz)
# 
# 
# dev.off()


png(file = "C:/Users/Yixuan/Documents/PhD/LISA/results/ET/psd_ET_order_255_plot.png",
    width = 5400, height = 3000, units = "px", res = 300)

(plot_xx + plot_xy + plot_xz) /
  (plot_yx + plot_yy + plot_yz) /
  (plot_zx + plot_zy + plot_zz)

dev.off()






#----------------------------coherence----------------------------------#

ET_1 <- read.table("ET.txt")
Peak10Hz <- read.table("Peak10Hz_new.txt")
Peak50Hz <- read.table("Peak50Hz_new.txt")
Peak90Hz <- read.table("Peak90Hz_new.txt")

ASD_X <- ET_1[,2] + Peak10Hz[,2] + Peak50Hz[,2]
ASD_Y <- ET_1[,2] + Peak10Hz[,2] + Peak90Hz[,2]
ASD_Z <- ET_1[,2] + Peak50Hz[,2] + Peak90Hz[,2]

coh_xy <- Peak10Hz[,2]^2 / (ASD_X * ASD_Y)
coh_xz <- Peak50Hz[,2]^2 / (ASD_X * ASD_Z)
coh_yz <- Peak90Hz[,2]^2 / (ASD_Y * ASD_Z)
freq_coh <- ET_1[,1]

coh_gg <- data.frame(
  coh_xy = coh_xy[526:1528],
  coh_xz = coh_xz[526:1528],
  coh_yz = coh_yz[526:1528],
  freq = freq_coh[526:1528]
)


plot_coh <- ggplot(coh_gg) +
  geom_line(aes(freq, y = coh_xy, colour = "True Coherence"), linewidth = 1) +
  geom_line(aes(freq, y = coh_xz, colour = "True Coherence"), linewidth = 1) +
  geom_line(aes(freq, y = coh_yz, colour = "True Coherence"), linewidth = 1) +
  geom_line(data = vnp_gg, aes(freq, y = coh_xy, colour = "VNP"), linewidth = .8) +
  geom_line(data = vnp_gg, aes(freq, y = coh_xz, colour = "VNP"), linewidth = .8) +
  geom_line(data = vnp_gg, aes(freq, y = coh_yz, colour = "VNP"), linewidth = .8) +
  geom_line(data = vnpc_gg, aes(freq, y = coh_xy, colour = "VNPC(255)"), linewidth = .8) +
  geom_line(data = vnpc_gg, aes(freq, y = coh_xz, colour = "VNPC(255)"), linewidth = .8) +
  geom_line(data = vnpc_gg, aes(freq, y = coh_yz, colour = "VNPC(255)"), linewidth = .8) +
  scale_colour_manual(values = c("True Coherence" = "black", 
                                 "VNP" = "red3", "VNPC(255)" = "cyan3")) +
  labs(x = NULL, y = "Coherence", title = expression(Coherence[XY]), colour = "") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22),
        plot.title = element_text(size = 20))


pdf(file = "coherence_ET_order_255_plot.pdf", 
    width = 10, height = 6)

plot_coh

dev.off()








# LISA noise --------------------------------------------------------------

np <- import("numpy")
data_npy <- np$load("data.npy") # original data
freq_true <- np$load("freq_true.npy") # true frequencies
true_matrix <- np$load("true_matrix.npy") # true psd

# Trim the plot as plotting the full true-psd can kill R!
# pick points evenly spaced on log scale of index
m <- 10000  # target number of points to plot
idx <- unique(round(exp(seq(log(1), log(length(freq_true)), length.out = m)))) # create indices that are evenly spaced on a log scale

true_matrix_sub <- true_matrix[idx,,]
freq_true_sub <- freq_true[idx]


true_psd_gg <- data.frame(
  fx = true_matrix_sub[1634:6321,1,1],
  fy = true_matrix_sub[1634:6321,2,2],
  fz = true_matrix_sub[1634:6321,3,3],
  fxy = true_matrix_sub[1634:6321,1,2],
  fxz = true_matrix_sub[1634:6321,1,3],
  fyz = true_matrix_sub[1634:6321,2,3],
  fyx = true_matrix_sub[1634:6321,2,1]/sd(true_matrix_sub[1634:6321,2,1]),
  fzx = true_matrix_sub[1634:6321,3,1]/sd(true_matrix_sub[1634:6321,3,1]),
  fzy = true_matrix_sub[1634:6321,3,2]/sd(true_matrix_sub[1634:6321,3,2]),
  cohxy = coh_mcmc(true_matrix_sub[1634:6321,1,1],
                   true_matrix_sub[1634:6321,2,2],
                   true_matrix_sub[1634:6321,1,2],
                   true_matrix_sub[1634:6321,2,1]),
  cohxz = coh_mcmc(true_matrix_sub[1634:6321,1,1],
                   true_matrix_sub[1634:6321,3,3],
                   true_matrix_sub[1634:6321,1,3],
                   true_matrix_sub[1634:6321,3,1]),
  cohyz = coh_mcmc(true_matrix_sub[1634:6321,2,2],
                   true_matrix_sub[1634:6321,3,3],
                   true_matrix_sub[1634:6321,2,3],
                   true_matrix_sub[1634:6321,3,2]),
  freq = freq_true_sub[1634:6321]
)


freq_vnp <- seq(0, .1, len = 8193) # frequencies for the estimates

sig <- apply(data_npy, 2, sd) 
hann_window <- hanning(16384)

load("mcmc_vnp_avg_lisa_simu_1yr_16384_N80000_70000.RData")

vnp_rescaled <- array(0, dim = c(3, 3, 8193))
for (i in 1:8193) {
  vnp_rescaled[,,i] <- mcmc_vnp_avg$psd.median[,,i]/mean(hann_window^2)*2*outer(sig, sig)
}



vnp_gg <- data.frame(
  fx = Re(vnp_rescaled[1,1,9:8193]),
  fy = Re(vnp_rescaled[2,2,9:8193]),
  fz = Re(vnp_rescaled[3,3,9:8193]),
  fxy = Re(vnp_rescaled[1,2,9:8193]),
  fxz = Re(vnp_rescaled[1,3,9:8193]),
  fyz = Re(vnp_rescaled[2,3,9:8193]),
  fyx = Im(vnp_rescaled[2,1,9:8193]),
  fzx = Im(vnp_rescaled[3,1,9:8193]),
  fzy = Im(vnp_rescaled[3,2,9:8193]),
  cohxy = coh_mcmc(vnp_rescaled[1,1,9:8193],
                   vnp_rescaled[2,2,9:8193],
                   vnp_rescaled[1,2,9:8193],
                   vnp_rescaled[2,1,9:8193]),
  cohxz = coh_mcmc(vnp_rescaled[1,1,9:8193],
                   vnp_rescaled[3,3,9:8193],
                   vnp_rescaled[1,3,9:8193],
                   vnp_rescaled[3,1,9:8193]),
  cohyz = coh_mcmc(vnp_rescaled[2,2,9:8193],
                   vnp_rescaled[3,3,9:8193],
                   vnp_rescaled[2,3,9:8193],
                   vnp_rescaled[3,2,9:8193]),
  freq = freq_vnp[9:8193]
)



load("mcmc_vnpc_avg_order_34_lisa_simu_1yr_16384_N50000_45000.RData")

vnpc_rescaled <- array(0, dim = c(3, 3, 8193))
for (i in 1:8193) {
  vnpc_rescaled[,,i] <- mcmc_vnpc_avg$psd.median[,,i]/mean(hann_window^2)*2*outer(sig, sig)
}

vnpc_34_gg <- data.frame(
  fx = Re(vnpc_rescaled[1,1,9:8193]),
  fy = Re(vnpc_rescaled[2,2,9:8193]),
  fz = Re(vnpc_rescaled[3,3,9:8193]),
  fxy = Re(vnpc_rescaled[1,2,9:8193]),
  fxz = Re(vnpc_rescaled[1,3,9:8193]),
  fyz = Re(vnpc_rescaled[2,3,9:8193]),
  fyx = Im(vnpc_rescaled[2,1,9:8193]),
  fzx = Im(vnpc_rescaled[3,1,9:8193]),
  fzy = Im(vnpc_rescaled[3,2,9:8193]),
  cohxy = coh_mcmc(Re(vnpc_rescaled[1,1,9:8193]),
                   Re(vnpc_rescaled[2,2,9:8193]),
                   vnpc_rescaled[1,2,9:8193],
                   vnpc_rescaled[2,1,9:8193]),
  cohxz = coh_mcmc(Re(vnpc_rescaled[1,1,9:8193]),
                   Re(vnpc_rescaled[3,3,9:8193]),
                   vnpc_rescaled[1,3,9:8193],
                   vnpc_rescaled[3,1,9:8193]),
  cohyz = coh_mcmc(Re(vnpc_rescaled[2,2,9:8193]),
                   Re(vnpc_rescaled[3,3,9:8193]),
                   vnpc_rescaled[2,3,9:8193],
                   vnpc_rescaled[3,2,9:8193]),
  freq = freq_vnp[9:8193]
)



load("mpgm_avg_16384.RData") # average periodogram

avg_mpg_rescaled <- array(0, dim = c(3, 3, 8193))
for (i in 1:8193) {
  avg_mpg_rescaled[,,i] <- mpgm_avg[,,i]/mean(hann_window^2)*2*outer(sig, sig)
}

avg_mpg_gg <- data.frame(
  fx = Re(avg_mpg_rescaled[1,1,9:8193]),
  fy = Re(avg_mpg_rescaled[2,2,9:8193]),
  fz = Re(avg_mpg_rescaled[3,3,9:8193]),
  fxy = Re(avg_mpg_rescaled[1,2,9:8193]),
  fxz = Re(avg_mpg_rescaled[1,3,9:8193]),
  fyz = Re(avg_mpg_rescaled[2,3,9:8193]),
  fyx = Im(avg_mpg_rescaled[2,1,9:8193]),
  fzx = Im(avg_mpg_rescaled[3,1,9:8193]),
  fzy = Im(avg_mpg_rescaled[3,2,9:8193]),
  cohxy = coh_mcmc(Re(avg_mpg_rescaled[1,1,9:8193]),
                   Re(avg_mpg_rescaled[2,2,9:8193]),
                   Re(avg_mpg_rescaled[1,2,9:8193]),
                   Re(avg_mpg_rescaled[2,1,9:8193])),
  cohxz = coh_mcmc(Re(avg_mpg_rescaled[1,1,9:8193]),
                   Re(avg_mpg_rescaled[3,3,9:8193]),
                   Re(avg_mpg_rescaled[1,3,9:8193]),
                   Re(avg_mpg_rescaled[3,1,9:8193])),
  cohyz = coh_mcmc(Re(avg_mpg_rescaled[2,2,9:8193]),
                   Re(avg_mpg_rescaled[3,3,9:8193]),
                   Re(avg_mpg_rescaled[2,3,9:8193]),
                   Re(avg_mpg_rescaled[3,2,9:8193])),
  freq = freq_vnp[9:8193]
)



#---------------------elbow criterion----------------------------#

load("var_order_16384.RData") # log-likelihood for the VAR fitting with different orders

plot_order <- ggplot(data.frame(llike_var), aes(3:102, -llike_var)) +
  geom_line(linewidth = 1) +
  geom_point(colour = 'red3', size = 1.5) +
  geom_vline(xintercept = seq(3,102,by=2), color = 'blue', alpha = .3) +
  scale_x_continuous(breaks = seq(3, 102, by = 6)) +
  scale_y_continuous(labels = NULL) +
  xlab("lag") +
  ylab("-log-likelihood") +
  theme(axis.title = element_text(size = 20)) +
  theme(axis.text = element_text(size = 15))


pdf(file = "fit_order_1yr_16384_plot.pdf",
    width = 14, height = 6)

plot_order

dev.off()



#---------------------------------psd---------------------------------------#

plot_xx <- ggplot(avg_mpg_gg) +
  geom_line(aes(x = freq, y = fx, colour = "Average Periodogram"), linewidth = .8) +
  geom_line(data = true_psd_gg, aes(x = freq, y = fx, colour = "True PSD"),
            linewidth = .8) +
  geom_line(data = vnp_gg, aes(x = freq, y = fx, colour = "VNP"),
            linewidth = .8) +
  # geom_line(data = vnpc_34_gg, aes(x = freq, y = fx, colour = "vnpc(34)"),
  #           linewidth = .8) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_colour_manual(values = c('True PSD' = 'black', #'S_true' = 'red3', 'Welch psd' = 'cyan3', 
                                 'VNP' = "red3",
                                 # 'vnpc(34)' = "cyan3",  
                                 'Average Periodogram' = "slategrey")) +
  labs(x = NULL, y = "PSD", title = expression(S[XX])) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none')


plot_yy <- ggplot(avg_mpg_gg) +
  geom_line(aes(x = freq, y = fy, colour = "Average Periodogram"), linewidth = .8) +
  geom_line(data = true_psd_gg, aes(x = freq, y = fy, colour = "True PSD"),
            linewidth = .8) +
  geom_line(data = vnp_gg, aes(x = freq, y = fy, colour = "VNP"),
            linewidth = .8) +
  # geom_line(data = vnpc_34_gg, aes(x = freq, y = fy, colour = "vnpc(34)"),
  #           linewidth = .8) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_colour_manual(values = c('True PSD' = 'black', #'S_true' = 'red3', 'Welch psd' = 'cyan3', 
                                 'VNP' = "red3", 
                                 # 'vnpc(34)' = "cyan3",  
                                 'Average Periodogram' = "slategrey")) +
  labs(x = NULL, y = NULL, title = expression(S[YY])) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none')


plot_zz <- ggplot(avg_mpg_gg) +
  geom_line(aes(x = freq, y = fz, colour = "Average Periodogram"), linewidth = .8) +
  geom_line(data = true_psd_gg, aes(x = freq, y = fz, colour = "True PSD"),
            linewidth = .8) +
  geom_line(data = vnp_gg, aes(x = freq, y = fz, colour = "VNP"),
            linewidth = .8) +
  # geom_line(data = vnpc_34_gg, aes(x = freq, y = fz, colour = "vnpc(34)"),
  #           linewidth = .8) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_colour_manual(values = c('True PSD' = 'black', # 'S_true' = 'red3', 'Welch psd' = 'cyan3',
                                 'VNP' = "red3", 
                                 # 'vnpc(34)' = "cyan3",  
                                 'Average Periodogram' = "slategrey")) +
  labs(x = NULL, y = NULL, title = expression(S[ZZ]), colour = "Method") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none')


plot_xy <- ggplot(avg_mpg_gg) +
  geom_line(aes(x = freq, y = fxy, colour = "Average Periodogram"), linewidth = .8) +
  geom_line(data = true_psd_gg, aes(x = freq, y = fxy, colour = "True PSD"),
            linewidth = .8) +
  geom_line(data = vnp_gg, aes(x = freq, y = fxy, colour = "VNP"),
            linewidth = 0.8) +
  # geom_line(data = vnpc_34_gg, aes(x = freq, y = fxy, colour = "vnpc(34)"),
  #           linewidth = .8) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(
    labels = label_scientific(digits = 1)
  ) +
  scale_colour_manual(values = c('True PSD' = 'black', #'S_true' = 'red3', 'Welch psd' = 'cyan3',
                                 'VNP' = "red3", 
                                 # 'vnpc(34)' = "cyan3",  
                                 'Average Periodogram' = "slategrey")) +
  labs(x = NULL, y = NULL, title = expression(Re(S[XY]))) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none')


plot_xz <- ggplot(avg_mpg_gg) +
  geom_line(aes(x = freq, y = fxz, colour = "Average Periodogram"), linewidth = .8) +
  geom_line(data = true_psd_gg, aes(x = freq, y = fxz, colour = "True PSD"),
            linewidth = .8) +
  geom_line(data = vnp_gg, aes(x = freq, y = fxz, colour = "VNP"),
            linewidth = 0.8) +
  # geom_line(data = vnpc_34_gg, aes(x = freq, y = fxz, colour = "vnpc(34)"),
  #           linewidth = .8) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(
    labels = label_scientific(digits = 1)
  ) +
  scale_colour_manual(values = c('True PSD' = 'black', #'S_true' = 'red3', 'Welch psd' = 'cyan3', 
                                 'VNP' = "red3", 
                                 # 'vnpc(34)' = "cyan3",
                                 'Average Periodogram' = "slategrey")) +
  labs(x = NULL, y = NULL, title = expression(Re(S[XZ]))) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none')


plot_yz <- ggplot(avg_mpg_gg) +
  geom_line(aes(x = freq, y = fyz, colour = "Average Periodogram"), linewidth = .8) +
  geom_line(data = true_psd_gg, aes(x = freq, y = fyz, colour = "True PSD"),
            linewidth = .8) +
  geom_line(data = vnp_gg, aes(x = freq, y = fyz, colour = "VNP"),
            linewidth = 0.8) +
  # geom_line(data = vnpc_34_gg, aes(x = freq, y = fyz, colour = "vnpc(34)"),
  #           linewidth = .8) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(
    labels = label_scientific(digits = 1)
  ) +
  scale_colour_manual(values = c('True PSD' = 'black', #'S_true' = 'red3', 'Welch psd' = 'cyan3', 
                                 'VNP' = "red3", 
                                 # 'vnpc(34)' = "cyan3",
                                 'Average Periodogram' = "slategrey")) +
  labs(x = NULL, y = NULL, title = expression(Re(S[YZ])), colour = "Method") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22),
        plot.title = element_text(size = 20))


plot_yx <- ggplot(avg_mpg_gg) +
  geom_line(aes(x = freq, y = fyx, colour = "Average Periodogram")) +
  geom_line(data = vnp_gg, aes(x = freq, y = fyx, colour = "VNP"),
            linewidth = 0.8) +
  # geom_line(data = vnpc_34_gg, aes(x = freq, y = fyx, colour = "vnpc(34)"),
  #           linewidth = .8) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(
    labels = label_scientific(digits = 1)
  ) +
  scale_colour_manual(values = c('VNP' = "red3", #'Welch psd' = 'cyan3',
                                 # 'vnpc(34)' = "cyan3", 
                                 'Average Periodogram' = "slategrey")) +
  labs(x = NULL, y = "PSD", title = expression(Im(S[YX]))) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none')


plot_zx <- ggplot(avg_mpg_gg) +
  geom_line(aes(x = freq, y = fzx, colour = "Average Periodogram")) +
  geom_line(data = vnp_gg, aes(x = freq, y = fzx, colour = "VNP"),
            linewidth = 0.8) +
  # geom_line(data = vnpc_34_gg, aes(x = freq, y = fzx, colour = "vnpc(34)"),
  #           linewidth = .8) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(
    labels = label_scientific(digits = 1)
  ) +
  scale_colour_manual(values = c('VNP' = "red3", 
                                 # 'vnpc(34)' = "cyan3", 
                                 'Average Periodogram' = "slategrey")) +
  labs(x = NULL, y = "PSD", title = expression(Im(S[ZX]))) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none')


plot_zy <- ggplot(avg_mpg_gg) +
  geom_line(aes(x = freq, y = fzy, colour = "Average Periodogram")) +
  geom_line(data = vnp_gg, aes(x = freq, y = fzy, colour = "VNP"),
            linewidth = 0.8) +
  # geom_line(data = vnpc_34_gg, aes(x = freq, y = fzy, colour = "vnpc(34)"),
  #           linewidth = .8) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(
    labels = label_scientific(digits = 1)
  ) +
  scale_colour_manual(values = c('VNP' = "red3", 
                                 # 'vnpc(34)' = "cyan3", 
                                 'Average Periodogram' = "slategrey")) +
  labs(x = NULL, y = NULL, title = expression(Im(S[ZY]))) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none')




# pdf(file = "PSD_1yr_16384_order_34_plot.pdf",
#     width = 18, height = 10)
# 
# (plot_xx | plot_xy | plot_xz) /
#   (plot_yx | plot_yy | plot_yz) /
#   (plot_zx | plot_zy | plot_zz)
# 
# dev.off()

pdf(file = "PSD_1yr_16384_vnp_plot.pdf",
    width = 18, height = 10)

(plot_xx | plot_xy | plot_xz) /
  (plot_yx | plot_yy | plot_yz) /
  (plot_zx | plot_zy | plot_zz)

dev.off()




#---------------------------------coherence---------------------------------------#

plot_coh_xy <- ggplot(true_psd_gg) +
  geom_line(aes(freq, y = cohxy, colour = "True PSD"), linewidth = 1) +
  geom_line(data = vnp_gg, aes(freq, y = cohxy, colour = "VNP"), linewidth = .8) +
  # geom_line(data = vnpc_34_gg, aes(freq, y = cohxy, colour = "vnpc(34)"), linewidth = .8) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_colour_manual(values = c("True PSD" = "black",
                                 # "vnpc(34)" = "cyan3",
                                 "VNP" = "red3")) +
  labs(x = NULL, y = "Coherence", title = expression(Coherence[XY])) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none')


plot_coh_xz <- ggplot(true_psd_gg) +
  geom_line(aes(freq, y = cohxz, colour = "True PSD"), linewidth = 1) +
  geom_line(data = vnp_gg, aes(freq, y = cohxz, colour = "VNP"), linewidth = .8) +
  # geom_line(data = vnpc_34_gg, aes(freq, y = cohxz, colour = "vnpc(34)"), linewidth = .8) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_colour_manual(values = c("True PSD" = "black", 
                                 # "vnpc(34)" = "cyan3",
                                 "VNP" = "red3")) +
  labs(x = NULL, y = NULL, title = expression(Coherence[XZ])) +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        legend.position = 'none')


plot_coh_yz <- ggplot(true_psd_gg) +
  geom_line(aes(freq, y = cohyz, colour = "True PSD"), linewidth = 1) +
  geom_line(data = vnp_gg, aes(freq, y = cohyz, colour = "VNP"), linewidth = .8) +
  # geom_line(data = vnpc_34_gg, aes(freq, y = cohyz, colour = "vnpc(34)"), linewidth = .8) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_colour_manual(values = c("True PSD" = "black",
                                 # "vnpc(34)" = "cyan3",
                                 "VNP" = "red3")) +
  labs(x = NULL, y = NULL, title = expression(Coherence[YZ]), colour = "Method") +
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 22),
        plot.title = element_text(size = 20))


# pdf(file = "Coherence_1yr_16384_order_34_plot.pdf",
#     width = 18, height = 10)
# 
# plot_coh_xy | plot_coh_xz | plot_coh_yz
# 
# dev.off()

pdf(file = "Coherence_1yr_16384_vnp_plot.pdf",
    width = 18, height = 10)

plot_coh_xy | plot_coh_xz | plot_coh_yz

dev.off()






