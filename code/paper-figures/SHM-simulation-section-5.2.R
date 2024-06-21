# Packages needed: --------------------------------------------------------
library(ggplot2)    # CRAN v3.3.5
library(data.table) # CRAN v1.14.2
library(tikzDevice) # CRAN v0.12.3.1
library(ggpubr)     # CRAN v0.4.0 # CRAN v0.4.0

options(tikzLatexPackages = c(getOption( "tikzLatexPackages" ),"\\usepackage{amsfonts}"))
options(tikzLatexPackages = c(getOption( "tikzLatexPackages" ),"\\usepackage{amsmath}"))
options(tikzLatexPackages = c(getOption( "tikzLatexPackages" ),"\\usepackage{amssymb}"))

# Some initial parameters: ------------------------------------------------
grid_range <- c(0, 2 * pi)
n_grid <- 100
grid_points <- seq(from = grid_range[1],
                   to = grid_range[2],
                   length.out = n_grid)

doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937

# Read in results of Simulation: ------------------------------------------
simulation_result_list_baseline <- readRDS(here::here("outputs",
                                             "SHM",
                                             "simulation-results",
                                             "simulation-section-4.2-updated.rds"))
simulation_result_list_extension_01 <- readRDS(here::here("outputs",
                                                      "SHM",
                                                      "simulation-results",
                                                      "simulation-section-4.2-extension-ic.rds"))

simulation_result_list_extension_02 <- readRDS(here::here("outputs",
                                                          "SHM",
                                                          "simulation-results",
                                                          "simulation-section-4.2-extension-intensity.rds"))

simulation_result_list_extension_03 <- readRDS(here::here("outputs",
                                                          "SHM",
                                                          "simulation-results",
                                                          "simulation-section-4.2-extension-forcing-variance.rds"))



# Wrangling of results: ---------------------------------------------------
results_array_baseline <- array(unlist(simulation_result_list_baseline$results_list), dim = c(n_grid, 4, 500))
results_array_extension_01 <- simulation_result_list_extension_01$results_array
results_array_extension_02 <- simulation_result_list_extension_02$results_array
results_array_extension_03 <- simulation_result_list_extension_03$results_array

settings_extension_01 <- simulation_result_list_extension_01$settings_dt
settings_extension_02 <- simulation_result_list_extension_02$settings_dt
settings_extension_03 <- simulation_result_list_extension_03$settings_dt

error_array_baseline <- results_array_baseline + 1
error_array_extension_01 <- results_array_extension_01 + 1
error_array_extension_02 <- results_array_extension_02 + 1
error_array_extension_03 <- results_array_extension_03 + 1


# Calculate Bias and SE of Bias for each scenario and do explorato --------
bias_est_baseline <- apply(error_array_baseline, 1:2, mean)
bias_se_baseline <- sqrt((1/dim(results_array_baseline)[3]) * apply(error_array_baseline, 1:2, var))
matplot(bias_est_baseline, type = "l", lty = 1)
matlines(bias_est_baseline - 2 * bias_se_baseline, type = "l", lty = 2)
matlines(bias_est_baseline +  2 * bias_se_baseline, type = "l", lty = 2)

bias_est_extension_01_low <- apply(error_array_extension_01[,,settings_extension_01$ic_setting_level == 1], 1:2, mean)
bias_se_extension_01_low <- sqrt((1/200) * apply(error_array_extension_01[,,settings_extension_01$ic_setting_level == 1], 1:2, var))
bias_est_extension_01_high <- apply(error_array_extension_01[,,settings_extension_01$ic_setting_level == 2], 1:2, mean)
bias_se_extension_01_high <- sqrt((1/200) * apply(error_array_extension_01[,,settings_extension_01$ic_setting_level == 2], 1:2, var))

par(mfrow = c(1, 2))
matplot(bias_est_extension_01_low, type = "l", lty = 1)
matlines(bias_est_extension_01_low - 2 * bias_se_extension_01_low, type = "l", lty = 2)
matlines(bias_est_extension_01_low + 2 * bias_se_extension_01_low, type = "l", lty = 2)
abline(h = 0, lwd =2 , col="darkgrey")

matplot(bias_est_extension_01_high, type = "l", lty = 1)
matlines(bias_est_extension_01_high - 2 * bias_se_extension_01_low, type = "l", lty = 2)
matlines(bias_est_extension_01_high + 2 * bias_se_extension_01_low, type = "l", lty = 2)
abline(h = 0, lwd =2 , col="darkgrey")


bias_est_extension_02_low <- apply(error_array_extension_02[,,settings_extension_02$intensity_setting_level == 1], 1:2, mean)
bias_se_extension_02_low <- sqrt((1/200) * apply(error_array_extension_02[,,settings_extension_02$intensity_setting_level == 1], 1:2, var))
bias_est_extension_02_high <- apply(error_array_extension_02[,,settings_extension_02$intensity_setting_level == 3], 1:2, mean)
bias_se_extension_02_high <- sqrt((1/200) * apply(error_array_extension_02[,,settings_extension_02$intensity_setting_level == 3], 1:2, var))

par(mfrow = c(1, 2))
matplot(bias_est_extension_02_low, type = "l", lty = 1)
matlines(bias_est_extension_02_low - 2 * bias_se_extension_02_low, type = "l", lty = 2)
matlines(bias_est_extension_02_low + 2 * bias_se_extension_02_low, type = "l", lty = 2)
abline(h = 0, lwd =2 , col="darkgrey")

matplot(bias_est_extension_02_high, type = "l", lty = 1)
matlines(bias_est_extension_02_high - 2 * bias_se_extension_02_low, type = "l", lty = 2)
matlines(bias_est_extension_02_high + 2 * bias_se_extension_02_low, type = "l", lty = 2)
abline(h = 0, lwd =2 , col="darkgrey")


bias_est_extension_03_low <- apply(error_array_extension_03[,,settings_extension_03$sigma_forcing_level == 0.15], 1:2, mean)
bias_se_extension_03_low <- sqrt((1/200) * apply(error_array_extension_03[,,settings_extension_03$sigma_forcing_level == 0.15], 1:2, var))
bias_est_extension_03_high <- apply(error_array_extension_03[,,settings_extension_03$sigma_forcing_level == 0.4], 1:2, mean)
bias_se_extension_03_high <- sqrt((1/200) * apply(error_array_extension_03[,,settings_extension_03$sigma_forcing_level == 0.4], 1:2, var))

par(mfrow = c(1, 2))
matplot(bias_est_extension_03_low, type = "l", lty = 1)
matlines(bias_est_extension_03_low - 2 * bias_se_extension_03_low, type = "l", lty = 2)
matlines(bias_est_extension_03_low + 2 * bias_se_extension_03_low, type = "l", lty = 2)
abline(h = 0, lwd =2 , col="darkgrey")

matplot(bias_est_extension_03_high, type = "l", lty = 1)
matlines(bias_est_extension_03_high - 2 * bias_se_extension_03_low, type = "l", lty = 2)
matlines(bias_est_extension_03_high + 2 * bias_se_extension_03_low, type = "l", lty = 2)
abline(h = 0, lwd =2 , col="darkgrey")





# ------------------------------------------------------------------------#
# ------------------------------------------------------------------------#
# Make Publication-Quality Plots: -----------------------------------------
# ------------------------------------------------------------------------#
# ------------------------------------------------------------------------#


# Baseline: ---------------------------------------------------------------
baseline_dt_est <- data.table(t = grid_points,
                          bias = bias_est_baseline)
baseline_dt_est_lng <- melt.data.table(baseline_dt_est,
                                       id.vars = c("t"),
                                       variable.name = "iteration",
                                       variable.factor = FALSE,
                                       value.name = "bias_est")
baseline_dt_est_lng[, iteration := as.numeric(stringr::str_remove(iteration, "bias.V"))]

baseline_dt_se <- data.table(t = grid_points,
                              bias_se = bias_se_baseline)
baseline_dt_se_lng <- melt.data.table(baseline_dt_se,
                                       id.vars = c("t"),
                                       variable.factor = FALSE,
                                       variable.name = "iteration",
                                       value.name = "bias_se")
baseline_dt_se_lng[, iteration := as.numeric(stringr::str_remove(iteration, "bias_se.V"))]
baseline_dt <- merge.data.table(x = baseline_dt_est_lng, y = baseline_dt_se_lng,
                                by = c("t", "iteration"))
baseline_dt[, ic_level := "baseline"]



# 01 - initial conditions -------------------------------------------------

extension_01_low_dt_est <- data.table(t = grid_points,
                              bias = bias_est_extension_01_low)
extension_01_low_dt_est_lng <- melt.data.table(extension_01_low_dt_est,
                                       id.vars = c("t"),
                                       variable.name = "iteration",
                                       variable.factor = FALSE,
                                       value.name = "bias_est")
extension_01_low_dt_est_lng[, iteration := as.numeric(stringr::str_remove(iteration, "bias.V"))]


extension_01_low_dt_se <- data.table(t = grid_points,
                             bias_se = bias_se_extension_01_low)
extension_01_low_dt_se_lng <- melt.data.table(extension_01_low_dt_se,
                                      id.vars = c("t"),
                                      variable.factor = FALSE,
                                      variable.name = "iteration",
                                      value.name = "bias_se")
extension_01_low_dt_se_lng[, iteration := as.numeric(stringr::str_remove(iteration, "bias_se.V"))]
extension_01_low_dt <- merge.data.table(x = extension_01_low_dt_est_lng, y = extension_01_low_dt_se_lng,
                                by = c("t", "iteration"))
extension_01_low_dt[, ic_level := "low"]


extension_01_high_dt_est <- data.table(t = grid_points,
                                      bias = bias_est_extension_01_high)
extension_01_high_dt_est_lng <- melt.data.table(extension_01_high_dt_est,
                                               id.vars = c("t"),
                                               variable.name = "iteration",
                                               variable.factor = FALSE,
                                               value.name = "bias_est")
extension_01_high_dt_est_lng[, iteration := as.numeric(stringr::str_remove(iteration, "bias.V"))]


extension_01_high_dt_se <- data.table(t = grid_points,
                                     bias_se = bias_se_extension_01_high)
extension_01_high_dt_se_lng <- melt.data.table(extension_01_high_dt_se,
                                              id.vars = c("t"),
                                              variable.factor = FALSE,
                                              variable.name = "iteration",
                                              value.name = "bias_se")
extension_01_high_dt_se_lng[, iteration := as.numeric(stringr::str_remove(iteration, "bias_se.V"))]
extension_01_high_dt <- merge.data.table(x = extension_01_high_dt_est_lng, y = extension_01_high_dt_se_lng,
                                        by = c("t", "iteration"))
extension_01_high_dt[, ic_level := "high"]


extension_01_dt <- rbind(
  baseline_dt, extension_01_low_dt, extension_01_high_dt
)


extension_01_dt[, bias_pw_lower := bias_est - 2 * bias_se]
extension_01_dt[, bias_pw_upper := bias_est + 2 * bias_se]
extension_01_dt[, iteration_coded := factor(iteration, levels = 1:4,
                                            labels = c("No correction",
                                                       "1 iteration",
                                                       paste0(2:3, " iterations")))]
extension_01_dt[, ic_level_label := factor(
  ic_level,
  levels = c("low", "baseline", "high"),
  labels = c("\\textbf{(a)}: $\\boldsymbol{\\mu}_0^\\top = (1, 0)^\\top$",
             "\\textbf{Baseline}: $\\boldsymbol{\\mu}_0^\\top = (0, 0)^\\top$",
             "\\textbf{(b)}: $\\boldsymbol{\\mu}_0^\\top = (0, 1)^\\top$")
)]



p01 <- ggplot(data=extension_01_dt) +
  theme_bw() +
  aes(x=t, colour = iteration_coded, group = iteration) +
  facet_wrap(~ ic_level_label, nrow = 1, ncol = 3) +
  geom_line(aes(y = bias_est))   +
  geom_line(aes(y = bias_pw_lower), lty = 2) +
  geom_line(aes(y = bias_pw_upper), lty = 2) +
  scale_x_continuous(breaks = c(0, pi, 2 * pi), 
                     labels = c("$0$", "$\\pi$", "$2 \\pi$"),
                     limits = c(0, 2 * pi), expand = rep(0.025, 2)) +
  labs(x = "$t$",
       title = "1. Initial Conditions",
       y = "$E[\\hat{\\beta}_0 (t)] - \\beta_0 (t)$")+
  theme(axis.title = element_text(size = 9),
        legend.title = element_blank(),
        legend.position = "none", # FOR NOW BUT WILL USED COMBINED LEGEND
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.title = element_text(size = 10, hjust = 0.5, face = "bold")) +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) 
p01


# tikz(here::here("outputs", "SHM", "paper-plots", "test-simulations-6.1.tex"), 
#      width = (1 * doc_width_inches),
#      height = (0.4 * doc_width_inches))
# print(p)
# dev.off()






# 02 Forcing Intensity ----------------------------------------------------
extension_02_low_dt_est <- data.table(t = grid_points,
                                      bias = bias_est_extension_02_low)
extension_02_low_dt_est_lng <- melt.data.table(extension_02_low_dt_est,
                                               id.vars = c("t"),
                                               variable.name = "iteration",
                                               variable.factor = FALSE,
                                               value.name = "bias_est")
extension_02_low_dt_est_lng[, iteration := as.numeric(stringr::str_remove(iteration, "bias.V"))]


extension_02_low_dt_se <- data.table(t = grid_points,
                                     bias_se = bias_se_extension_02_low)
extension_02_low_dt_se_lng <- melt.data.table(extension_02_low_dt_se,
                                              id.vars = c("t"),
                                              variable.factor = FALSE,
                                              variable.name = "iteration",
                                              value.name = "bias_se")
extension_02_low_dt_se_lng[, iteration := as.numeric(stringr::str_remove(iteration, "bias_se.V"))]
extension_02_low_dt <- merge.data.table(x = extension_02_low_dt_est_lng, y = extension_02_low_dt_se_lng,
                                        by = c("t", "iteration"))
extension_02_low_dt[, intensity_level := 1]


extension_02_high_dt_est <- data.table(t = grid_points,
                                       bias = bias_est_extension_02_high)
extension_02_high_dt_est_lng <- melt.data.table(extension_02_high_dt_est,
                                                id.vars = c("t"),
                                                variable.name = "iteration",
                                                variable.factor = FALSE,
                                                value.name = "bias_est")
extension_02_high_dt_est_lng[, iteration := as.numeric(stringr::str_remove(iteration, "bias.V"))]


extension_02_high_dt_se <- data.table(t = grid_points,
                                      bias_se = bias_se_extension_02_high)
extension_02_high_dt_se_lng <- melt.data.table(extension_02_high_dt_se,
                                               id.vars = c("t"),
                                               variable.factor = FALSE,
                                               variable.name = "iteration",
                                               value.name = "bias_se")
extension_02_high_dt_se_lng[, iteration := as.numeric(stringr::str_remove(iteration, "bias_se.V"))]
extension_02_high_dt <- merge.data.table(x = extension_02_high_dt_est_lng, y = extension_02_high_dt_se_lng,
                                         by = c("t", "iteration"))
extension_02_high_dt[, intensity_level := 3]

baseline_dt[, intensity_level := 2]
extension_02_dt <- rbind(
  baseline_dt[, -c("ic_level")], extension_02_low_dt, extension_02_high_dt
)


extension_02_dt[, bias_pw_lower := bias_est - 2 * bias_se]
extension_02_dt[, bias_pw_upper := bias_est + 2 * bias_se]
extension_02_dt[, iteration_coded := factor(iteration, levels = 1:4,
                                            labels = c("No correction",
                                                       "1 iteration",
                                                       paste0(2:3, " iterations")))]
extension_02_dt[, intensity_level_label := factor(
  intensity_level,
  levels = c(1:3),
  labels = c("\\textbf{(a)}: $l=1$",
             "\\textbf{Baseline}: $l=2$",
             "\\textbf{(b)}: $l=3$")
)]

p02 <- ggplot(data=extension_02_dt) +
  theme_bw() +
  aes(x=t, colour = iteration_coded, group = iteration) +
  facet_wrap(~ intensity_level_label, nrow = 1, ncol = 3) +
  geom_line(aes(y = bias_est))   +
  geom_line(aes(y = bias_pw_lower), lty = 2) +
  geom_line(aes(y = bias_pw_upper), lty = 2) +
  scale_x_continuous(breaks = c(0, pi, 2 * pi), 
                     labels = c("$0$", "$\\pi$", "$2 \\pi$"),
                     limits = c(0, 2 * pi), expand = rep(0.025, 2)) +
  labs(x = "$t$",
       title = "2. Lengthscale of Stochastic Disturbance",
       y = "$E[\\hat{\\beta}_0 (t)] - \\beta_0 (t)$")+
  theme(axis.title = element_text(size = 9),
        legend.title = element_blank(),
        legend.position = "none", # FOR NOW BUT WILL USED COMBINED LEGEND
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.title = element_text(size = 10, hjust = 0.5, face = "bold")) +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) 
p02


# tikz(here::here("outputs", "SHM", "paper-plots", "test-simulations-6.1-02.tex"), 
#      width = (1 * doc_width_inches),
#      height = (0.4 * doc_width_inches))
# print(p)
# dev.off()






# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
extension_03_low_dt_est <- data.table(t = grid_points,
                                      bias = bias_est_extension_03_low)
extension_03_low_dt_est_lng <- melt.data.table(extension_03_low_dt_est,
                                               id.vars = c("t"),
                                               variable.name = "iteration",
                                               variable.factor = FALSE,
                                               value.name = "bias_est")
extension_03_low_dt_est_lng[, iteration := as.numeric(stringr::str_remove(iteration, "bias.V"))]


extension_03_low_dt_se <- data.table(t = grid_points,
                                     bias_se = bias_se_extension_03_low)
extension_03_low_dt_se_lng <- melt.data.table(extension_03_low_dt_se,
                                              id.vars = c("t"),
                                              variable.factor = FALSE,
                                              variable.name = "iteration",
                                              value.name = "bias_se")
extension_03_low_dt_se_lng[, iteration := as.numeric(stringr::str_remove(iteration, "bias_se.V"))]
extension_03_low_dt <- merge.data.table(x = extension_03_low_dt_est_lng, y = extension_03_low_dt_se_lng,
                                        by = c("t", "iteration"))
extension_03_low_dt[, sigma_forcing := 0.15]


extension_03_high_dt_est <- data.table(t = grid_points,
                                       bias = bias_est_extension_03_high)
extension_03_high_dt_est_lng <- melt.data.table(extension_03_high_dt_est,
                                                id.vars = c("t"),
                                                variable.name = "iteration",
                                                variable.factor = FALSE,
                                                value.name = "bias_est")
extension_03_high_dt_est_lng[, iteration := as.numeric(stringr::str_remove(iteration, "bias.V"))]


extension_03_high_dt_se <- data.table(t = grid_points,
                                      bias_se = bias_se_extension_03_high)
extension_03_high_dt_se_lng <- melt.data.table(extension_03_high_dt_se,
                                               id.vars = c("t"),
                                               variable.factor = FALSE,
                                               variable.name = "iteration",
                                               value.name = "bias_se")
extension_03_high_dt_se_lng[, iteration := as.numeric(stringr::str_remove(iteration, "bias_se.V"))]
extension_03_high_dt <- merge.data.table(x = extension_03_high_dt_est_lng, y = extension_03_high_dt_se_lng,
                                         by = c("t", "iteration"))
extension_03_high_dt[, sigma_forcing := 0.4]

baseline_dt[, sigma_forcing := 0.25]

extension_03_dt <- rbind(
  baseline_dt[, -c("ic_level", "intensity_level")], extension_03_low_dt, extension_03_high_dt
)


extension_03_dt[, bias_pw_lower := bias_est - 2 * bias_se]
extension_03_dt[, bias_pw_upper := bias_est + 2 * bias_se]
extension_03_dt[, iteration_coded := factor(iteration, levels = 1:4,
                                            labels = c("No correction",
                                                       "1 iteration",
                                                       paste0(2:3, " iterations")))]
extension_03_dt[, sigma_forcing_label := factor(
  sigma_forcing,
  levels = c(0.15, 0.25, 0.4),
  labels = c("\\textbf{(a)}: $\\sigma=0.15$",
             "\\textbf{Baseline}: $\\sigma=0.25$",
             "\\textbf{(b)}: $\\sigma=0.4$")
)]





# ------------------------------------------------------------------------#
# ------------------------------------------------------------------------#
# ------------------------------------------------------------------------#




theme_set(theme_bw())
theme_update(axis.title = element_text(size = 7.5),
             strip.text = element_text(size = 7.5),
             legend.text = element_text(size = 7),
             axis.text = element_text(size = 7),
             legend.title = element_blank(),
             legend.position = "none", # FOR NOW BUT WILL USED COMBINED LEGEND
             axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
             plot.title = element_text(size = 9, hjust = 0.5, face = "bold"))

p01 <- ggplot(data=extension_01_dt) +
  aes(x=t, colour = iteration_coded, group = iteration) +
  facet_wrap(~ ic_level_label, nrow = 1, ncol = 3) +
  geom_line(aes(y = bias_est))   +
  geom_line(aes(y = bias_pw_lower), lty = 2) +
  geom_line(aes(y = bias_pw_upper), lty = 2) +
  scale_x_continuous(breaks = c(0, pi, 2 * pi), 
                     labels = c("$0$", "$\\pi$", "$2 \\pi$"),
                     limits = c(0, 2 * pi), expand = rep(0.025, 2)) +
  labs(x = "$t$",
       title = "1. Initial Conditions",
       y = "$\\mathbb{E}[\\hat{\\beta}_0 (t)] - \\beta_0 (t)$") +
  guides(colour = guide_legend(override.aes = list(size = 2))) 
p01


p02 <- ggplot(data=extension_02_dt) +
  aes(x=t, colour = iteration_coded, group = iteration) +
  facet_wrap(~ intensity_level_label, nrow = 1, ncol = 3) +
  geom_line(aes(y = bias_est))   +
  geom_line(aes(y = bias_pw_lower), lty = 2) +
  geom_line(aes(y = bias_pw_upper), lty = 2) +
  scale_x_continuous(breaks = c(0, pi, 2 * pi), 
                     labels = c("$0$", "$\\pi$", "$2 \\pi$"),
                     limits = c(0, 2 * pi), expand = rep(0.025, 2)) +
  labs(x = "$t$",
       title = "2. Lengthscale of Stochastic Disturbance",
       y = "$\\mathbb{E}[\\hat{\\beta}_0 (t)] - \\beta_0 (t)$")+
  guides(colour = guide_legend(override.aes = list(size = 2))) 
p02

p03 <- ggplot(data=extension_03_dt) +
  aes(x=t, colour = iteration_coded, group = iteration) +
  facet_wrap(~ sigma_forcing_label, nrow = 1, ncol = 3) +
  geom_line(aes(y = bias_est))   +
  geom_line(aes(y = bias_pw_lower), lty = 2) +
  geom_line(aes(y = bias_pw_upper), lty = 2) +
  scale_x_continuous(breaks = c(0, pi, 2 * pi), 
                     labels = c("$0$", "$\\pi$", "$2 \\pi$"),
                     limits = c(0, 2 * pi), expand = rep(0.025, 2)) +
  labs(x = "$t$",
       title = "3. Amplitude of Stochastic Disturbance",
       y = "$\\mathbb{E}[\\hat{\\beta}_0 (t)] - \\beta_0 (t)$")+
  guides(colour = guide_legend(override.aes = list(size = 2))) 
p03

ggarrange(p01, p02, p03, ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom")
# tikz(here::here("outputs", "SHM", "paper-plots", "test-simulations-6.1-03.tex"), 
#      width = (1 * doc_width_inches),
#      height = (0.4 * doc_width_inches))
# print(p)
# dev.off()



# Combine plots: ----------------------------------------------------------

tikz(here::here("outputs", "SHM", "paper-plots", "simulation-results-6.1.tex"),
     width = (0.7 * doc_width_inches),
     height = (1 * doc_width_inches))
ggarrange(p01, p02, p03, ncol = 1, nrow = 3, common.legend = TRUE, legend = "bottom") +
  guides(colour = guide_legend(override.aes = list(size = 2)))
dev.off()


