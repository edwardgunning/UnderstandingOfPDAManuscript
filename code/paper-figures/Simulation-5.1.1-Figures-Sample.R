# ------------------------------------------------------------------------#
# Plots of example datasets with varying sigma parameters and also 
# the associated signal to noise ratio.
# ------------------------------------------------------------------------#

# Packages and Functions: -------------------------------------------------
source(here::here("code", "VDP-helper-functions.R"))
library(data.table)  # CRAN v1.14.2
library(ggplot2)     # CRAN v3.4.0
library(tikzDevice)  # CRAN v0.12.3.1
library(randomcoloR) # CRAN v1.1.0.1

# Plot settings -----------------------------------------------------------
source(here::here("code", "theme_gunning.R"))
theme_gunning()
theme_update(legend.position = "bottom", 
             legend.title = element_text(face = "bold", size = 13),
             strip.text = element_text(size = 12),
             axis.text = element_text(size = 11),
             axis.title = element_text(size = 12),
             legend.text = element_text(size = 13))
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937


# Import data: ------------------------------------------------------------
simulation_results <- readRDS(
  file = here::here("outputs", "VDP", 
                    "vdp-simulation-section-5.1-vary-sigma-forcing.rds"))
grid_points <- simulation_results$grid_points
results_list <- simulation_results$results_list
mu <- simulation_results$mu
settings <- simulation_results$settings


unique(settings$sigma_forcing)
                     # [1] 0.1 0.2 0.4

# Generate Example Datasets: ----------------------------------------------
# we want tp show what datasets "look" like:
mu_ic <- c(1.9922213675594, -0.910974076470711)
Sigma_ic <- matrix(c(0.025, 0, 0, 0.025),
                   nrow = 2, 
                   ncol = 2) # Not too much IC noise
n_for_example <- 20

set.seed(1)
data_set_sigma_0.1 <- generate_vdp_curves_stochastic_ic_and_forcing(n = n_for_example,
                                              mu = mu,
                                              mu_ic = mu_ic,
                                              grid_points = grid_points,
                                              Sigma_ic = Sigma_ic,
                                              sigma_forcing = 0.1,
                                              intensity = 2)

data_set_sigma_0.2 <- generate_vdp_curves_stochastic_ic_and_forcing(n = n_for_example,
                                                                    mu = mu,
                                                                    mu_ic = mu_ic,
                                                                    grid_points = grid_points,
                                                                    Sigma_ic = Sigma_ic,
                                                                    sigma_forcing = 0.2,
                                                                    intensity = 2)
data_set_sigma_0.4 <- generate_vdp_curves_stochastic_ic_and_forcing(n = n_for_example,
                                                                    mu = mu,
                                                                    mu_ic = mu_ic,
                                                                    grid_points = grid_points,
                                                                    Sigma_ic = Sigma_ic,
                                                                    sigma_forcing = 0.4,
                                                                    intensity = 2)

par(mfrow = c(2, 3))
matplot(data_set_sigma_0.1$data[, 1,], type = "l")
matplot(data_set_sigma_0.2$data[, 1,], type = "l")
matplot(data_set_sigma_0.4$data[, 1,], type = "l")

matplot(data_set_sigma_0.1$data[, 2,], type = "l")
matplot(data_set_sigma_0.2$data[, 2,], type = "l")
matplot(data_set_sigma_0.4$data[, 2,], type = "l")




plot_datasets_x_dt <- data.table(
  t = rep(grid_points, times = 3),
  sigma = rep(c(0.1, 0.2, 0.4), each = length(grid_points)),
  x = rbind(
    data_set_sigma_0.1$data[,1,],
    data_set_sigma_0.2$data[,1,],
    data_set_sigma_0.4$data[,1,]
  ))

names(plot_datasets_x_dt)[-c(1:2)] <- paste0("obs_", seq_len(n_for_example))

plot_datasets_x_dt_lng <- melt.data.table(data = plot_datasets_x_dt,
                                          id.vars = c("t", "sigma"),
                                          measure.vars = paste0("obs_", seq_len(n_for_example)),
                                          variable.factor = FALSE, 
                                          value.factor = FALSE,
                                          variable.name = "obs",
                                          value.name = "x")

plot_datasets_y_dt <- data.table(
  t = rep(grid_points, times = 3),
  sigma = rep(c(0.1, 0.2, 0.4), each = length(grid_points)),
  x = rbind(
    data_set_sigma_0.1$data[,2,],
    data_set_sigma_0.2$data[,2,],
    data_set_sigma_0.4$data[,2,]
  ))

names(plot_datasets_y_dt)[-c(1:2)] <- paste0("obs_", seq_len(n_for_example))

plot_datasets_y_dt_lng <- melt.data.table(data = plot_datasets_y_dt,
                                          id.vars = c("t", "sigma"),
                                          measure.vars = paste0("obs_", seq_len(n_for_example)),
                                          variable.factor = FALSE, 
                                          value.factor = FALSE,
                                          variable.name = "obs",
                                          value.name = "y")

plot_datasets_dt <- merge.data.table(x = plot_datasets_x_dt_lng,
                                     y = plot_datasets_y_dt_lng, 
                                     by = c("t", "sigma", "obs"), 
                                     all = TRUE
                                     )

plot_datasets_dt_lng <- melt.data.table(data = plot_datasets_dt, 
                                        id.vars =  c("t", "sigma", "obs"),
                                        measure.vars = c("x", "y"),
                                        variable.name = "dimension",
                                        value.name = "value"
                                        )
plot_datasets_dt_lng[, sigma := factor(sigma, 
                                               levels = c(0.1, 0.2, 0.4),
                                               labels = paste0("$\\sigma = ", c(0.1, 0.2, 0.4), "$"))]

plot_datasets_dt_lng[, dimension := factor(dimension, 
                                           levels = c("x", "y"),
                                           labels = c("$x(t)$", "$y(t)$"))]

# set.seed(2)
# palette <- randomColor(n_for_example) # simulate random colours to plot each curve
# dput(palette)
palette <- c("#6ae2d2", "#96ce42", "#e9c6ff", "#c3a9f9", "#e58edc", "#9185e5",
"#ed8c6f", "#e7c9ff", "#d9f79e", "#1bc6ad", "#a899e8", "#ffd7cc",
"#70ead0", "#fcc962", "#477ea5", "#baffaf", "#cc005b", "#a9dd37",
"#fcc2c5", "#151291")
p <- ggplot(plot_datasets_dt_lng) +
  aes(x = t, group = obs, y = value, color = obs) +
  facet_wrap(dimension ~ sigma, nrow = 2, ncol = 3) +
  geom_line() +
  labs(x = "$t$", y = "Function Value") +
  scale_color_manual(guide = "none", values = palette)

p

tikz(here::here("outputs", "VDP", "paper-plots", "vary-sigma-dataset-demo.tex"),
     width = (3/2) * doc_width_inches, 
     height =  doc_width_inches)
p
dev.off()



# Plot of Signal to Noise Ratio:  ------------------------------------------

x_sigma_0.1 <- data_set_sigma_0.1$data[,1,]
y_sigma_0.1 <- data_set_sigma_0.1$data[,2,]
eps_x_sigma_0.1 <- data_set_sigma_0.1$forcing[,1,]
eps_y_sigma_0.1 <- data_set_sigma_0.1$forcing[,2,]
dx_sigma_0.1_linear_pred <- mu * (x_sigma_0.1 - (1/3) * (x_sigma_0.1^3) - y_sigma_0.1)
dy_sigma_0.1_linear_pred <- (1 / mu) * x_sigma_0.1



x_sigma_0.2 <- data_set_sigma_0.2$data[,1,]
y_sigma_0.2 <- data_set_sigma_0.2$data[,2,]
eps_x_sigma_0.2 <- data_set_sigma_0.2$forcing[,1,]
eps_y_sigma_0.2 <- data_set_sigma_0.2$forcing[,2,]
dx_sigma_0.2_linear_pred <- mu * (x_sigma_0.2 - (1/3) * (x_sigma_0.2^3) - y_sigma_0.2)
dy_sigma_0.2_linear_pred <- (1 / mu) * x_sigma_0.2

x_sigma_0.4 <- data_set_sigma_0.4$data[,1,]
y_sigma_0.4 <- data_set_sigma_0.4$data[,2,]
eps_x_sigma_0.4 <- data_set_sigma_0.4$forcing[,1,]
eps_y_sigma_0.4 <- data_set_sigma_0.4$forcing[,2,]
dx_sigma_0.4_linear_pred <- mu * (x_sigma_0.4 - (1/3) * (x_sigma_0.4^3) - y_sigma_0.4)
dy_sigma_0.4_linear_pred <- (1 / mu) * x_sigma_0.4



#                    # Create datasets for plotting -------------------------------------------


### Dx Linear Predictor: --------------------------------------------------
plot_predictor_dx_dt <- data.table(
  t = rep(grid_points, times = 3),
  sigma = rep(c(0.1, 0.2, 0.4), each = length(grid_points)),
  dx_pred = rbind(
    dx_sigma_0.1_linear_pred,
    dx_sigma_0.2_linear_pred,
    dx_sigma_0.4_linear_pred
  ))

names(plot_predictor_dx_dt)[-c(1:2)] <- paste0("obs_", seq_len(n_for_example))

plot_predictor_dx_dt_lng <- melt.data.table(data = plot_predictor_dx_dt,
                                          id.vars = c("t", "sigma"),
                                          measure.vars = paste0("obs_", seq_len(n_for_example)),
                                          variable.factor = FALSE, 
                                          value.factor = FALSE,
                                          variable.name = "obs",
                                          value.name = "dx_pred")

plot_predictor_dy_dt <- data.table(
  t = rep(grid_points, times = 3),
  sigma = rep(c(0.1, 0.2, 0.4), each = length(grid_points)),
  dy_pred = rbind(
    dy_sigma_0.1_linear_pred,
    dy_sigma_0.2_linear_pred,
    dy_sigma_0.4_linear_pred
  ))

names(plot_predictor_dy_dt)[-c(1:2)] <- paste0("obs_", seq_len(n_for_example))


plot_predictor_dy_dt_lng <- melt.data.table(data = plot_predictor_dy_dt,
                                            id.vars = c("t", "sigma"),
                                            measure.vars = paste0("obs_", seq_len(n_for_example)),
                                            variable.factor = FALSE, 
                                            value.factor = FALSE,
                                            variable.name = "obs",
                                            value.name = "dy_pred")

plot_predictor_dt <- merge.data.table(x = plot_predictor_dx_dt_lng,
                                     y = plot_predictor_dy_dt_lng, 
                                     by = c("t", "sigma", "obs"), 
                                     all = TRUE)

plot_predictor_dt_lng <- melt.data.table(data = plot_predictor_dt, 
                                        id.vars =  c("t", "sigma", "obs"),
                                        measure.vars = c("dx_pred", "dy_pred"),
                                        variable.name = "dimension",
                                        value.name = "value")
plot_predictor_dt_lng[, sigma := factor(sigma, 
                                       levels = c(0.1, 0.2, 0.4),
                                       labels = paste0("$\\sigma = ", c(0.1, 0.2, 0.4), "$"))]

plot_predictor_dt_lng[, dimension := factor(dimension, 
                                           levels = c("dx_pred", "dy_pred"),
                                           labels = c("$Dx(t)$", "$Dy(t)$"))]



# Now plot generated noise: -----------------------------------------------

plot_eps_dx_dt <- data.table(
  t = rep(grid_points, times = 3),
  sigma = rep(c(0.1, 0.2, 0.4), each = length(grid_points)),
  rbind(
    eps_x_sigma_0.1,
    eps_x_sigma_0.2,
    eps_x_sigma_0.4
  ))

names(plot_eps_dx_dt)[-c(1:2)] <- paste0("obs_", seq_len(n_for_example))

plot_eps_dx_dt_lng <- melt.data.table(data = plot_eps_dx_dt,
                                            id.vars = c("t", "sigma"),
                                            measure.vars = paste0("obs_", seq_len(n_for_example)),
                                            variable.factor = FALSE, 
                                            value.factor = FALSE,
                                            variable.name = "obs",
                                            value.name = "eps_x")

plot_eps_dy_dt <- data.table(
  t = rep(grid_points, times = 3),
  sigma = rep(c(0.1, 0.2, 0.4), each = length(grid_points)),
  rbind(
    eps_y_sigma_0.1,
    eps_y_sigma_0.2,
    eps_y_sigma_0.4
  ))

names(plot_eps_dy_dt)[-c(1:2)] <- paste0("obs_", seq_len(n_for_example))

plot_eps_dy_dt_lng <- melt.data.table(data = plot_eps_dy_dt,
                                      id.vars = c("t", "sigma"),
                                      measure.vars = paste0("obs_", seq_len(n_for_example)),
                                      variable.factor = FALSE, 
                                      value.factor = FALSE,
                                      variable.name = "obs",
                                      value.name = "eps_y")

plot_eps_dt <- merge.data.table(x = plot_eps_dx_dt_lng,
                                      y = plot_eps_dy_dt_lng, 
                                      by = c("t", "sigma", "obs"), 
                                      all = TRUE)

plot_eps_dt_lng <- melt.data.table(data = plot_eps_dt, 
                                         id.vars =  c("t", "sigma", "obs"),
                                         measure.vars = c("eps_x", "eps_y"),
                                         variable.name = "dimension",
                                         value.name = "value")
plot_eps_dt_lng[, sigma := factor(sigma, 
                                        levels = c(0.1, 0.2, 0.4),
                                        labels = paste0("$\\sigma = ", c(0.1, 0.2, 0.4), "$"))]

plot_eps_dt_lng[, dimension := factor(dimension, 
                                            levels = c("eps_x", "eps_y"),
                                            labels = c("$Dx(t)$", "$Dy(t)$"))]


plot_eps_dt_lng$quantity = "Stochastic Disturbance: $\\epsilon (t)$"
plot_predictor_dt_lng$quantity <- "Predictor Function: $g(x(t), y(t))$"
plot_snr_dt <- rbind(plot_eps_dt_lng, plot_predictor_dt_lng)



p2 <- ggplot(plot_snr_dt) +
  aes(x = t, group = interaction(obs, quantity), y = value, color = quantity) +
  facet_wrap(dimension ~ sigma, nrow = 2, ncol = 3) +
  geom_line() +
  labs(x = "$t$", y = "Function Value") +
  scale_color_manual(values = c(
    "#D16103", "#4E84C4"
  )) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.5))) +
  theme(legend.title = element_blank())
p2


tikz(here::here("outputs", "VDP", "paper-plots", "vary-sigma-snr-demo.tex"),
     width = (3/2) * doc_width_inches, 
     height =  1.1 * doc_width_inches)
p2
dev.off()
