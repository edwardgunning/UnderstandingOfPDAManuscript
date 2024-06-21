# ------------------------------------------------------------------------#
# Plots of example datasets with varying 'mu' parameter
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
                    "vdp-simulation-section-5.1-vary-mu.rds"))
grid_points <- simulation_results$grid_points
results_list <- simulation_results$results_list
settings <- simulation_results$settings


unique(settings$intensity_forcing)
unique(settings$intensity_forcing)
unique(settings$mu)


# Generate Example Datasets: ----------------------------------------------
# we want tp show what datasets "look" like:
mu_ic <- c(1.9922213675594, -0.910974076470711)
Sigma_ic <- matrix(c(0.025, 0, 0, 0.025),
                   nrow = 2, 
                   ncol = 2) # Not too much IC noise
n_for_example <- 20

set.seed(1)


Sigma_ic_0.5 <- Sigma_ic
Sigma_ic_0.5[2,2] <- 2^2 * Sigma_ic_0.5[2,2]




Sigma_ic_2 <- Sigma_ic
Sigma_ic_2[2,2] <- 0.5^2 * Sigma_ic_2[2,2]

data_set_mu_0.5 <- generate_vdp_curves_stochastic_ic_and_forcing(n = n_for_example,
                                                              mu = 0.5,
                                                              mu_ic = mu_ic,
                                                              grid_points = grid_points,
                                                              Sigma_ic = Sigma_ic_0.5,
                                                              sigma_forcing = 0.1,
                                                              intensity = 2)

data_set_mu_1 <- generate_vdp_curves_stochastic_ic_and_forcing(n = n_for_example,
                                                              mu = 1,
                                                              mu_ic = mu_ic,
                                                              grid_points = grid_points,
                                                              Sigma_ic = Sigma_ic,
                                                              sigma_forcing = 0.1,
                                                              intensity = 2)
data_set_mu_2 <- generate_vdp_curves_stochastic_ic_and_forcing(n = n_for_example,
                                                              mu = 2,
                                                              mu_ic = mu_ic,
                                                              grid_points = grid_points,
                                                              Sigma_ic = Sigma_ic_2,
                                                              sigma_forcing = 0.1,
                                                              intensity = 2)




par(mfrow = c(2, 3))
matplot(data_set_mu_0.5$data[, 1,], type = "l")
matplot(data_set_mu_1$data[, 1,], type = "l")
matplot(data_set_mu_2$data[, 1,], type = "l")

matplot(data_set_mu_0.5$data[, 2,], type = "l")
matplot(data_set_mu_1$data[, 2,], type = "l")
matplot(data_set_mu_2$data[, 2,], type = "l")



# -------------------------------------------------------------------------


plot_datasets_x_dt <- data.table(
  t = rep(grid_points, times = 3),
  mu = rep(c(0.5, 1, 2), each = length(grid_points)),
  x = rbind(
    data_set_mu_0.5$data[,1,],
    data_set_mu_1$data[,1,],
    data_set_mu_2$data[,1,]
  ))

names(plot_datasets_x_dt)[-c(1:2)] <- paste0("obs_", seq_len(n_for_example))

plot_datasets_x_dt_lng <- melt.data.table(data = plot_datasets_x_dt,
                                          id.vars = c("t", "mu"),
                                          measure.vars = paste0("obs_", seq_len(n_for_example)),
                                          variable.factor = FALSE, 
                                          value.factor = FALSE,
                                          variable.name = "obs",
                                          value.name = "x")

plot_datasets_y_dt <- data.table(
  t = rep(grid_points, times = 3),
  mu = rep(c(0.5, 1, 2), each = length(grid_points)),
  x = rbind(
    data_set_mu_0.5$data[,2,],
    data_set_mu_1$data[,2,],
    data_set_mu_2$data[,2,]
  ))

names(plot_datasets_y_dt)[-c(1:2)] <- paste0("obs_", seq_len(n_for_example))

plot_datasets_y_dt_lng <- melt.data.table(data = plot_datasets_y_dt,
                                          id.vars = c("t", "mu"),
                                          measure.vars = paste0("obs_", seq_len(n_for_example)),
                                          variable.factor = FALSE, 
                                          value.factor = FALSE,
                                          variable.name = "obs",
                                          value.name = "y")

plot_datasets_dt <- merge.data.table(x = plot_datasets_x_dt_lng,
                                     y = plot_datasets_y_dt_lng, 
                                     by = c("t", "mu", "obs"), 
                                     all = TRUE
)

plot_datasets_dt_lng <- melt.data.table(data = plot_datasets_dt, 
                                        id.vars =  c("t", "mu", "obs"),
                                        measure.vars = c("x", "y"),
                                        variable.name = "dimension",
                                        value.name = "value"
)
plot_datasets_dt_lng[, mu := factor(mu, 
                                   levels = c(0.5, 1, 2),
                                   labels = paste0("$\\mu = ", c(0.5, 1, 2), "$"))]

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
  facet_wrap(dimension ~ mu, nrow = 2, ncol = 3, "free_y") +
  geom_line() +
  labs(x = "$t$", y = "Function Value") +
  scale_color_manual(guide = "none", values = palette)

p

tikz(here::here("outputs", "VDP", "paper-plots", "vary-mu-dataset-demo.tex"),
     width = (3/2) * doc_width_inches, 
     height =  doc_width_inches)
p
dev.off()




