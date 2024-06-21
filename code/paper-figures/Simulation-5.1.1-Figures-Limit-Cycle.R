# ------------------------------------------------------------------------#
# at sigma = 0.4, the bias correction seemingly doesn't work well.
# in that the estimated parameters do begin failing in resembling
# the Jacobian of the system about the limit cycle.
# There are a number of reasons for this:
# 1) the bias simply increased by increasing sigma,
# and is hence harder to correct
# 2) the linearised approximation is less likely to be valid for
# larger perturbations
# 3) due to phase variation introduced the mean is no longer a point on the
# limit cycle
# we illustrate the third point in this script
# ------------------------------------------------------------------------#

# Packages and Functions: -------------------------------------------------
source(here::here("code", "VDP-helper-functions.R"))
library(data.table)  # CRAN v1.14.2
library(ggplot2)     # CRAN v3.4.0
library(tikzDevice)  # CRAN v0.12.3.1
library(purrr)       # CRAN v0.3.4

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
                    "vdp-simulation-section-5.1-vary-sigma-forcing-updated.rds"))
grid_points <- simulation_results$grid_points
results_list <- simulation_results$results_list
mu <- simulation_results$mu
settings <- simulation_results$settings

# Parameters that will be fixed: ------------------------------------------
mu_ic <- c(1.9922213675594, -0.910974076470711)
Sigma_ic <- matrix(c(0.025, 0, 0, 0.025),
                   nrow = 2, 
                   ncol = 2) # Not too much IC noise


# -------------------------------------------------------------------------
# Construct ground truth to compare to by samplign data points 
# determinsitically from limit cycle
# -------------------------------------------------------------------------
limit_cycle_data <- generate_vdp_deterministic(mu = 1, grid_points = grid_points, y0 = mu_ic[2], x0 = mu_ic[1])
limit_x <- limit_cycle_data[,2]
limit_y <- limit_cycle_data[,3]

b_0x_limit <- mu * (limit_x - (1/3)*limit_x^3 - limit_y)  - mu * (1 - limit_x^2) * limit_x + mu * limit_y
beta_xx_limit <- mu * (1-limit_x^2)

truth_dt <- data.table(
  t = grid_points,
  b0_x = b_0x_limit,
  b0_y = 0,
  beta_xx = beta_xx_limit,
  beta_yx = 1/mu,
  beta_xy = - mu,
  beta_yy = 0)
truth_dt_lng <- melt.data.table(truth_dt, 
                                id.vars = "t",
                                variable.name = "parameter",
                                variable.factor = FALSE, 
                                value.factor = FALSE,
                                value.name = "value")

ggplot(truth_dt_lng) +
  aes(x = t, y = value) +
  facet_wrap(~ parameter, scales = "free") +
  geom_line()
dev.off()



# -------------------------------------------------------------------------
# calculate mean x and y empircally using a very big sample:
mu_ic <- c(1.9922213675594, -0.910974076470711)
Sigma_ic <- matrix(c(0.025, 0, 0, 0.025),
                   nrow = 2, 
                   ncol = 2) # Not too much IC noise
n_for_example <- 10000 # large N
data_set_sigma_0.4 <- generate_vdp_curves_stochastic_ic_and_forcing(n = n_for_example,
                                                                    mu = mu,
                                                                    mu_ic = mu_ic,
                                                                    grid_points = grid_points,
                                                                    Sigma_ic = Sigma_ic,
                                                                    sigma_forcing = 0.4,
                                                                    intensity = 2)
dim(data_set_sigma_0.4$data)
limit_x_hat <- apply(data_set_sigma_0.4$data[,1,], 1, mean)
limit_y_hat <- apply(data_set_sigma_0.4$data[,2,], 1, mean)

plot(limit_y_hat)
lines(limit_y)

plot(limit_x_hat)
lines(limit_x)

beta_xx_hat <- (1-limit_x_hat^2)
empirical_dt <- data.table(t = grid_points,
                           beta_xx = beta_xx_hat)


inds_sigma_0.4 <- which(settings$sigma_forcing == 0.4)
results_sigma_0.4 <- results_list[inds_sigma_0.4]
sum(sapply(results_sigma_0.4, function(x){is.null(x$beta_array)})) # Need tp check these out:
estimates_sigma_0.4 <- map_dfr(.x = results_sigma_0.4, function(x) {
  
  if(is.null(x$beta_array)) {
    NULL
  } else{
    ols <- data.frame(t = grid_points,
                      cbind(x$beta_array[,,1,1],
                            x$beta_array[,,2,1]),
                      estimator = "ols")
    bc <- data.frame(t = grid_points,
                     cbind(x$beta_array[,,1,11],
                           x$beta_array[,,2,11]),
                     estimator = "bc")
    full <- rbind(ols, bc)
    names(full) <- c("t", "b0_x", "beta_xx", "beta_xy", "b0_y", "beta_yx", "beta_yy", "estimator")
    full
  }
  
}, .id = "sim_rep")

estimates_sigma_0.4 <- data.table(estimates_sigma_0.4)
estimates_sigma_0.4[, estimator := factor(estimator, levels = c("ols", "bc"), labels = c("OLS", "Bias-Corrected"))]
# -------------------------------------------------------------------------
head(estimates_sigma_0.4)


limit_cycle_plot <- ggplot(data = estimates_sigma_0.4) +
  aes(x = t, y = beta_xx) +
  geom_line(aes(group = interaction(sim_rep, estimator), color = estimator), alpha = 0.1) +
  geom_line(data = truth_dt, aes(linetype = "Limit Cycle")) +
  geom_line(data = empirical_dt, aes(linetype = "Empirical")) +
  theme(legend.position = "bottom",
        legend.title = element_blank()) +
  scale_linetype_manual(values = c(2, 1)) +
  ylim(c(-3.5, 1.25)) +
  labs(x = "$t$", y = "$\\beta_{xx}(t)$") +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 0.75))) +
  scale_color_manual(values = c("#CC79A7", "#009E73")) 


tikz(here::here("outputs", "VDP", "paper-plots", "vary-sigma-limit_cycle_plot.tex"),
     width = 1 * doc_width_inches, 
     height =  0.85 * doc_width_inches)
limit_cycle_plot
dev.off()


