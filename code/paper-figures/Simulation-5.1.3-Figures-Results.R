# ------------------------------------------------------------------------#
# Plotting results of simulation changing the sigma parameter in the
# Van der Pol Model.
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
                    "vdp-simulation-section-5.1-vary-mu-updated.rds"))
grid_points <- simulation_results$grid_points
results_list <- simulation_results$results_list
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
limit_cycle_data_mu_0.5 <- generate_vdp_deterministic(mu = 0.5, grid_points = grid_points, y0 = mu_ic[2], x0 = mu_ic[1])
limit_x_mu_0.5 <- limit_cycle_data_mu_0.5[,2]
limit_y_mu_0.5 <- limit_cycle_data_mu_0.5[,3]


b_0x_limit_mu_0.5 <- 0.5 * (limit_x_mu_0.5 - (1/3)*limit_x_mu_0.5^3 - limit_y_mu_0.5)  - 0.5 * (1 - limit_x_mu_0.5^2) * limit_x_mu_0.5 + 0.5 * limit_y_mu_0.5
beta_xx_limit_mu_0.5 <- 0.5 * (1-limit_x_mu_0.5^2)

truth_dt_mu_0.5 <- data.table(
  t = grid_points,
  b0_x = b_0x_limit_mu_0.5,
  b0_y = 0,
  beta_xx = beta_xx_limit_mu_0.5,
  beta_yx = 1/0.5,
  beta_xy = - 0.5,
  beta_yy = 0)

limit_cycle_data_mu_1 <- generate_vdp_deterministic(mu = 1, grid_points = grid_points, y0 = mu_ic[2], x0 = mu_ic[1])
limit_x_mu_1 <- limit_cycle_data_mu_1[,2]
limit_y_mu_1 <- limit_cycle_data_mu_1[,3]


b_0x_limit_mu_1 <- 1 * (limit_x_mu_1 - (1/3)*limit_x_mu_1^3 - limit_y_mu_1)  - 1 * (1 - limit_x_mu_1^2) * limit_x_mu_1 + 1 * limit_y_mu_1
beta_xx_limit_mu_1 <- 1 * (1-limit_x_mu_1^2)

truth_dt_mu_1 <- data.table(
  t = grid_points,
  b0_x = b_0x_limit_mu_1,
  b0_y = 0,
  beta_xx = beta_xx_limit_mu_1,
  beta_yx = 1/1,
  beta_xy = - 1,
  beta_yy = 0)

limit_cycle_data_mu_2 <- generate_vdp_deterministic(mu = 2, grid_points = grid_points, y0 = mu_ic[2], x0 = mu_ic[1])
limit_x_mu_2 <- limit_cycle_data_mu_2[,2]
limit_y_mu_2 <- limit_cycle_data_mu_2[,3]


b_0x_limit_mu_2 <- 2 * (limit_x_mu_2 - (1/3)*limit_x_mu_2^3 - limit_y_mu_2)  - 2 * (1 - limit_x_mu_2^2) * limit_x_mu_2 + 2 * limit_y_mu_2
beta_xx_limit_mu_2 <- 2 * (1-limit_x_mu_2^2)

truth_dt_mu_2 <- data.table(
  t = grid_points,
  b0_x = b_0x_limit_mu_2,
  b0_y = 0,
  beta_xx = beta_xx_limit_mu_2,
  beta_yx = 1/2,
  beta_xy = - 2,
  beta_yy = 0)

truth_dt_mu_0.5$mu <- 0.5
truth_dt_mu_1$mu <- 1
truth_dt_mu_2$mu <- 2

truth_dt <- rbind(truth_dt_mu_0.5,
                  truth_dt_mu_1,
                  truth_dt_mu_2)

truth_dt_lng <- melt.data.table(truth_dt, 
                                id.vars = c("t", "mu"),
                                variable.name = "parameter",
                                variable.factor = FALSE, 
                                value.factor = FALSE,
                                value.name = "value")

ggplot(truth_dt_lng) +
  aes(x = t, y = value) +
  facet_grid(parameter~mu, scales = "free") +
  geom_line()
dev.off()


sum(sapply(results_list, function(x){is.null(x$beta_array)})) # Need tp check these out:


# -------------------------------------------------------------------------

# Sigma = 0.1: ------------------------------------------------------------
inds_mu_0.5 <- which(settings$mu == 0.5)
results_mu_0.5 <- results_list[inds_mu_0.5]
estimates_mu_0.5 <- map_dfr(.x = results_mu_0.5, function(x) {
  
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
}, .id = "sim_rep")

inds_mu_1 <- which(settings$mu == 1)
results_mu_1 <- results_list[inds_mu_1]
estimates_mu_1 <- map_dfr(.x = results_mu_1, function(x) {
  
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
}, .id = "sim_rep")


inds_mu_2 <- which(settings$mu == 2)
results_mu_2 <- results_list[inds_mu_2]

estimates_mu_2 <- map_dfr(.x = results_mu_2, function(x) {
  
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

estimates_mu_0.5$mu <- 0.5
estimates_mu_1$mu <- 1
estimates_mu_2$mu <- 2

estimates_dt <- data.table(rbind(estimates_mu_0.5,
                                 estimates_mu_1,
                                 estimates_mu_2))
estimates_dt_lng <- melt.data.table(estimates_dt,
                                    id.vars = c("t", "mu", "estimator", "sim_rep"),
                                    variable.name = "parameter",
                                    value.name = "value",
                                    variable.factor = FALSE,
                                    value.factor = FALSE,
                                    verbose = TRUE)


estimates_dt_lng[, mu := factor(mu,
                               levels = c(0.5, 1, 2),
                               labels = paste0("$\\mu = ", c(0.5, 1, 2), "$"))]
truth_dt_lng[, mu := factor(mu,
                                levels = c(0.5, 1, 2),
                                labels = paste0("$\\mu = ", c(0.5, 1, 2), "$"))]
estimates_dt_lng[, parameter := factor(parameter,
                                       levels = c("b0_x", "beta_xx", "beta_xy", "b0_y", "beta_yx", "beta_yy"),
                                       labels = c("$\\alpha^{(x)} (t)$",
                                                  "$\\beta_{xx} (t)$",
                                                  "$\\beta_{xy} (t)$",
                                                  "$\\alpha^{(y)} (t)$",
                                                  "$\\beta_{yx} (t)$",
                                                  "$\\beta_{yy} (t)$"))]
truth_dt_lng[, parameter := factor(parameter,
                                   levels = c("b0_x", "beta_xx", "beta_xy", "b0_y", "beta_yx", "beta_yy"),
                                   labels = c("$\\alpha^{(x)} (t)$",
                                              "$\\beta_{xx} (t)$",
                                              "$\\beta_{xy} (t)$",
                                              "$\\alpha^{(y)} (t)$",
                                              "$\\beta_{yx} (t)$",
                                              "$\\beta_{yy} (t)$"))]

estimates_dt_lng[, estimator := factor(estimator, levels = c("ols", "bc"), labels = c("OLS", "Bias-Corrected"))]
head(truth_dt_lng)


p <- ggplot(estimates_dt_lng) +
  aes(x = t, y = value) +
  facet_wrap(parameter~mu, scales = "free_y", ncol =3) +
  geom_line(aes( group = interaction(estimator, sim_rep), color = estimator), alpha = 0.4) +
  geom_line(data = truth_dt_lng, linewidth = 1) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.y = element_blank()) +
  labs(x = "$t$") +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 0.75))) +
  scale_color_manual(values = c("#CC79A7", "#009E73")) 
p


tikz(here::here("outputs", "VDP", "paper-plots", "vary-mu-results.tex"),
     width = 1.5 * doc_width_inches, 
     height =  2.25 * doc_width_inches,
     standAlone = TRUE)
p
dev.off()


tinytex::lualatex(here::here("outputs", "VDP", "paper-plots", "vary-mu-results.tex"))

