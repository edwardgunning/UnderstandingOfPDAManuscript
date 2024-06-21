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
                    "vdp-simulation-section-5.1-vary-intensity-forcing.rds"))
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


sum(sapply(results_list, function(x){is.null(x$beta_array)})) # Need tp check these out:


# Sigma = 0.1: ------------------------------------------------------------
inds_l_1 <- which(settings$intensity_forcing == 1)
results_l_1 <- results_list[inds_l_1]
estimates_l_1 <- map_dfr(.x = results_l_1, function(x) {
  
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

inds_l_2 <- which(settings$intensity_forcing == 2)
results_l_2 <- results_list[inds_l_2]
estimates_l_2 <- map_dfr(.x = results_l_2, function(x) {
  
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


inds_l_3 <- which(settings$intensity_forcing == 3)
results_l_3 <- results_list[inds_l_3]

estimates_l_3 <- map_dfr(.x = results_l_3, function(x) {
  
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

estimates_l_1$l <- 1
estimates_l_2$l <- 2
estimates_l_3$l <- 3

estimates_dt <- data.table(rbind(estimates_l_1,
                                 estimates_l_2,
                                 estimates_l_3))
estimates_dt_lng <- melt.data.table(estimates_dt,
                                    id.vars = c("t", "l", "estimator", "sim_rep"),
                                    variable.name = "parameter",
                                    value.name = "value",
                                    variable.factor = FALSE,
                                    value.factor = FALSE,
                                    verbose = TRUE)


estimates_dt_lng[, l := factor(l,
                               levels = c(1:3),
                               labels = paste0("$l = ", c(1:3), "$"))]
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
  facet_grid(parameter~l, scales = "free_y") +
  geom_line(aes( group = interaction(estimator, sim_rep), color = estimator), alpha = 0.4) +
  geom_line(data = truth_dt_lng, linewidth = 0.8) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.y = element_blank()) +
  labs(x = "$t$") +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 0.75))) +
  scale_color_manual(values = c("#CC79A7", "#009E73")) 
p

tikz(here::here("outputs", "VDP", "paper-plots", "vary-l-results.tex"),
     width = 1.5 * doc_width_inches, 
     height =  2.25 * doc_width_inches,
     standAlone = TRUE)
p
dev.off()

tinytex::lualatex(here::here("outputs", "VDP", "paper-plots", "vary-l-results.tex"))




