# ------------------------------------------------------------------------#
# Simple simulation to demonstrate general proposed bias correction on
# SHM model.
# ------------------------------------------------------------------------#

# Load necessary functions and packages: ----------------------------------
source(here::here("code", "SHM-helper-functions.R"))
source(here::here("DHM-TV-general-data-generation.R"))
source(here::here("PDA-2-pw-general.R"))
library(fda) # Functional Data Analysis, CRAN v5.5.1

# Path to store results: --------------------------------------------------
results_path <- here::here("outputs", "SHM", "simulation-results", "general-DHM-simulation.rds")


# Settings for simulation: ------------------------------------------------
grid_range <- c(0, 4 * pi)
n_grid <- 200
grid_points_i <- seq(from = grid_range[1],
                     to = grid_range[2],
                     length.out = n_grid)
n_i <- 200
num_iter <- 10
set.seed(1996)

damping_fun_test <- function(t) {
  0.01 *  (t - 2 * pi) ^2
}


dataset_i <- generate_DHM_TV_dataset_general(n = n_i,
                                          grid_points = grid_points_i,
                                          sigma = 0.4,
                                          mu_init = c(1, 0),
                                          intensity = 1,
                                          damping_fun = damping_fun_test,
                                          sigma_init = diag(rep(0.05, 2)))



x_i_test <- dataset_i$x
Dx_i_test <- dataset_i$dx
eps_i_test <- dataset_i$noise
D2x_i_test <- - x_i_test + damping_fun_test(grid_points_i) * Dx_i_test + eps_i_test

par(mfrow = c(1, 3))
matplot(x_i_test, type = "l")
matplot(Dx_i_test, type = "l")
matplot(D2x_i_test, type = "l")
# plot(damping_fun_test(grid_points_i))


try(expr = {pda_test <- do_pda_iterative_br_parallel(grid_points = grid_points_i,
                                                                         x = x_i_test,
                                                                         Dx = Dx_i_test,
                                                                         D2x = D2x_i_test,
                                                                         n = n_i, 
                                                                         num_iter = num_iter,
                                                                         verbose = TRUE, 
                                                                         n_cores = 8)})



par(mfrow = c(1, 3))
matplot(pda_test$beta[,1,], type = "l")
matplot(pda_test$beta[,2,], type = "l")
matplot(pda_test$beta[,3,], type = "l")
lines(damping_fun_test(grid_points_i))




