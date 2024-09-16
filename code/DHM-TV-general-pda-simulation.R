# ------------------------------------------------------------------------#
# Simple simulation to demonstrate general proposed bias correction on
# SHM model.
# ------------------------------------------------------------------------#

# Load necessary functions and packages: ----------------------------------
source(here::here("code", "SHM-helper-functions.R"))
source(here::here("code", "DHM-TV-general-data-generation.R"))
source(here::here("code", "PDA-2-pw-general.R"))
library(fda) # Functional Data Analysis, CRAN v5.5.1

# Path to store results: --------------------------------------------------
results_path <- here::here("outputs",
                           "SHM", 
                           "simulation-results",
                           "general-DHM-TV-simulation.rds")


# Settings for simulation: ------------------------------------------------
grid_range <- c(0, 4 * pi)
n_grid <- 200
grid_points_i <- seq(from = grid_range[1],
                     to = grid_range[2],
                     length.out = n_grid)
n_i <- 200
num_simulations <- 50
simulation_seeds <- vector(mode = "list", length = num_simulations)
num_iter <- 10
beta_array_simulations <- array(NA, dim = c(length(grid_points_i),3, num_iter + 1, num_simulations))
set.seed(1996)
error_catch <- vector("list", length = num_simulations)
damping_fun_test <- function(t) {
  0.01 *  (t - 2 * pi) ^2
}
# Run simulations ---------------------------------------------------------
for(i in seq_len(num_simulations)) {
  print(paste("Simulation Number", i))
  simulation_seeds[[i]] <- .Random.seed
  
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
  D2x_i_test <- - x_i_test + damping_fun_test(grid_points_i)* Dx_i_test + eps_i_test
  
  error_catch[[i]] <- try(expr = {pda_test <- do_pda_iterative_br_parallel(grid_points = grid_points_i,
                                                                           x = x_i_test,
                                                                           Dx = Dx_i_test,
                                                                           D2x = D2x_i_test,
                                                                           n = n_i, 
                                                                           num_iter = num_iter,
                                                                           verbose = TRUE, 
                                                                           n_cores = 8)
  beta_array_simulations[,,,i] <- pda_test$beta})
}

# Check if there were any fails: ------------------------------------------
stopifnot(sapply(error_catch, class) == "array")





# Rough, quick and dirty plots of results: --------------------------------
par(mfrow = c(1, 3))
matplot(beta_array_simulations[,1,1,], type = "l", col = 1,
        ylim = range(beta_array_simulations[,1,,]))
# matlines(beta_array_simulations[,1,2,], col = 2)
# matlines(beta_array_simulations[,1,5,], col = 3)
matlines(beta_array_simulations[,1,11,], col = 4)
lines(apply(beta_array_simulations[,1,11,], 1, mean), lwd = 2)
abline(h = 0)

matplot(beta_array_simulations[,2,1,], type = "l", col = 1)
# matlines(beta_array_simulations[,2,2,], col = 2)
# matlines(beta_array_simulations[,2,5,], col = 3)
matlines(beta_array_simulations[,2,11,], col = 4)
lines(apply(beta_array_simulations[,2,11,], 1, mean), lwd = 2)

matplot(beta_array_simulations[,3,1,], type = "l", col = 1,
        ylim = range(beta_array_simulations[,3,,], na.rm = TRUE))
# matlines(beta_array_simulations[,3,2,], col = 2)
# matlines(beta_array_simulations[,3,5,], col = 3)
matlines(beta_array_simulations[,3,11,], col = 4)
lines(apply(beta_array_simulations[,3,11,], 1, mean), lwd = 2)

saveRDS(object = list(
  beta = beta_array_simulations,
  grid_points = grid_points_i,
  num_iter = num_iter,
  num_simulations = num_simulations,
  sigma = 0.25,
  mu_init = c(0, 2),
  intensity = 0.5,
  sigma_init = diag(rep(0.05, 2)),
  n_i = n_i,
  num_simulations = num_simulations,
  num_iter = num_iter,
  simulation_seeds = simulation_seeds
), file = results_path)


