# Re-run replicates of simulation using "slower" functions 
# for the iterations that failed because of numerical issue. 
# in this setting, only 2 replciates failed!

# Load packages and required functions: -----------------------------------
source(here::here("code", "VDP-helper-functions.R"))
source(here::here("code", "VDP-Expectation-Functions-Optimised.R"))
source(here::here("code", "VDP-bias-correction-function.R"))
library(parallel)
library(fda) # CRAN v5.5.1
ncores <- detectCores()

# Read in simulation results: ---------------------------------------------
simulation_results <- readRDS(
  file = here::here("outputs", "VDP", 
                    "vdp-simulation-section-5.1-vary-sigma-forcing.rds"))

grid_points <- simulation_results$grid_points
results_list <- simulation_results$results_list
settings <- simulation_results$settings
Sigma_ic <- simulation_results$Sigma_ic
n <- 200 # hard coded because forgot to store this (checked against original script and is the same for all of these simulations)
mu_ic <- c(1.9922213675594, -0.910974076470711) # same
n_iter <- 10L
(null_result_inds <- which(sapply(results_list, function(x){is.null(x$beta_array)})))
# [1] 124 130


simulation_results_updated <- simulation_results

for(i in null_result_inds) {
  print('----------------------------------------')
  print(paste0("Simulation rep ", i))
  # we have stored the simulation seed, so this ok:
  .Random.seed <- simulation_results$simulation_seeds[[i]]
  mu <- settings[i, "mu"]
  
  if(mu == 0.5) {
    Sigma_ic_i <- Sigma_ic
    Sigma_ic_i[2,2] <- 2^2 * Sigma_ic_i[2,2]
  }
  if(mu == 1) {
    Sigma_ic_i <- Sigma_ic
  }
  if(mu == 2) {
    Sigma_ic_i <- Sigma_ic
    Sigma_ic_i[2,2] <- 0.5^2 * Sigma_ic_i[2,2]
  }
  
  dataset_i <- generate_vdp_curves_stochastic_ic_and_forcing(n = n,
                                                             mu = mu,
                                                             mu_ic = mu_ic,
                                                             grid_points = grid_points,
                                                             Sigma_ic = Sigma_ic_i,
                                                             sigma_forcing = settings$sigma_forcing[i],
                                                             intensity = settings$intensity_forcing[i])
  x_i <- dataset_i$data[,1,]
  y_i <- dataset_i$data[,2,]
  eps_x_i <- dataset_i$forcing[,1,]
  eps_y_i <- dataset_i$forcing[,2,]
  dx_i <- mu * (x_i - (1/3) * (x_i^3) - y_i) + eps_x_i
  dy_i <- (1 / mu) * x_i + eps_y_i
  
  # Just to double check theis is the replicate that threw an error, we go and recreate the error
  # by re-running the optimising code and saving the result/ error
  retry_optimised_code <- try(do_pda_vdp_multiple_iterations_opt(x = x_i, 
                                          y = y_i,
                                          dx = dx_i, 
                                          dy = dy_i,
                                          grid_points = grid_points,
                                          num_iter = n_iter,
                                          silent = FALSE,
                                          n_cores = ncores))
  stopifnot(class(retry_optimised_code) == "try-error")
  
  try_not_optimised_code <- try(do_pda_vdp_multiple_iterations(x = x_i, 
                                                              y = y_i,
                                                              dx = dx_i, 
                                                              dy = dy_i,
                                                              grid_points = grid_points,
                                                              num_iter = n_iter,
                                                              silent = FALSE,
                                                              n_cores = ncores))
  # Make sure it didn't throw error again
  if(class(try_not_optimised_code) != "try-error") {
    print("Succesful with non-optimised code")
  } else if(class(try_not_optimised_code) == "try-error") {
    print("Failed again with non-otimised code")
  }
  
  simulation_results_updated$results_list[[i]] <- try_not_optimised_code
}


saveRDS(object = simulation_results_updated,
        file = here::here(
          "outputs", "VDP",
          "vdp-simulation-section-5.1-vary-sigma-forcing-updated.rds")


