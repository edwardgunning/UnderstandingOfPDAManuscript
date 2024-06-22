# Load packages and required functions: -----------------------------------
source(here::here("code", "VDP-helper-functions.R"))
source(here::here("code", "VDP-Expectation-Functions-Optimised.R"))
library(parallel)
library(fda)
ncores <- detectCores()


# Set up settings: --------------------------------------------------------
Nrep <- 50
levels_sigma_forcing <- c(0.1, 0.2, 0.4)
levels_intensity_forcing <- c(1, 2, 3)
levels_mu <- c(0.5, 1, 2)


# Vary different parameters one at a time: --------------------------------
# settings_vary_sigma_forcing <- expand.grid(rep = seq_len(Nrep),
#                                            mu = 1,
#                                            intensity_forcing = 2,
#                                            sigma_forcing = levels_sigma_forcing)

settings_vary_intensity_forcing <- expand.grid(rep = seq_len(Nrep),
                                           mu = 1,
                                           intensity_forcing = levels_intensity_forcing,
                                           sigma_forcing = 0.1)
# settings_vary_mu <- expand.grid(rep = seq_len(Nrep),
#                                                mu = levels_mu,
#                                                intensity_forcing = 2,
#                                                sigma_forcing = 0.1)



settings <- settings_vary_intensity_forcing

# Parameters that will be fixed: ------------------------------------------
mu_ic <- c(1.9922213675594, -0.910974076470711)
Sigma_ic <- matrix(c(0.025, 0, 0, 0.025),
                   nrow = 2, 
                   ncol = 2) # Not too much IC noise
grid_points <- seq(0, 13, length.out = 200) 
n <- 200 # number of curves
results_list <- vector(mode = "list", length = nrow(settings))
n_iter <- 10L
simulation_seeds <- vector(mode = "list", length = nrow(settings))



# Loop through simulation replicates: -------------------------------------
for(i in 1:nrow(settings)) {
  print('----------------------------------------')
  print(paste0("Simulation rep ", i))
  simulation_seeds[[i]] <- .Random.seed
  mu <- settings[i, "mu"]
  # Do some manual adjustment of initial conditions based on  value of mu
  # to make datasets realistic.
  
  # ^ not needed here because \mu fixed = 1
  
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
  
  try(results_list[[i]] <- do_pda_vdp_multiple_iterations_opt(x = x_i, 
                                                              y = y_i,
                                                              dx = dx_i, 
                                                              dy = dy_i,
                                                              grid_points = grid_points,
                                                              num_iter = n_iter,
                                                              silent = FALSE,
                                                              n_cores = ncores))
}


# Save results: -----------------------------------------------------------
saveRDS(object = list(
  timestamp = timestamp(),
  R.version = R.version,
  simulation_seeds = simulation_seeds,
  sessionInfo = sessionInfo(),
  grid_points = grid_points,
  mu = mu,
  Sigma_ic = Sigma_ic,
  settings = settings,
  results_list = results_list), 
  file =  here::here("outputs", "VDP", "vdp-simulation-section-5.1-vary-intensity-forcing.rds"))
