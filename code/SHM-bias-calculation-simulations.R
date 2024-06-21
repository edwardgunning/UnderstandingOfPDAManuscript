source(here::here("code", "SHM-bias-calculation-function.R"))
source(here::here("code", "SHM-helper-functions.R"))
library(data.table) # CRAN v1.14.2

grid_range <- c(0, 2 * pi)
n_grid <- 100
grid_points <- seq(from = grid_range[1],
                   to = grid_range[2],
                   length.out = n_grid)

n_sim <- 250
sigma_vals <- c(0.05, 0.2)
sigma_init_vals <- c(0.05, 0.2)
lengthscale_vals <- c(1, 2)
init_x0_ind_vals <- c(1, 2)
n_obs_vals <- c(100, 1000)
rep_vals <- seq_len(n_sim)

init_x0_mat <- matrix(data = c(
  1, 0,
  0, 1
), nrow = 2, ncol = 2, byrow = TRUE)

settings <- expand.grid(sigma = sigma_vals,
                        sigma_init = sigma_init_vals,
                        lengthscale = lengthscale_vals,
                        n_obs = n_obs_vals,
                        init_x0_ind = init_x0_ind_vals,
                        rep = rep_vals)
settings <- as.data.table(settings)

n_iter <- 3

unique_settings <- unique(settings[, -c("rep" , "n_obs")])

par(mfrow = c(4, 4))

for(i in seq_len(nrow(unique_settings))) {
  dataset_i <- generate_SHM_dataset(n = 50, 
                       grid_points = grid_points,
                       sigma = unique_settings[i, ]$sigma,
                       mu_init = init_x0_mat[unique_settings[i,]$init_x0, ],
                       intensity = unique_settings[i,]$lengthscale,
                       sigma_init = diag(rep(unique_settings[i, ]$sigma_init, 2)))
  
  matplot(dataset_i$x, type = "l")
  title(paste("settting", i))
}


# Simulation --------------------------------------------------------------
# settings <- as.data.table(settings)

results_array <- array(NA,
                       dim = c(length(grid_points), 
                               (n_iter + 1),
                               nrow(settings)))
results_seed <- vector(mode = "list", length = nrow(settings))
for(i in seq_len(nrow(settings))) { 
  if(i == 1) start_time <- proc.time()
  print(paste("Iteration", i, "of", nrow(settings)))
  print(paste("Time elapsed since start:", (proc.time() - start_time)["elapsed"]))
  
  
  results_seed[[i]] <- .Random.seed
  dataset_i <- generate_SHM_dataset(n = settings[i, ]$n_obs, 
                                    grid_points = grid_points,
                                    sigma = settings[i, ]$sigma,
                                    init_x0_mat[settings[i,]$init_x0, ], 
                                    intensity = settings[i,]$lengthscale,
                                    sigma_init = diag(rep(settings[i, ]$sigma_init, 2)))
  
  pda_i <- fit_SHM_pda_pointwise_bc(grid_points = grid_points, x = dataset_i$x, niter = n_iter)
  results_array[,, i] <- pda_i$beta
  }

saveRDS(object = list(array = results_array, seed = results_seed),
        file = here::here("outputs", "SHM", "simulation-results", "bias-correction-simulation.R"),
        )


