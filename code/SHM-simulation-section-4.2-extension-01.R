# Extend to 2 other PDA scenarios with different average initial conditions
source(here::here("code", "SHM-bias-calculation-function.R"))
source(here::here("code", "SHM-helper-functions.R"))
library(data.table) # CRAN v1.14.2

grid_range <- c(0, 2 * pi)
n_grid <- 100
grid_points <- seq(from = grid_range[1],
                   to = grid_range[2],
                   length.out = n_grid)
n_sim <- 200
n_iter <- 3

ic_settings <- matrix(c(0, 1, 1, 0), ncol = 2, nrow = 2, byrow = 2)
ic_setting_levels <- c(1, 2)
rep <- seq_len(n_sim)
settings_dt <- data.table(expand.grid(rep = rep, ic_setting_level = ic_setting_levels))

results_array <- array(NA,
                       dim = c(length(grid_points), 
                               (n_iter + 1),
                               nrow(settings_dt)))
results_seed <- vector(mode = "list", length = nrow(settings_dt))
set.seed(1996)
for(i in seq_len(nrow(settings_dt))) { 
  print(paste("Iteration", i, "of", nrow(settings_dt)))
  results_seed[[i]] <- .Random.seed
  dataset_i <- generate_SHM_dataset(n = 500, 
                                    grid_points = grid_points,
                                    sigma = 0.25,
                                    mu_init = ic_settings[settings_dt[i, ]$ic_setting_level, ], 
                                    intensity = 2,
                                    sigma_init = diag(rep(0.05, 2)))
  
  dataset_i$D2x <- - 1 * dataset_i$x + dataset_i$noise
  
  pda_i <- fit_SHM_pda_pointwise_bc_known_D2x_tm(grid_points = grid_points,
                                              x = dataset_i$x, 
                                              D2x = dataset_i$D2x,
                                              niter = n_iter)
  results_array[,, i] <- pda_i$beta
}



saveRDS(object = list(time_stamp = timestamp(),
                      results_array = results_array,
                      settings_dt = settings_dt,
                      results_seed = results_seed,
                      session_info = sessionInfo()),
        file = here::here("outputs", "SHM",
                          "simulation-results",
                          "simulation-section-4.2-extension-ic.rds"))
