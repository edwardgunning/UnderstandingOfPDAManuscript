# Extend to 2 other PDA scenarios with different average initial conditions
source(here::here("code", "SHM-bias-calculation-function.R"))
source(here::here("code", "SHM-helper-functions.R"))
source(here::here("code", "SHM-covariance-calculation-function.R"))
library(data.table) # CRAN v1.14.2
library(fda)        # CRAN v5.5.1
library(parallel)

ncores <- detectCores()
grid_range <- c(0, 2 * pi)
# resudced grid size:
n_grid <- 50
grid_points <- seq(from = grid_range[1],
                   to = grid_range[2],
                   length.out = n_grid)

# Only compute upper-diagonal + diagonal elements because 
# Covariance function is symmetric.
test_grid <- rbind(
  cbind(grid_points, grid_points),
  t(combn(grid_points, 2))
  )

n_sim <- 10
n_iter <- 3
results_array <- array(NA,
                       dim = c(length(grid_points), 
                               (n_iter + 1),
                               n_sim))
results_seed <- vector(mode = "list", length = n_sim)
set.seed(1996)

final_resid_array <- array(NA, dim = c(length(grid_points), 500, n_sim))
covariance_est_list <- vector(mode = "list", length = n_sim)

for(i in seq_len(n_sim)) { 
  print(paste("Iteration", i, "of", n_sim))
  results_seed[[i]] <- .Random.seed
  dataset_i <- generate_SHM_dataset(n = 500, 
                                    grid_points = grid_points,
                                    sigma = 0.25,
                                    mu_init = c(0, 0), 
                                    intensity = 2,
                                    sigma_init = diag(rep(0.05, 2)))
  
  dataset_i$D2x <- - 1 * dataset_i$x + dataset_i$noise
  
  pda_i <- fit_SHM_pda_pointwise_bc_known_D2x_tm(grid_points = grid_points,
                                                 x = dataset_i$x, 
                                                 D2x = dataset_i$D2x,
                                                 niter = n_iter)
  results_array[,, i] <- pda_i$beta
  final_resid_mat <- pda_i$final_resid
  final_resid_array[,,i] <- final_resid_mat
  final_resid_fd <- Data2fd(argvals = grid_points, y = final_resid_mat)
  K_st_bifd <- var.fd(fdobj1 = final_resid_fd, fdobj2 = final_resid_fd)
  K_st_hat <- function(u, v) {
    eval.bifd(sevalarg = u, tevalarg = v, bifd = K_st_bifd)
  }
  beta_hat <- approxfun(x = grid_points, y = pda_i$beta[, n_iter+1])
    
  covariance_est_list[[i]] <- mclapply(seq_len(nrow(test_grid)), function (i) {
    calculate_covariance_shm_est(K_st_est = K_st_hat,
                                 beta_u_est = beta_hat,
                                  s = test_grid[i, 1],
                                  t = test_grid[i, 2])},
    mc.cores = ncores)
}

saveRDS(object = list(time_stamp = timestamp(),
                      results_array = results_array,
                      final_resid_mat= final_resid_mat,
                      covariance_est_list = covariance_est_list,
                      results_seed = results_seed,
                      session_info = sessionInfo()),
        file = here::here("outputs", "SHM",
                          "simulation-results",
                          "simulation-covariance-estimation.rds"))
