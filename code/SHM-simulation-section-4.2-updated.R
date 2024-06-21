source(here::here("code", "SHM-bias-calculation-function.R"))
source(here::here("code", "SHM-helper-functions.R"))
library(data.table) # CRAN v1.14.2
library(foreach)    # CRAN v1.5.2
library(doParallel) # CRAN v1.0.17
library(abind)      # CRAN v1.4-5

grid_range <- c(0, 2 * pi)
n_grid <- 100
grid_points <- seq(from = grid_range[1],
                   to = grid_range[2],
                   length.out = n_grid)
n_sim <- 500
n_iter <- 3
results_seed <- vector(mode = "list", length = n_sim)
set.seed(1996)


cl <- makePSOCKcluster(8)
registerDoParallel(cl)
acomb <- function(...) abind(..., along=3)

system.time(
results_list <- foreach(i = seq_len(n_sim)#,
                         # .combine='acomb',
                         # .multicombine=TRUE
                         ) %dopar% {
  print(paste("Iteration", i))
  results_seed[[i]] <- .Random.seed
  dataset_i <- generate_SHM_dataset(n = 500,
                                    grid_points = grid_points,
                                    sigma = 0.25,
                                    mu_init = c(0, 0),
                                    intensity = 2,
                                    sigma_init = diag(rep(0.05, 2)))
  dataset_i$D2x <- - 1 * dataset_i$x + dataset_i$noise
  fit_SHM_pda_pointwise_bc_known_D2x_tm(grid_points = grid_points,
                                                 x = dataset_i$x,
                                                 D2x = dataset_i$D2x,
                                                 niter = n_iter)$beta
})




saveRDS(object = list(time_stamp = timestamp(), session_info = sessionInfo(), 
                      results_list = results_list, results_seed = results_seed),
        file = here::here("outputs", "SHM", "simulation-results", "simulation-section-4.2-updated.rds"))
