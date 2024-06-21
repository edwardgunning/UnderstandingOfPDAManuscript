source(here::here("code", "SHM-covariance-calculation-function.R"))
source(here::here("code", "SHM-helper-functions.R"))
library(parallel)
K_st_true <- function(s, t) {
  0.25 ^ 2 * dnorm(x = 2 * abs(t - s))
}

grid_range <- c(0, 2 * pi)
n_grid <- 100
grid_points <- seq(from = grid_range[1],
                   to = grid_range[2],
                   length.out = n_grid)

test_grid <- expand.grid(grid_points, grid_points)

ncores <- detectCores()
cov_xs_xt <- mclapply(seq_len(nrow(test_grid)), function (i) {
  calculate_covariance_shm_true(K_st_est = K_st_true,
                           s = test_grid[i, 1],
                           t = test_grid[i, 2])},
  mc.cores = ncores)


saveRDS(object = cov_xs_xt, 
        file = here::here("outputs", "SHM", "simulation-results", "cov_xs_xt_new.rds"))

