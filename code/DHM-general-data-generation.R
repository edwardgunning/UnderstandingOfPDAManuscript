# -------------------------------------------------------------------------
generate_DHM_general <- function(grid_points, sigma, mu_init, sigma_init, intensity = 2, damping_param) {
  stopifnot(is.numeric(damping_param))
  noise <- generate_random_noise_curves(grid_points = grid_points, n = 1, sigma = sigma, intensity = intensity)
  noise_as_fun <- approxfun(x = grid_points,
                            y = noise,
                            method = "linear",
                            yleft = NA,
                            yright = NA)
  
  init_draw <- c(mvtnorm::rmvnorm(mean = mu_init, sigma = sigma_init, n = 1))
  names(init_draw) <- c("X", "Y")
  
  dynamics_equations <- function(t, y, ...) {
    with(as.list(c(y)),{
      # rate of change
      dX <- Y
      dY <-  - 1 * X + damping_param * Y + noise_as_fun(t)
      list(c(dX, dY))
    }
    )
  }
  
  result <- deSolve::lsoda(y = init_draw,
                           times = grid_points,
                           func = dynamics_equations,
                           tcrit = max(grid_points))[,2:3]
  
  list(x = result[,1], dx = result[,2], noise = c(noise), init = init_draw)
}


generate_DHM_dataset_general <- function(n, grid_points, sigma, mu_init, sigma_init, intensity = 2, damping_param) {
  
  dx_mat <- x_mat <- noise_mat <- matrix(NA, ncol = n, nrow = length(grid_points))
  init_mat <- matrix(NA, ncol = 2, nrow = n)
  
  for(i in seq_len(n)) {
    DHM <- generate_DHM_general(grid_points = grid_points,
                                sigma = sigma,
                                intensity = intensity,
                                mu_init = mu_init, 
                                sigma_init = sigma_init,
                                damping_param = damping_param)
    x_mat[, i] <- DHM$x
    dx_mat[, i] <- DHM$dx
    noise_mat[, i] <- DHM$noise
    init_mat[i, ] <- DHM$init
  }
  
  list(x = x_mat, dx = dx_mat, noise = noise_mat, init = init_mat)
  
}
