generate_empirical_random_noise_curve <- function(C_grid) {
  mvtnorm::rmvnorm(n = 1, sigma = C_grid)[1,]
}

generate_com <-  function(b0_fun,
                          beta_0_fun,
                          beta_1_fun,
                          mu_init_vec,
                          sigma_init_mat,
                          C_grid,
                          grid_points) {
  
  noise <- generate_empirical_random_noise_curve(C_grid = C_grid)
  noise_as_fun <- approxfun(x = grid_points,
                            y = noise,
                            method = "linear",
                            yleft = NA,
                            yright = NA)
  
  init_draw <- c(mvtnorm::rmvnorm(mean = mu_init_vec, sigma = sigma_init_mat, n = 1))
  names(init_draw) <- c("X", "Y")
  
  dynamics_equations <- function(t, y, ...) {
    with(as.list(c(y)),{
      # rate of change
      dX <- Y
      dY <-  b0_fun(t) + beta_0_fun(t) * X + beta_1_fun(t) * Y + noise_as_fun(t)
      list(c(dX, dY))
    }
    )
  }
  
  sol <- deSolve::lsoda(y = init_draw,
                        times = grid_points,
                        func = dynamics_equations,
                        tcrit = max(grid_points))[,2:3]
  
  list(x = sol[, 1], dx = sol[, 2], noise = c(noise), init = init_draw)
}


generate_com_dataset_general <- function(n,
                                         b0_fun,
                                         beta_0_fun,
                                         beta_1_fun,
                                         mu_init_vec,
                                         sigma_init_mat,
                                         C_grid,
                                         grid_points) {
  
  dx_mat <- x_mat <- noise_mat <- matrix(NA, ncol = n, nrow = length(grid_points))
  init_mat <- matrix(NA, ncol = 2, nrow = n)
  
  for(i in seq_len(n)) {
    COM <- generate_com(b0_fun = b0_fun,
                         beta_0_fun = beta_0_fun,
                         beta_1_fun = beta_1_fun,
                         mu_init_vec = mu_init_vec,
                         sigma_init_mat = sigma_init_mat,
                         C_grid = C_grid,
                         grid_points = grid_points)
    x_mat[, i] <- COM$x
    dx_mat[, i] <- COM$dx
    noise_mat[, i] <- COM$noise
    init_mat[i, ] <- COM$init
  }
  
  list(x = x_mat, dx = dx_mat, noise = noise_mat, init = init_mat)
}
