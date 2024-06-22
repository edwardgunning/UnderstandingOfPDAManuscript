# This is a classic non-linear system: see its wikipedia papge, where we might take the
# two-dimensional form and force it in the same way that we did Simple Harmonic Motion.
# We’d need to tune the the length of the covariance kernel and σ along similar lines – have
# correlation die out at a fraction of the period, and a signal to noise ratio of 10 to 20. I’d
# add noise to both state variables, but wouldn’t worry about correlating the nois

source(here::here("code", "SHM-helper-functions.R"))

# Van der Pol Simulation

generate_vdp_deterministic <- function(mu, grid_points, y0, x0) {
  
  dynamics_equations <- function(t, y, ...) {
    with(as.list(c(y)),{
      # rate of change
      dX <-  mu * (X - (1/3)*(X^3) - Y)
      dY <- (1 / mu) * X
      list(c(dX, dY))
    }
    )
  }
  
  ic <- c(X = x0, Y = y0)
  deSolve::lsoda(y = ic,
                 times = grid_points,
                 func = dynamics_equations,
                 tcrit = max(grid_points))
  
}


generate_vdp_stochastic_ic <- function(mu, grid_points, mu_ic, Sigma_ic) {
  
  ic <- as.vector(mvtnorm::rmvnorm(n = 1, mean = mu_ic, sigma = Sigma_ic))
  names(ic) <- c("X", "Y")
  
  dynamics_equations <- function(t, y, ...) {
    with(as.list(c(y)),{
      # rate of change
      dX <-  mu * (X - (1/3)*(X^3) - Y)
      dY <- (1 / mu) * X
      list(c(dX, dY))
    }
    )
  }
  
  list(data = deSolve::lsoda(y = ic,
                             times = grid_points,
                             func = dynamics_equations,
                             tcrit = max(grid_points)),
       ic = ic)
  
}



generate_vdp_curves_stochastic_ic <- function(n, mu, grid_points, mu_ic, Sigma_ic) {
  data_array <- array(NA, dim = c(length(grid_points), 2, n))
  ic_mat <- matrix(NA, nrow = n, ncol = 2)
  
  for(i in seq_len(n)) {
    result_i <- generate_vdp_stochastic_ic(mu = mu, grid_points = grid_points, mu_ic = mu_ic, Sigma_ic = Sigma_ic)
    data_array[,,i] <- result_i$data[, 2:3]
    ic_mat[i, ] <- result_i$ic
  }
  
  list(data = data_array, 
       ic = ic_mat, 
       grid_points = grid_points)
}







# -------------------------------------------------------------------------

generate_vdp_stochastic_ic_and_forcing <- function(mu, grid_points, mu_ic, Sigma_ic, sigma_forcing, intensity) {
  
  ic <- as.vector(mvtnorm::rmvnorm(n = 1, mean = mu_ic, sigma = Sigma_ic))
  names(ic) <- c("X", "Y")
  
  smooth_noise <- generate_random_noise_curves(n = 2, grid_points = grid_points, sigma = sigma_forcing, intensity = intensity)
  smooth_noise_x_fun <- approxfun(x = grid_points, y = smooth_noise[1,])
  smooth_noise_y_fun <- approxfun(x = grid_points, y = smooth_noise[2,])
  
  dynamics_equations <- function(t, y, ...) {
    with(as.list(c(y)),{
      # rate of change
      dX <-  mu * (X - (1/3)*(X^3) - Y) + smooth_noise_x_fun(t)
      dY <- (1 / mu) * X + smooth_noise_y_fun(t)
      list(c(dX, dY))
    }
    )
  }
  
  list(data = deSolve::lsoda(y = ic,
                             times = grid_points,
                             func = dynamics_equations,
                             tcrit = max(grid_points)),
       ic = ic,
       forcing = smooth_noise)
  
}

generate_vdp_curves_stochastic_ic_and_forcing <- function(n, mu, grid_points, mu_ic, Sigma_ic, sigma_forcing, intensity) {
  data_array <- forcing_array <- array(NA, dim = c(length(grid_points), 2, n))
  ic_mat <- matrix(NA, nrow = n, ncol = 2)
  
  for(i in seq_len(n)) {
    result_i <- generate_vdp_stochastic_ic_and_forcing(mu = mu,
                                                       grid_points = grid_points,
                                                       mu_ic = mu_ic, 
                                                       Sigma_ic = Sigma_ic, sigma_forcing = sigma_forcing, intensity = intensity)
    data_array[,,i] <- result_i$data[, 2:3]
    ic_mat[i, ] <- result_i$ic
    forcing_array[,, i] <- t(result_i$forcing)
  }
  
  list(data = data_array, 
       ic = ic_mat, 
       forcing = forcing_array,
       grid_points = grid_points)
}
