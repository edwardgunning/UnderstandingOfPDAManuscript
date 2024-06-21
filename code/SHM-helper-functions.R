# Helper Functions for Simple Harmonic Motion Simulation: -----------------

generate_random_noise_curves <- function(grid_points, n, sigma, intensity = 2) {
  n <- as.integer(n)
  stopifnot(sigma > 0 & n > 0 & is.numeric(grid_points))
  
  # Matrix of outer product |s-t|
  s_minus_t_mat <- outer(X = grid_points, Y = grid_points, FUN = "-")
  abs_s_minus_t_mat <- abs(s_minus_t_mat)
  
  # $K(|t-s|) = \sigma^2\phi(2*|t-s|)$:
  K <- (sigma ^ 2) * dnorm(x = intensity * abs_s_minus_t_mat)
  
  # Random draw of n curves at gridpoints
  mvtnorm::rmvnorm(n = n, sigma = K)
}


# Test out:
# set.seed(1996)
# matplot(t(generate_random_noise_curves(grid_points = seq(0, 2 * pi, length.out = 100),
#                                        n = 100, sigma = 0.1)),type = "l", ylab = "")

# Further: Maybe add option to set random seed within.
# Also add extra tests
testthat::expect_error(generate_random_noise_curves(grid_points = seq(0, 100), n = 5, sigma = - 1))
testthat::expect_error(generate_random_noise_curves(grid_points = seq(0, 100), n = 0.4, sigma = 0.1))




# -------------------------------------------------------------------------


generate_SHM <- function(grid_points, sigma, mu_init, sigma_init, intensity = 2) {
  
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
      dY <-  - 1 * X + noise_as_fun(t)
      list(c(dX, dY))
    }
    )
  }
  
  x <- deSolve::lsoda(y = init_draw,
                      times = grid_points,
                      func = dynamics_equations,
                      tcrit = max(grid_points))[,2]
  
  list(x = x, noise = c(noise), init = init_draw)
}


generate_SHM_dataset <- function(n, grid_points, sigma, mu_init, sigma_init, intensity = 2) {
  
  x_mat <- noise_mat <- matrix(NA, ncol = n, nrow = length(grid_points))
  init_mat <- matrix(NA, ncol = 2, nrow = n)
  
  for(i in seq_len(n)) {
    shm <- generate_SHM(grid_points = grid_points,
                        sigma = sigma,
                        intensity = intensity,
                        mu_init = mu_init, 
                        sigma_init = sigma_init)
    x_mat[, i] <- shm$x
    noise_mat[, i] <- shm$noise
    init_mat[i, ] <- shm$init
  }
  
  list(x = x_mat, noise = noise_mat, init = init_mat)
  
}




calculate_deriv <- function(grid_points, x, norder) {
  
  stopifnot(norder %in% c(1, 2))
  
  order_6_bspline <- fda::create.bspline.basis(rangeval = range(grid_points),
                                               norder = 6,
                                               breaks = grid_points)
  x_fd <- fda::Data2fd(argvals = grid_points, y = x, basisobj = order_6_bspline)
  
  fda::eval.fd(evalarg = grid_points, fdobj = x_fd, Lfdobj = fda::int2Lfd(norder))
}


fit_SHM_pda_pointwise <- function(grid_points, x) {
  
  D2x <- calculate_deriv(grid_points = grid_points, x = x, norder = 2)
  
  beta <- vector(mode = "numeric", length = length(grid_points))
  noise_resid <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
  
  for(t in seq_along(grid_points)) {
    pda_t <- lm(D2x[t, ] ~ 0 + x[t, ])
    beta[t] <-  coef(pda_t)
    noise_resid[t, ] <- c(resid(pda_t))
  }
  
  list(beta = beta, noise_resid = noise_resid, x = x, D2x = D2x)
  
}




# PDA With bias correction ------------------------------------------------

fit_SHM_pda_pointwise_bc <- function(grid_points, x, niter = 1) {
  
  pda_init <- fit_SHM_pda_pointwise(grid_points = grid_points, x = x)
  D2x <- pda_init$D2x
  resid <- pda_init$noise_resid
  bias_denom <- apply(X = x ^2, 1, mean)
  
  beta_array <- matrix(data = NA, nrow = length(grid_points), ncol = (niter + 1))
  beta_array[, 1] <- pda_init$beta
  
  
  for(iter in seq_len(niter)) {
    # Turn latest beta into fda object and then function:
    beta_i <- beta_array[, iter]
    beta_u_fd <- fda::Data2fd(argvals = grid_points, y = beta_i)
    beta_u_est <- function(u) {
      fda::eval.fd(evalarg = u, fdobj = beta_u_fd)
    }
    # Now get estimate of residuals:
    resid <- D2x - sweep(x = x, MARGIN = 1, STATS = beta_i, FUN = "*")
    if(iter == 1) { # Sense check we calculate residuals correctly
      stopifnot(abs(resid - pda_init$noise_resid) < 10^-10)
    }
    
    # Set up residual covariance as function
    resid_fd <-  fda::Data2fd(argvals = grid_points, y = resid)
    cov_fd <- fda::var.fd(fdobj1 = resid_fd, fdobj2 = resid_fd)
    K_st_est <- function(s, t) {
      fda::eval.bifd(sevalarg = s, tevalarg = t, bifd = cov_fd)
    }
    
    # Calculate bias:
    Exe_calc <- bias_calculation(beta_u_est = beta_u_est,
                                 K_st_est = K_st_est,
                                 grid_points = grid_points)
    bias_est <- Exe_calc[1,2,] / bias_denom
    
    # Subtract updated estimate of bias from initial beta to get estimate:
    beta_array[, iter + 1] <- pda_init$beta - bias_est
  }
  
  list(beta = beta_array, final_resid = resid, pda_init = pda_init)
  
}





# -------------------------------------------------------------------------
# ASSUME D2x Perfectly OBSERVED:

fit_SHM_pda_pointwise_known_D2x <- function(grid_points, x, D2x) {
  
  stopifnot(dim(x) == dim(D2x))
  
  beta <- vector(mode = "numeric", length = length(grid_points))
  noise_resid <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
  
  for(t in seq_along(grid_points)) {
    pda_t <- lm(D2x[t, ] ~ 0 + x[t, ])
    beta[t] <-  coef(pda_t)
    noise_resid[t, ] <- c(resid(pda_t))
  }
  
  list(beta = beta, noise_resid = noise_resid, x = x, D2x = D2x)
  
}



fit_SHM_pda_pointwise_bc_known_D2x <- function(grid_points, x, niter = 1, D2x) {
  
  pda_init <- fit_SHM_pda_pointwise_known_D2x(grid_points = grid_points, x = x, D2x = D2x)
  D2x <- pda_init$D2x
  resid <- pda_init$noise_resid
  bias_denom <- apply(X = x ^2, 1, mean)
  
  beta_array <- matrix(data = NA, nrow = length(grid_points), ncol = (niter + 1))
  beta_array[, 1] <- pda_init$beta
  
  
  for(iter in seq_len(niter)) {
    # Turn latest beta into fda object and then function:
    beta_i <- beta_array[, iter]
    beta_u_fd <- fda::Data2fd(argvals = grid_points, y = beta_i)
    beta_u_est <- function(u) {
      fda::eval.fd(evalarg = u, fdobj = beta_u_fd)
    }
    # Now get estimate of residuals:
    resid <- D2x - sweep(x = x, MARGIN = 1, STATS = beta_i, FUN = "*")
    if(iter == 1) { # Sense check we calculate residuals correctly
      stopifnot(abs(resid - pda_init$noise_resid) < 10^-10)
    }
    
    # Set up residual covariance as function
    resid_fd <-  fda::Data2fd(argvals = grid_points, y = resid)
    cov_fd <- fda::var.fd(fdobj1 = resid_fd, fdobj2 = resid_fd)
    K_st_est <- function(s, t) {
      fda::eval.bifd(sevalarg = s, tevalarg = t, bifd = cov_fd)
    }
    
    # Calculate bias:
    Exe_calc <- bias_calculation(beta_u_est = beta_u_est,
                                 K_st_est = K_st_est,
                                 grid_points = grid_points)
    bias_est <- Exe_calc[1,2,] / bias_denom
    
    # Subtract updated estimate of bias from initial beta to get estimate:
    beta_array[, iter + 1] <- pda_init$beta - bias_est
  }
  
  list(beta = beta_array, final_resid = resid, pda_init = pda_init)
  
}





# -------------------------------------------------------------------------
# Add updated version with transition matrix...

fit_SHM_pda_pointwise_bc_tm <- function(grid_points, x, niter = 1) {
  
  pda_init <- fit_SHM_pda_pointwise(grid_points = grid_points, x = x)
  D2x <- pda_init$D2x
  resid <- pda_init$noise_resid
  bias_denom <- apply(X = x ^2, 1, mean)
  
  beta_array <- matrix(data = NA, nrow = length(grid_points), ncol = (niter + 1))
  beta_array[, 1] <- pda_init$beta
  
  
  for(iter in seq_len(niter)) {
    # Turn latest beta into fda object and then function:
    beta_i <- beta_array[, iter]
    beta_u_fd <- fda::Data2fd(argvals = grid_points, y = beta_i)
    beta_u_est <- function(u) {
      fda::eval.fd(evalarg = u, fdobj = beta_u_fd)
    }
    # Now get estimate of residuals:
    resid <- D2x - sweep(x = x, MARGIN = 1, STATS = beta_i, FUN = "*")
    if(iter == 1) { # Sense check we calculate residuals correctly
      stopifnot(abs(resid - pda_init$noise_resid) < 10^-10)
    }
    
    # Set up residual covariance as function
    resid_fd <-  fda::Data2fd(argvals = grid_points, y = resid)
    cov_fd <- fda::var.fd(fdobj1 = resid_fd, fdobj2 = resid_fd)
    K_st_est <- function(s, t) {
      fda::eval.bifd(sevalarg = s, tevalarg = t, bifd = cov_fd)
    }
    
    # Calculate bias:
    Exe_calc <- bias_calculation_with_transition_matrix(beta_u_est = beta_u_est,
                                 K_st_est = K_st_est,
                                 grid_points = grid_points)
    bias_est <- Exe_calc[1,2,] / bias_denom
    
    # Subtract updated estimate of bias from initial beta to get estimate:
    beta_array[, iter + 1] <- pda_init$beta - bias_est
  }
  
  list(beta = beta_array, final_resid = resid, pda_init = pda_init)
  
}


fit_SHM_pda_pointwise_bc_known_D2x_tm <- function(grid_points, x, niter = 1, D2x) {
  
  pda_init <- fit_SHM_pda_pointwise_known_D2x(grid_points = grid_points, x = x, D2x = D2x)
  D2x <- pda_init$D2x
  resid <- pda_init$noise_resid
  bias_denom <- apply(X = x ^2, 1, mean)
  
  beta_array <- matrix(data = NA, nrow = length(grid_points), ncol = (niter + 1))
  beta_array[, 1] <- pda_init$beta
  
  
  for(iter in seq_len(niter)) {
    # Turn latest beta into fda object and then function:
    beta_i <- beta_array[, iter]
    beta_u_fd <- fda::Data2fd(argvals = grid_points, y = beta_i)
    beta_u_est <- function(u) {
      fda::eval.fd(evalarg = u, fdobj = beta_u_fd)
    }
    # Now get estimate of residuals:
    resid <- D2x - sweep(x = x, MARGIN = 1, STATS = beta_i, FUN = "*")
    if(iter == 1) { # Sense check we calculate residuals correctly
      stopifnot(abs(resid - pda_init$noise_resid) < 10^-10)
    }
    
    # Set up residual covariance as function
    resid_fd <-  fda::Data2fd(argvals = grid_points, y = resid)
    cov_fd <- fda::var.fd(fdobj1 = resid_fd, fdobj2 = resid_fd)
    K_st_est <- function(s, t) {
      fda::eval.bifd(sevalarg = s, tevalarg = t, bifd = cov_fd)
    }
    
    # Calculate bias:
    Exe_calc <- bias_calculation_with_transition_matrix(beta_u_est = beta_u_est,
                                 K_st_est = K_st_est,
                                 grid_points = grid_points)
    bias_est <- Exe_calc[1,2,] / bias_denom
    
    # Subtract updated estimate of bias from initial beta to get estimate:
    beta_array[, iter + 1] <- pda_init$beta - bias_est
  }
  
  list(beta = beta_array, final_resid = resid, pda_init = pda_init)
  
}
