require(progress)
require(progress)
require(mgcv)
require(fda)
require(deSolve)
require(refund)


fit_pda_pointwise_general <- function(grid_points, x, Dx, D2x) {
  stopifnot(dim(x) == dim(D2x))
  stopifnot(dim(x) == dim(Dx))
  stopifnot(length(grid_points) == nrow(x))
  beta <- matrix(NA, nrow = nrow(x), ncol = 3)
  noise_resid <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
  
  for(t in seq_along(grid_points)) {
    pda_t <- lm(D2x[t, ] ~ x[t, ] + Dx[t, ])
    beta[t, ] <-  coef(pda_t)
    noise_resid[t, ] <- c(resid(pda_t))
  }
  list(beta = beta, noise_resid = noise_resid, x = x, Dx = Dx, D2x = D2x)
}



compute_pda_bias <- function(grid_points, beta_0_grid, beta_1_grid, resid_grid, n, x, Dx) {
  # Some routine checks to try to catch bugs.
  stopifnot(length(grid_points) == length(beta_0_grid))
  stopifnot(length(grid_points) == length(beta_1_grid))
  stopifnot(nrow(resid_grid) == length(grid_points))
  stopifnot(dim(x) == dim(resid_grid))
  stopifnot(dim(x) == dim(Dx))
  stopifnot(is.numeric(n) & (length(n)==1))
  stopifnot(ncol(resid_grid) == n)

  # Step 1: Convert grid-valued betas to functions:
  beta_0_fun <- approxfun(grid_points, beta_0_grid)
  beta_1_fun <- approxfun(grid_points, beta_1_grid)

  # Step 2: Create state transition matrix function using beta functions:
  transition_matrix_sarg_to_targ <- function(sarg, targ) {

    dynamics_equations <- function(t, y, ...) {
      with(as.list(c(y)),{
        # rate of change
        dX <- Y
        dY <-  beta_0_fun(t) * X + beta_1_fun(t) * Y
        list(c(dX, dY))
      }
      )
    }

    col1 <- deSolve::lsoda(y = c(X = 1, Y = 0),
                           times = c(sarg, targ),
                           func = dynamics_equations,
                           tcrit = targ)#only take 2nd row not interested in initial value

    col2 <- deSolve::lsoda(y = c(X = 0, Y = 1),
                           times = c(sarg, targ),
                           func = dynamics_equations,
                           tcrit = targ)#only take 2nd row not interested in initial value

    cbind(col1[2,-1], col2[2,-1])

  }

  # Step 3: Compute empirical estimate of residual covariance and convert to function:
  resid_fd <- Data2fd(argvals = grid_points, y = resid_grid)
  resid_cov_bifd <- var.fd(fdobj1 = resid_fd, fdobj2 = resid_fd)
  K_mat_fun <- function(s, t) {
    matrix(c(0, 0,
             0, eval.bifd(sevalarg = s, tevalarg = t, bifd = resid_cov_bifd)),
           byrow = TRUE, nrow = 2, ncol = 2
    )
  }

  # Step 4: Create integrand function (which is itself a matrix):
  integrand_st <- function(s, t) {
    transition_matrix_sarg_to_targ(sarg = s, targ = t) %*% K_mat_fun(s = s, t = t)
  }

  # Step 5: Split integrand matrix into entries so we can integrate separ --------
  # THESE MUST BE VECTORISED IN S SO THAT THEY CAN BE INTEGRATED ALONG S
  integrand_st_11 <- function(s, t) {
    return_vector <- vector(mode = "numeric", length = length(s))
    for(i in seq_along(s)) {
      return_vector[i] <- as.matrix(integrand_st(s = s[i], t = t))[1, 1]
    }
    return_vector
  }

  integrand_st_12 <- function(s, t) {
    return_vector <- vector(mode = "numeric", length = length(s))
    for(i in seq_along(s)) {
      return_vector[i] <- as.matrix(integrand_st(s = s[i], t = t))[1, 2]
    }
    return_vector
  }

  integrand_st_21 <- function(s, t) {
    return_vector <- vector(mode = "numeric", length = length(s))
    for(i in seq_along(s)) {
      return_vector[i] <- as.matrix(integrand_st(s = s[i], t = t))[2, 1]
    }
    return_vector
  }

  integrand_st_22 <- function(s, t) {
    return_vector <- vector(mode = "numeric", length = length(s))
    for(i in seq_along(s)) {
      return_vector[i] <- as.matrix(integrand_st(s = s[i], t = t))[2, 2]
    }
    return_vector
  }


  # Step 6: Create function to integrate the entries:
  integral_t_11 <- function(t) {
    integrate(f = integrand_st_11, lower = 0, upper = t, t = t)$value
  }

  integral_t_12 <- function(t) {
    integrate(f = integrand_st_12, lower = 0, upper = t, t = t)$value
  }

  integral_t_21 <- function(t) {
    integrate(f = integrand_st_21, lower = 0, upper = t, t = t)$value
  }

  integral_t_22 <- function(t) {
    integrate(f = integrand_st_22, lower = 0, upper = t, t = t)$value
  }

  # Step 7: Create function that evaluates makes matrix function out of 4 functions:
  integral_matrix <- function(t) {
    matrix(data = c(
      integral_t_11(t = t), integral_t_12(t = t),
      integral_t_21(t = t), integral_t_22(t = t)
    ), nrow = 2, ncol = 2, byrow = TRUE)
  }

  # Step 8: Compute Bias
  print("Computing bias:")
  bias <- matrix(data = NA, nrow = length(grid_points), ncol = 3)
  pb <- progress_bar$new(total = length(grid_points))
  for(t in seq_along(grid_points)) {
    pb$tick()
    Z <- cbind(1, x[t, ], Dx[t, ])
    ZtZinv <- solve(t(Z) %*% Z)
    Eexp_Z_eps <- c(0, integral_matrix(grid_points[t])[, 2])
    bias[t, ] <- ZtZinv %*% (n * Eexp_Z_eps)
  }
  bias
}



compute_pda_bias_parallel <- function(grid_points, beta_0_grid, beta_1_grid, resid_grid, n, x, Dx, n_cores) {
  # Some routine checks to try to catch bugs.
  stopifnot(length(grid_points) == length(beta_0_grid))
  stopifnot(length(grid_points) == length(beta_1_grid))
  stopifnot(nrow(resid_grid) == length(grid_points))
  stopifnot(dim(x) == dim(resid_grid))
  stopifnot(dim(x) == dim(Dx))
  stopifnot(is.numeric(n) & (length(n)==1))
  stopifnot(ncol(resid_grid) == n)
  
  # Step 1: Convert grid-valued betas to functions:
  beta_0_fun <- approxfun(grid_points, beta_0_grid)
  beta_1_fun <- approxfun(grid_points, beta_1_grid)
  
  # Step 2: Create state transition matrix function using beta functions:
  transition_matrix_sarg_to_targ <- function(sarg, targ) {
    
    dynamics_equations <- function(t, y, ...) {
      with(as.list(c(y)),{
        # rate of change
        dX <- Y
        dY <-  beta_0_fun(t) * X + beta_1_fun(t) * Y
        list(c(dX, dY))
      }
      )
    }
    
    col1 <- deSolve::lsoda(y = c(X = 1, Y = 0),
                           times = c(sarg, targ),
                           func = dynamics_equations,
                           tcrit = targ)#only take 2nd row not interested in initial value
    
    col2 <- deSolve::lsoda(y = c(X = 0, Y = 1),
                           times = c(sarg, targ),
                           func = dynamics_equations,
                           tcrit = targ)#only take 2nd row not interested in initial value
    
    cbind(col1[2,-1], col2[2,-1])
    
  }
  
  # Step 3: Compute empirical estimate of residual covariance and convert to function:
  resid_fd <- Data2fd(argvals = grid_points, y = resid_grid)
  resid_cov_bifd <- var.fd(fdobj1 = resid_fd, fdobj2 = resid_fd)
  K_mat_fun <- function(s, t) {
    matrix(c(0, 0,
             0, eval.bifd(sevalarg = s, tevalarg = t, bifd = resid_cov_bifd)),
           byrow = TRUE, nrow = 2, ncol = 2
    )
  }
  
  # Step 4: Create integrand function (which is itself a matrix):
  integrand_st <- function(s, t) {
    transition_matrix_sarg_to_targ(sarg = s, targ = t) %*% K_mat_fun(s = s, t = t)
  }
  
  # Step 5: Split integrand matrix into entries so we can integrate separ --------
  # THESE MUST BE VECTORISED IN S SO THAT THEY CAN BE INTEGRATED ALONG S
  integrand_st_11 <- function(s, t) {
    return_vector <- vector(mode = "numeric", length = length(s))
    for(i in seq_along(s)) {
      return_vector[i] <- as.matrix(integrand_st(s = s[i], t = t))[1, 1]
    }
    return_vector
  }
  
  integrand_st_12 <- function(s, t) {
    return_vector <- vector(mode = "numeric", length = length(s))
    for(i in seq_along(s)) {
      return_vector[i] <- as.matrix(integrand_st(s = s[i], t = t))[1, 2]
    }
    return_vector
  }
  
  integrand_st_21 <- function(s, t) {
    return_vector <- vector(mode = "numeric", length = length(s))
    for(i in seq_along(s)) {
      return_vector[i] <- as.matrix(integrand_st(s = s[i], t = t))[2, 1]
    }
    return_vector
  }
  
  integrand_st_22 <- function(s, t) {
    return_vector <- vector(mode = "numeric", length = length(s))
    for(i in seq_along(s)) {
      return_vector[i] <- as.matrix(integrand_st(s = s[i], t = t))[2, 2]
    }
    return_vector
  }
  
  
  # Step 6: Create function to integrate the entries:
  integral_t_11 <- function(t) {
    integrate(f = integrand_st_11, lower = 0, upper = t, t = t)$value
  }
  
  integral_t_12 <- function(t) {
    integrate(f = integrand_st_12, lower = 0, upper = t, t = t)$value
  }
  
  integral_t_21 <- function(t) {
    integrate(f = integrand_st_21, lower = 0, upper = t, t = t)$value
  }
  
  integral_t_22 <- function(t) {
    integrate(f = integrand_st_22, lower = 0, upper = t, t = t)$value
  }
  
  # Step 7: Create function that evaluates makes matrix function out of 4 functions:
  integral_matrix <- function(t) {
    matrix(data = c(
      integral_t_11(t = t), integral_t_12(t = t),
      integral_t_21(t = t), integral_t_22(t = t)
    ), nrow = 2, ncol = 2, byrow = TRUE)
  }
  
  # Step 8: Compute Bias
  print("Computing bias:")
  bias <- matrix(data = NA, nrow = length(grid_points), ncol = 3)
  pb <- progress_bar$new(total = length(grid_points))
  
  bias_list <- parallel::mclapply(
    X = seq_along(grid_points),
    FUN = function(t) {
      Z <- cbind(1, x[t, ], Dx[t, ])
      ZtZinv <- solve(t(Z) %*% Z)
      Eexp_Z_eps <- c(0, integral_matrix(grid_points[t])[, 2])
      ZtZinv %*% (n * Eexp_Z_eps)
    }, mc.cores = n_cores)
  
  for(t in seq_along(grid_points)) {
    bias[t, ] <- bias_list[[t]]
    }
  
  bias
}



do_pda_iterative_br <- function(grid_points, x, Dx, D2x, n, num_iterations, verbose = TRUE) {
  # Some routine checks:
  stopifnot(is.numeric(n) & (length(n)==1))
  stopifnot(is.numeric(grid_points))
  grid_length <- length(grid_points)
  stopifnot(dim(x) == c(grid_length, n))
  stopifnot(dim(x) == dim(Dx))
  stopifnot(dim(x) == dim(D2x))
  num_iterations <- as.integer(num_iterations)
  if(!(num_iterations>=1)) stop("num_iterations must be an integer greater than or equal to 1.")
  
  # Set up objects to store results:
  beta_array <- array(NA, c(grid_length, 3, num_iterations + 1))
  bias_array <- array(NA, c(grid_length, 3, num_iterations))
  
  # Do initial PDA (pointwise by OLS) fit.
  pda_initial_fit <- fit_pda_pointwise_general(grid_points = grid_points,
                                               x = x,
                                               Dx = Dx,
                                               D2x = D2x)
  stopifnot(dim(beta_array[,,1]) == dim(pda_initial_fit$beta))
  
  # Initialise beta estimate at initial OLS estimate:
  beta_array[,,1] <- pda_initial_fit$beta
  
  for(it in seq_len(num_iterations)) {
    if(verbose) {print(paste("Iteration", it, "of", num_iterations))}
    # Update Betas Being Used in Bias Correction:
    beta_current <- beta_array[,, it]
    beta_0_current_grid <- beta_current[, 2]
    beta_1_current_grid <- beta_current[, 3]
    
    # Update Residuals Being Used
    if(it == 1) { # if first iteration just use residuals from pda.
      resid_current <- pda_initial_fit$noise_resid
    } else { # otherwise calculate new residuals based on new beta parameters.
      for(t in seq_along(grid_points)) {
        Z <- cbind(1, x[t, ], Dx[t, ])
        resid_current[t, ] <- D2x[t, ] - Z %*% beta_current[t, ]
      }
    }
    # Compute bias using current values 
    bias_current <- compute_pda_bias(grid_points = grid_points,
                                     n = n,
                                     beta_0_grid = beta_0_current_grid,
                                     beta_1_grid = beta_1_current_grid,
                                     x = x,
                                     Dx = Dx,
                                     resid_grid = resid_current)
    # Stores bias:
    bias_array[,,it] <- bias_current
    # Update beta: Initial OLS estimate - current bias estimate
    beta_array[,,it + 1] <- beta_array[,,1] - bias_current
  }
  
  # Return parameter and bias estimates:
  list(
    beta = beta_array,
    bias = bias_array,
    resid = list(
      initial = pda_initial_fit$noise_resid,
      final = resid_current
    )
  )
  
}



# Parallell Version -------------------------------------------------------

do_pda_iterative_br_parallel <- function(grid_points, x, Dx, D2x, n, num_iterations, verbose = TRUE, n_cores) {
  # Some routine checks:
  stopifnot(is.numeric(n) & (length(n)==1))
  stopifnot(is.numeric(grid_points))
  grid_length <- length(grid_points)
  stopifnot(dim(x) == c(grid_length, n))
  stopifnot(dim(x) == dim(Dx))
  stopifnot(dim(x) == dim(D2x))
  num_iterations <- as.integer(num_iterations)
  if(!(num_iterations>=1)) stop("num_iterations must be an integer greater than or equal to 1.")
  
  # Set up objects to store results:
  beta_array <- array(NA, c(grid_length, 3, num_iterations + 1))
  bias_array <- array(NA, c(grid_length, 3, num_iterations))
  
  # Do initial PDA (pointwise by OLS) fit.
  pda_initial_fit <- fit_pda_pointwise_general(grid_points = grid_points,
                                               x = x,
                                               Dx = Dx,
                                               D2x = D2x)
  stopifnot(dim(beta_array[,,1]) == dim(pda_initial_fit$beta))
  
  # Initialise beta estimate at initial OLS estimate:
  beta_array[,,1] <- pda_initial_fit$beta
  
  for(it in seq_len(num_iterations)) {
    if(verbose) {print(paste("Iteration", it, "of", num_iterations))}
    # Update Betas Being Used in Bias Correction:
    beta_current <- beta_array[,, it]
    beta_0_current_grid <- beta_current[, 2]
    beta_1_current_grid <- beta_current[, 3]
    
    # Update Residuals Being Used
    if(it == 1) { # if first iteration just use residuals from pda.
      resid_current <- pda_initial_fit$noise_resid
    } else { # otherwise calculate new residuals based on new beta parameters.
      for(t in seq_along(grid_points)) {
        Z <- cbind(1, x[t, ], Dx[t, ])
        resid_current[t, ] <- D2x[t, ] - Z %*% beta_current[t, ]
      }
    }
    # Compute bias using current values 
    bias_current <- compute_pda_bias_parallel(grid_points = grid_points,
                                     n = n,
                                     beta_0_grid = beta_0_current_grid,
                                     beta_1_grid = beta_1_current_grid,
                                     x = x,
                                     Dx = Dx,
                                     resid_grid = resid_current,
                                     n_cores = n_cores)
    # Stores bias:
    bias_array[,,it] <- bias_current
    # Update beta: Initial OLS estimate - current bias estimate
    beta_array[,,it + 1] <- beta_array[,,1] - bias_current
  }
  
  # Return parameter and bias estimates:
  list(
    beta = beta_array,
    bias = bias_array,
    resid = list(
      initial = pda_initial_fit$noise_resid,
      final = resid_current
    )
  )
  
}




# -------------------------------------------------------------------------



# Parallel Version with smoothing of beta coefficients:
do_pda_iterative_br_parallel_post_smooth <- function(grid_points, x, Dx, D2x, n, num_iterations, verbose = TRUE, n_cores, k, method = "REML", resid_smooth = FALSE, resid_k = 35, resid_pve = 0.99) {
  # Some routine checks:
  if(resid_smooth) {
    stopifnot(!is.null(resid_k))
    stopifnot(!is.null(resid_pve))
  }
  stopifnot(is.numeric(n) & (length(n)==1))
  stopifnot(is.numeric(grid_points))
  grid_length <- length(grid_points)
  stopifnot(dim(x) == c(grid_length, n))
  stopifnot(dim(x) == dim(Dx))
  stopifnot(dim(x) == dim(D2x))
  num_iterations <- as.integer(num_iterations)
  if(!(num_iterations>=1)) stop("num_iterations must be an integer greater than or equal to 1.")
  
  # Set up objects to store results:
  beta_array <- array(NA, c(grid_length, 3, num_iterations + 1))
  bias_array <- array(NA, c(grid_length, 3, num_iterations))
  
  # Do initial PDA (pointwise by OLS) fit.
  pda_initial_fit <- fit_pda_pointwise_general(grid_points = grid_points,
                                               x = x,
                                               Dx = Dx,
                                               D2x = D2x)
  stopifnot(dim(beta_array[,,1]) == dim(pda_initial_fit$beta))
  
  # Initialise beta estimate at initial OLS estimate:
  beta_array[,,1] <- pda_initial_fit$beta
  
  for(it in seq_len(num_iterations)) {
    if(verbose) {print(paste("Iteration", it, "of", num_iterations))}
    # Update Betas Being Used in Bias Correction:
    print("Post-smoothing beta functions")
    
    # Final postsmooth of functions:
    
    
    beta_current <- beta_array[,, it]
    for(j in 1:3) {
      beta_current[, j] <- predict(gam(beta_current[, j] ~ s(grid_points, bs = "bs", k = k), method = method))
    }
    beta_0_current_grid <- beta_current[, 2]
    beta_1_current_grid <- beta_current[, 3]
    
    # Update Residuals Being Used
    if(it == 1) { # if first iteration just use residuals from pda.
      resid_current <- pda_initial_fit$noise_resid
    } else { # otherwise calculate new residuals based on new beta parameters.
      for(t in seq_along(grid_points)) {
        Z <- cbind(1, x[t, ], Dx[t, ])
        resid_current[t, ] <- D2x[t, ] - Z %*% beta_current[t, ]
      }
    }
    
    if(resid_smooth) {
      face_smooth_resid <- fpca.face(Y = t(resid_current), 
                                     pve = resid_pve, 
                                     knots = resid_k)
      resid_current <- t(face_smooth_resid$Yhat)
    }
    
    # Compute bias using current values 
    bias_current <- compute_pda_bias_parallel(grid_points = grid_points,
                                              n = n,
                                              beta_0_grid = beta_0_current_grid,
                                              beta_1_grid = beta_1_current_grid,
                                              x = x,
                                              Dx = Dx,
                                              resid_grid = resid_current,
                                              n_cores = n_cores)
    
    # Stores bias:
    bias_array[,,it] <- bias_current
    for(j in 1:3) {
      bias_current[, j] <- predict(gam(bias_current[, j] ~ s(grid_points, bs = "bs", k = k), method = method))
    }
    
    # Update beta: Initial OLS estimate - current bias estimate
    beta_array[,,it + 1] <- beta_array[,,1] - bias_current
  }
  
  print("Final Postsmooth")
  # Final postsmooth of functions:
  for(j in 1:3) {
    beta_array[,j,num_iterations + 1] <- predict(gam(beta_array[,j,num_iterations + 1] ~ s(grid_points, bs = "bs", k = k), method = method))
  }
  
  
  # Return parameter and bias estimates:
  list(
    beta = beta_array,
    bias = bias_array,
    resid = list(
      initial = pda_initial_fit$noise_resid,
      final = resid_current
    )
  )
  
}

