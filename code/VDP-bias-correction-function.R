VDP_Expectation_function <- function(t, 
                                 beta_xx,
                                 beta_xy,
                                 beta_yx,
                                 beta_yy,
                                 K_st_mat_fun) {
  
  # -------------------------------------------------------------------------
  transition_matrix_sarg_to_targ <- function(sarg, targ) {
    
    dynamics_equations <- function(t, y, ...) {
      with(as.list(c(y)),{
        # rate of change
        dX <- beta_xx(t) * X + beta_xy(t) * Y
        dY <- beta_yx(t) * X + beta_yy(t) * Y
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
  
  
  
  integrand_matrix <- function(sevalarg, tevalarg) {
    transition_matrix_sarg_to_targ(sarg = sevalarg, targ = tevalarg) %*%
      K_st_mat_fun(sevalarg = sevalarg, tevalarg = tevalarg)
  }
  
  
  # -------------------------------------------------------------------------
  integrand_11 <- function(sevalarg, tevalarg) {
    integrand_matrix(sevalarg = sevalarg, tevalarg = tevalarg)[1,1]
  }
  
  integrand_12 <- function(sevalarg, tevalarg) {
    integrand_matrix(sevalarg = sevalarg, tevalarg = tevalarg)[1,2]
  }
  
  integrand_21 <- function(sevalarg, tevalarg) {
    integrand_matrix(sevalarg = sevalarg, tevalarg = tevalarg)[2,1]
  }
  
  integrand_22 <- function(sevalarg, tevalarg) {
    integrand_matrix(sevalarg = sevalarg, tevalarg = tevalarg)[2,2]
  }
  
  
  # Vectorise --------------------------------------------------------------
  integrand_11_vectorised <- Vectorize(integrand_11, vectorize.args = "sevalarg")
  integrand_12_vectorised <- Vectorize(integrand_12, vectorize.args = "sevalarg")
  integrand_21_vectorised <- Vectorize(integrand_21, vectorize.args = "sevalarg")
  integrand_22_vectorised <- Vectorize(integrand_22, vectorize.args = "sevalarg")
  
  
  # -------------------------------------------------------------------------
    matrix(data = c(
      integrate(f = integrand_11_vectorised, lower = 0, upper = t, tevalarg = t)$value,
      integrate(f = integrand_12_vectorised, lower = 0, upper = t, tevalarg = t)$value,
      integrate(f = integrand_21_vectorised, lower = 0, upper = t, tevalarg = t)$value,
      integrate(f = integrand_22_vectorised, lower = 0, upper = t, tevalarg = t)$value
    ), nrow = 2, ncol = 2, byrow = TRUE)
  }


do_pda_vdp_multiple_iterations <- function(x, y, dx, dy, grid_points, num_iter, silent = FALSE, n_cores = 8) {
  
  # Routine Checks on Input
  stopifnot(dim(x) == dim(y))
  stopifnot(dim(x) == dim(dx))
  stopifnot(dim(x) == dim(dy))
  stopifnot(nrow(x) == length(grid_points))
  stopifnot(is.integer(num_iter) & num_iter > 0)
  N <- ncol(x)
  D <- length(grid_points)
  
  # Set up storage for restults:
  beta_array <- bias_est_array <- array(NA, dim = c(D, 3, 2, num_iter + 1))
  expectation_est_array <- array(NA, dim = c(D, 2, 2, num_iter))
  resid_array <- array(NA, dim = c(D, N, 2, num_iter))
  ZtZ_inv_array  <- array(NA, dim = c(D, 3, 3))
  Z_array <- array(data = NA, c(D, N, 3))
  
  if(!silent) {
    print("Fitting Initial Pointwise Models")
  }
  for(tind in seq_along(grid_points)) {
    # Exy part
    Z_array[tind,,] <- Z <- cbind(1, x[tind,], y[tind,])
    ZtZ_inv_array[tind,,] <- ZtZ_inv <- solve(t(Z) %*% Z)
    
    # PDA/ regression part:
    lm_xtind <- lm(dx[tind, ] ~ x[tind, ] + y[tind, ])
    lm_ytind <- lm(dy[tind, ] ~ x[tind, ] + y[tind, ])
    
    beta_array[tind,,1,1] <- coef(lm_xtind)
    beta_array[tind,,2,1] <- coef(lm_ytind)
    
    resid_array[tind,,1,1] <- resid(lm_xtind)
    resid_array[tind,,2,1] <- resid(lm_ytind)
  }  
  
  # Start bias correction ---------------------------------------------------
  # Turn beta estimates into functions:
  beta_xx_i <- approxfun(x = grid_points, y = beta_array[,2,1,1])
  beta_xy_i <- approxfun(x = grid_points, y = beta_array[,3,1,1])
  beta_yx_i <- approxfun(x = grid_points, y = beta_array[,2,2,1])
  beta_yy_i <- approxfun(x = grid_points, y = beta_array[,3,2,1])
  
  resid_dx_fd <- Data2fd(argvals = grid_points,
                         y = resid_array[,,1,1])
  resid_dy_fd <- Data2fd(argvals = grid_points,
                         y = resid_array[,,2,1])
  
  K_xx_bifd <- var.fd(fdobj1 = resid_dx_fd, fdobj2 = resid_dx_fd)
  K_xx_i <- function(sevalarg, tevalarg) {
    eval.bifd(sevalarg = sevalarg, tevalarg = tevalarg, 
              bifd = K_xx_bifd)
  }
  k_yy_bifd <- var.fd(fdobj1 = resid_dy_fd, fdobj2 = resid_dy_fd)
  K_yy_i <- function(sevalarg, tevalarg) {
    eval.bifd(sevalarg = sevalarg, tevalarg = tevalarg, 
              bifd = k_yy_bifd)
  }
  
  K_mat_i <- function(sevalarg, tevalarg) {
    matrix(data = c(
      K_xx_i(sevalarg = sevalarg, tevalarg = tevalarg), 0,
      0, K_yy_i(sevalarg = sevalarg, tevalarg = tevalarg)
    ), nrow = 2, ncol = 2, byrow = TRUE)
  }
  if(!silent) {
    print("Doing Iteration 1 of Bias Correction")
  }
  VDP_Expectation_list_i <- mclapply(X = grid_points,
                                     FUN = VDP_Expectation_function,
                                     beta_xx = beta_xx_i,
                                     beta_xy = beta_xy_i,
                                     beta_yx = beta_yx_i,
                                     beta_yy = beta_yy_i, 
                                     K_st_mat_fun = K_mat_i,
                                     mc.cores = n_cores)
  for(tind in seq_along(grid_points)) {
    expectation_est_array[tind,,,1] <- VDP_Expectation_list_i[[tind]]
    bias_est_array[tind,,,1] <- ZtZ_inv_array[tind,,] %*% (N * (rbind(0, expectation_est_array[tind,,,1])))
  }
  for(j in 2:num_iter) {
    if(!silent) {
      print(paste0("Iteration ", j))
      print("Updating betas")
    }
    beta_array[,,,j] <- beta_array[,,,1] - bias_est_array[,,, (j - 1)]
    beta_xx_i <- approxfun(x = grid_points, y = beta_array[,2,1,j])
    beta_xy_i <- approxfun(x = grid_points, y = beta_array[,3,1,j])
    beta_yx_i <- approxfun(x = grid_points, y = beta_array[,2,2,j])
    beta_yy_i <- approxfun(x = grid_points, y = beta_array[,3,2,j])
    
    if(!silent) {print("Updating residuals")}
    # Update residuals:
    for(tind in seq_along(grid_points)) {
      design_matrix_tind <- model.matrix(~ x[tind,] + y[tind,])
      resid_array[tind,,,j] <- cbind(dx[tind,], dy[tind,]) -  design_matrix_tind %*% beta_array[tind,,,j]
    }
    
    # And add into functions
    resid_dx_fd <- Data2fd(argvals = grid_points,
                           y = resid_array[,,1,j])
    resid_dy_fd <- Data2fd(argvals = grid_points,
                           y = resid_array[,,2,j])
    K_xx_bifd <- var.fd(fdobj1 = resid_dx_fd, fdobj2 = resid_dx_fd)
    K_xy_bifd <- var.fd(fdobj1 = resid_dx_fd, fdobj2 = resid_dy_fd)
    K_yx_bifd <- var.fd(fdobj1 = resid_dy_fd, fdobj2 = resid_dx_fd)
    K_yy_bifd <- var.fd(fdobj1 = resid_dy_fd, fdobj2 = resid_dy_fd)
    
    
    if(!silent) {print("Estimating bias (in paralell)")}
    system.time(VDP_Expectation_list_ij <- mclapply(X = grid_points,
                                                    FUN = VDP_Expectation_function,
                                                    beta_xx = beta_xx_i,
                                                    beta_xy = beta_xy_i,
                                                    beta_yx = beta_yx_i,
                                                    beta_yy = beta_yy_i, 
                                                    K_st_mat_fun = K_mat_i,
                                                    mc.cores = n_cores))
    
    for(tind in seq_along(grid_points)) {
      expectation_est_array[tind,,,j] <- VDP_Expectation_list_ij[[tind]]
      bias_est_array[tind,,,j] <- ZtZ_inv_array[tind,,] %*% (N * (rbind(0, expectation_est_array[tind,,,j])))
    }
    
  }
  
  # Final beta subtracts final bias estimate:  
  beta_array[,,,num_iter + 1] <- beta_array[,,,1] - bias_est_array[,,, num_iter]
    
  
  
  list(beta_array = beta_array, bias_est_array = bias_est_array, resid_array = resid_array)
  }



