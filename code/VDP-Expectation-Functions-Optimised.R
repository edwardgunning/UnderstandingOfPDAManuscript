require(fda)

VDP_Expectation_function_opt <- function(t, 
                                        ST_11,
                                        ST_12,
                                        ST_21, 
                                        ST_22,
                                        K_11,
                                        K_12,
                                        K_21,
                                        K_22
                                        ) {
  
  
  
  # -------------------------------------------------------------------------
  integrand_11 <- function(sevalarg, tevalarg) {
    ST_11(sevalarg, tevalarg) * K_11(sevalarg, tevalarg) +
      ST_12(sevalarg, tevalarg) * K_21(sevalarg, tevalarg)
  }
  
  integrand_12 <- function(sevalarg, tevalarg) {
    ST_11(sevalarg, tevalarg) * K_12(sevalarg, tevalarg) +
      ST_12(sevalarg, tevalarg) * K_22(sevalarg, tevalarg)
  }
  
  integrand_21 <- function(sevalarg, tevalarg) {
    ST_21(sevalarg, tevalarg) * K_11(sevalarg, tevalarg) +
      ST_22(sevalarg, tevalarg) * K_21(sevalarg, tevalarg)
  }
  
  integrand_22 <- function(sevalarg, tevalarg) {
    ST_21(sevalarg, tevalarg) * K_12(sevalarg, tevalarg) +
      ST_22(sevalarg, tevalarg) * K_22(sevalarg, tevalarg)
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



# -------------------------------------------------------------------------

VDP_Expectation_function_biv_opt <- function(s, 
                                            t, 
                                            ST_11,
                                            ST_12,
                                            ST_21, 
                                            ST_22,
                                            K_11,
                                            K_12,
                                            K_21,
                                            K_22) {
  # -------------------------------------------------------------------------
  integrand_11 <- function(uevalarg, sevalarg, tevalarg) {
    ST_11(uevalarg, sevalarg) * K_11(uevalarg, tevalarg) +
      ST_12(uevalarg, sevalarg) * K_21(uevalarg, tevalarg)
  }
  
  integrand_12 <- function(uevalarg, sevalarg, tevalarg) {
    ST_11(uevalarg, sevalarg) * K_12(uevalarg, tevalarg) +
      ST_12(uevalarg, sevalarg) * K_22(uevalarg, tevalarg)
  }
  
  integrand_21 <- function(uevalarg, sevalarg, tevalarg) {
    ST_21(uevalarg, sevalarg) * K_11(uevalarg, tevalarg) +
      ST_22(uevalarg, sevalarg) * K_21(uevalarg, tevalarg)
  }
  
  integrand_22 <- function(uevalarg, sevalarg, tevalarg) {
    ST_21(uevalarg, sevalarg) * K_12(uevalarg, tevalarg) +
      ST_22(uevalarg, sevalarg) * K_22(uevalarg, tevalarg)
  }
  
  
  # Vectorise --------------------------------------------------------------
  integrand_11_vectorised <- Vectorize(integrand_11, vectorize.args = "uevalarg")
  integrand_12_vectorised <- Vectorize(integrand_12, vectorize.args = "uevalarg")
  integrand_21_vectorised <- Vectorize(integrand_21, vectorize.args = "uevalarg")
  integrand_22_vectorised <- Vectorize(integrand_22, vectorize.args = "uevalarg")
  
  
  # -------------------------------------------------------------------------
  matrix(data = c(
    integrate(f = integrand_11_vectorised, lower = 0, upper = s, sevalarg = s, tevalarg = t)$value,
    integrate(f = integrand_12_vectorised, lower = 0, upper = s, sevalarg = s, tevalarg = t)$value,
    integrate(f = integrand_21_vectorised, lower = 0, upper = s, sevalarg = s, tevalarg = t)$value,
    integrate(f = integrand_22_vectorised, lower = 0, upper = s, sevalarg = s, tevalarg = t)$value
  ), nrow = 2, ncol = 2, byrow = TRUE)
}



# -------------------------------------------------------------------------

do_pda_vdp_multiple_iterations_opt <- function(x, y, dx, dy, grid_points, num_iter, silent = FALSE, n_cores = 8) {
  
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
  K_yx_i <- K_xy_i <- function(sevalarg, tevalarg) {0}
  
  
  transition_matrix_sarg_to_targs_i <- function(sarg, targs) {
    
    dynamics_equations <- function(t, y, ...) {
      with(as.list(c(y)),{
        # rate of change
        dX <- beta_xx_i(t) * X + beta_xy_i(t) * Y
        dY <- beta_yx_i(t) * X + beta_yy_i(t) * Y
        list(c(dX, dY))
      }
      )
    }
    
    col1 <- deSolve::lsoda(y = c(X = 1, Y = 0),
                           times = c(sarg, targs),
                           func = dynamics_equations,
                           tcrit = max(targs)) #only take 2nd row not interested in initial value
    
    col2 <- deSolve::lsoda(y = c(X = 0, Y = 1),
                           times = c(sarg, targs),
                           func = dynamics_equations,
                           tcrit = max(targs)) #only take 2nd row not interested in initial value
    
    stopifnot(nrow(col1) == nrow(col2))
    
    
    col1 <- col1[-1,, drop = FALSE]
    col2 <- col2[-1,, drop = FALSE]
    
    return_array <- array(NA, dim = c(nrow(col1), 2, 2))
    
    for(tind in seq_len(nrow(col1))) {
      return_array[tind,,] <- cbind(col1[tind,-1], col2[tind,-1])
    }
    
    return_array
    
  }
  
  transition_array_sarg_to_targ_i <- array(NA, dim = c(
    D, D, 2, 2
  ))
  
  for(i in seq_len(D)) {
    jinds <- i:(D)
    transition_array_sarg_to_targ_i[i,jinds,,] <- transition_matrix_sarg_to_targs_i(sarg = grid_points[i], targs = grid_points[jinds])
  }
  
  for(j in seq_len(length.out = D - 1)) {
    for(i in (j+1):D) {
      transition_array_sarg_to_targ_i[i, j,,] <- solve(transition_array_sarg_to_targ_i[j,i,,])
    }
  }
  
  STM11_i <- transition_array_sarg_to_targ_i[,,1,1]
  STM12_i <- transition_array_sarg_to_targ_i[,,1,2]
  STM21_i <- transition_array_sarg_to_targ_i[,,2,1]
  STM22_i <- transition_array_sarg_to_targ_i[,,2,2]
  
  STM11_st_i <- function(sevalarg, tevalarg) {
    pracma::interp2(grid_points, grid_points, t(STM11_i), sevalarg, tevalarg)
  }
  
  STM12_st_i <- function(sevalarg, tevalarg ) {
    pracma::interp2(grid_points, grid_points, t(STM12_i), sevalarg, tevalarg)
  }
  
  STM21_st_i <- function(sevalarg, tevalarg ) {
    pracma::interp2(grid_points, grid_points, t(STM21_i), sevalarg, tevalarg)
  }
  
  STM22_st_i <- function(sevalarg, tevalarg ) {
    pracma::interp2(grid_points, grid_points, t(STM22_i), sevalarg, tevalarg)
  }
  
  if(!silent) {
    print("Doing Iteration 1 of Bias Correction")
  }
  VDP_Expectation_list_i <- mclapply(X = grid_points,
                                     FUN = VDP_Expectation_function_opt,
                                     ST_11 = STM11_st_i,
                                     ST_12 = STM12_st_i,
                                     ST_21 = STM21_st_i,
                                     ST_22 = STM22_st_i,
                                     K_11 = K_xx_i,
                                     K_12 = K_xy_i,
                                     K_21 = K_yx_i,
                                     K_22 = K_yy_i,
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
    
    transition_matrix_sarg_to_targs_i <- function(sarg, targs) {
      
      dynamics_equations <- function(t, y, ...) {
        with(as.list(c(y)),{
          # rate of change
          dX <- beta_xx_i(t) * X + beta_xy_i(t) * Y
          dY <- beta_yx_i(t) * X + beta_yy_i(t) * Y
          list(c(dX, dY))
        }
        )
      }
      
      col1 <- deSolve::lsoda(y = c(X = 1, Y = 0),
                             times = c(sarg, targs),
                             func = dynamics_equations,
                             tcrit = max(targs)) #only take 2nd row not interested in initial value
      
      col2 <- deSolve::lsoda(y = c(X = 0, Y = 1),
                             times = c(sarg, targs),
                             func = dynamics_equations,
                             tcrit = max(targs)) #only take 2nd row not interested in initial value
      
      stopifnot(nrow(col1) == nrow(col2))
      
      
      col1 <- col1[-1,, drop = FALSE]
      col2 <- col2[-1,, drop = FALSE]
      
      return_array <- array(NA, dim = c(nrow(col1), 2, 2))
      
      for(tind in seq_len(nrow(col1))) {
        return_array[tind,,] <- cbind(col1[tind,-1], col2[tind,-1])
      }
      
      return_array
      
    }
    
    
    transition_array_sarg_to_targ_i <- array(NA, dim = c(
      length(grid_points), length(grid_points), 2, 2
    ))
  
    
    for(i in seq_len(length(grid_points))) {
      jinds <- i:(length(grid_points))
      transition_array_sarg_to_targ_i[i,jinds,,] <- transition_matrix_sarg_to_targs_i(sarg = grid_points[i], targs = grid_points[jinds])
    }
    
    
    
    for(i2 in seq_len(length.out = length(grid_points) - 1)) {
      for(i in (i2+1):length(grid_points)) {
        transition_array_sarg_to_targ_i[i, i2,,] <- solve(transition_array_sarg_to_targ_i[i2,i,,])
      }
    }
    
    
    STM11_i <- transition_array_sarg_to_targ_i[,,1,1]
    STM12_i <- transition_array_sarg_to_targ_i[,,1,2]
    STM21_i <- transition_array_sarg_to_targ_i[,,2,1]
    STM22_i <- transition_array_sarg_to_targ_i[,,2,2]
    
    STM11_st_i <- function(sevalarg, tevalarg) {
      pracma::interp2(grid_points, grid_points, t(STM11_i), sevalarg, tevalarg)
    }
    
    STM12_st_i <- function(sevalarg, tevalarg ) {
      pracma::interp2(grid_points, grid_points, t(STM12_i), sevalarg, tevalarg)
    }
    
    STM21_st_i <- function(sevalarg, tevalarg ) {
      pracma::interp2(grid_points, grid_points, t(STM21_i), sevalarg, tevalarg)
    }
    
    STM22_st_i <- function(sevalarg, tevalarg ) {
      pracma::interp2(grid_points, grid_points, t(STM22_i), sevalarg, tevalarg)
    }
    
    
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
    VDP_Expectation_list_ij <- mclapply(X = grid_points,
                                        FUN = VDP_Expectation_function_opt,
                                        ST_11 = STM11_st_i,
                                        ST_12 = STM12_st_i,
                                        ST_21 = STM21_st_i,
                                        ST_22 = STM22_st_i,
                                        K_11 = K_xx_i,
                                        K_12 = K_xy_i,
                                        K_21 = K_yx_i,
                                        K_22 = K_yy_i,
                                        mc.cores = n_cores)
    
    for(tind in seq_along(grid_points)) {
      expectation_est_array[tind,,,j] <- VDP_Expectation_list_ij[[tind]]
      bias_est_array[tind,,,j] <- ZtZ_inv_array[tind,,] %*% (N * (rbind(0, expectation_est_array[tind,,,j])))
    }
    
  }
  
  # Final beta subtracts final bias estimate:  
  beta_array[,,,num_iter + 1] <- beta_array[,,,1] - bias_est_array[,,, num_iter]
  
  list(beta_array = beta_array, bias_est_array = bias_est_array, resid_array = resid_array)
}






