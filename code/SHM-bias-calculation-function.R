bias_calculation <- function(beta_u_est, K_st_est, grid_points) {
  
  K_st <- K_st_est
  beta_u <- beta_u_est
  # Some routine checks:
  stopifnot(is.numeric(grid_points) & length(grid_points >= 1))
  # check if beta_u can be evaluated at all gridpoints
  stopifnot(sapply(grid_points, FUN = function(y) {
    beta_u_eval <- beta_u(u = y)
    is.numeric(beta_u_eval) & (length(beta_u_eval) == 1)
  }))
  # and do same for K(s, t)
  grid_test <- expand.grid(grid_points, grid_points)
  stopifnot(all(do.call(mapply, c(function(a, b) {
    K_st_eval <- K_st(s = a, t = b)
    is.numeric(K_st_eval) & (length(K_st_eval) == 1)
  },
  unname(grid_test)))))
  
  # 1) $\bar{B} (s)$  the antiderivative of coefficient matrix
  B_bar_s <- function(s) {
    int_beta_t <- ifelse(test = s==0,
                         yes =  0,
                         no = integrate(f = beta_u, lower = 0, upper = s)$value)
    
    matrix(data = c(0,  int_beta_t, s, 0), 
           nrow = 2, 
           ncol = 2, 
           byrow = FALSE)
  }
  
  expand.grid(grid_points, grid_points)
  
  # 2) $\Sigma (s, t)$ the multivariate error covariance ------------
  Sigma_st <- function(s, t) {
    matrix(data = c(0, 0,
                    0, K_st(s = s, t = t)),
           nrow = 2, ncol = 2, 
           byrow = TRUE)
  }
  
  
  # 3)The integrand $exp(- \bar{B} (s)) \Sigma (s, t) ------------
  integrand_st <- function(s, t) {
    as.matrix(Matrix::expm(x = - B_bar_s(s = s)) %*% Sigma_st(s = s, t = t))
  }
  
  
  # 4) Split integrand matrix into entries so we can integrate separ --------
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
  
  
  # 5) Integrate the entries
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
  
  
  
  # Function to compute integral matrix -------------------------------------
  integral_matrix <- function(t) {
    matrix(data = c(
      integral_t_11(t = t), integral_t_12(t = t),
      integral_t_21(t = t), integral_t_22(t = t)
    ), nrow = 2, ncol = 2, byrow = TRUE)
  }
  
  
  
  # And compute full bias ---------------------------------------------------
  bias_matrix <- function(t) {
    as.matrix(Matrix::expm(x = B_bar_s(s = t))) %*% integral_matrix(t = t)
  }
  
  bias_matrix_return_array <- function(t_grid) {
    return_array <- array(NA, dim = c(2, 2, length(t_grid)))
    for(tind in seq_along(t_grid)) {
      return_array[,,tind] <- bias_matrix(t = t_grid[tind])
    }
    return_array
  }
  
  
  # final return: -----------------------------------------------------------
  bias_matrix_return_array(t_grid = grid_points)
}





# Now for general transition matrix: --------------------------------------

bias_calculation_with_transition_matrix <- function(beta_u_est, K_st_est, grid_points) {
  
  K_st <- K_st_est
  beta_u <- beta_u_est
  # Some routine checks:
  stopifnot(is.numeric(grid_points) & length(grid_points >= 1))
  # check if beta_u can be evaluated at all gridpoints
  stopifnot(sapply(grid_points, FUN = function(y) {
    beta_u_eval <- beta_u(u = y)
    is.numeric(beta_u_eval) & (length(beta_u_eval) == 1)
  }))
  # and do same for K(s, t)
  grid_test <- expand.grid(grid_points, grid_points)
  stopifnot(all(do.call(mapply, c(function(a, b) {
    K_st_eval <- K_st(s = a, t = b)
    is.numeric(K_st_eval) & (length(K_st_eval) == 1)
  },
  unname(grid_test)))))
  
  # 1) TRANSIOTION MATRIX
  transition_matrix_sarg_to_targ <- function(sarg, targ) {
    
    dynamics_equations <- function(t, y, ...) {
      with(as.list(c(y)),{
        # rate of change
        dX <- Y
        dY <-  beta_u(t) * X  
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
  
  
  
  # 2) $\Sigma (s, t)$ the multivariate error covariance ------------
  Sigma_st <- function(s, t) {
    matrix(data = c(0, 0,
                    0, K_st(s = s, t = t)),
           nrow = 2, ncol = 2, 
           byrow = TRUE)
  }
  
  
  # 3)The integrand $exp(- \bar{B} (s)) \Sigma (s, t) ------------
  integrand_st <- function(s, t) {
    transition_matrix_sarg_to_targ(sarg = s, targ = t) %*% Sigma_st(s = s, t = t)
  }
  
  
  # 4) Split integrand matrix into entries so we can integrate separ --------
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
  
  
  # 5) Integrate the entries
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
  
  
  
  # Function to compute integral matrix -------------------------------------
  integral_matrix <- function(t) {
    matrix(data = c(
      integral_t_11(t = t), integral_t_12(t = t),
      integral_t_21(t = t), integral_t_22(t = t)
    ), nrow = 2, ncol = 2, byrow = TRUE)
  }
  
  
  
  # And compute full bias ---------------------------------------------------
  
  bias_matrix_return_array <- function(t_grid) {
    return_array <- array(NA, dim = c(2, 2, length(t_grid)))
    for(tind in seq_along(t_grid)) {
      return_array[,,tind] <- integral_matrix(t = t_grid[tind])
    }
    return_array
  }
  
  
  # final return: -----------------------------------------------------------
  bias_matrix_return_array(t_grid = grid_points)
}



# bias_calculation_with_transition_matrix_2 <- function(beta_u_est, K_st_est, grid_points) {
#   
#   K_st <- K_st_est
#   beta_u <- beta_u_est
#   # Some routine checks:
#   stopifnot(is.numeric(grid_points) & length(grid_points >= 1))
#   # check if beta_u can be evaluated at all gridpoints
#   stopifnot(sapply(grid_points, FUN = function(y) {
#     beta_u_eval <- beta_u(u = y)
#     is.numeric(beta_u_eval) & (length(beta_u_eval) == 1)
#   }))
#   # and do same for K(s, t)
#   grid_test <- expand.grid(grid_points, grid_points)
#   stopifnot(all(do.call(mapply, c(function(a, b) {
#     K_st_eval <- K_st(s = a, t = b)
#     is.numeric(K_st_eval) & (length(K_st_eval) == 1)
#   },
#   unname(grid_test)))))
#   
#   # 1) TRANSIOTION MATRIX
#   transition_matrix_sarg_to_targ <- function(sarg, targ) {
#     
#     dynamics_equations <- function(t, y, ...) {
#       with(as.list(c(y)),{
#         # rate of change
#         dX <- Y
#         dY <-  beta_u(t) * X  
#         list(c(dX, dY))
#       }
#       )
#     }
#     
#     col1 <- deSolve::lsoda(y = c(X = 1, Y = 0),
#                            times = c(sarg, targ),
#                            func = dynamics_equations,
#                            tcrit = targ)#only take 2nd row not interested in initial value
#     
#     col2 <- deSolve::lsoda(y = c(X = 0, Y = 1),
#                            times = c(sarg, targ),
#                            func = dynamics_equations,
#                            tcrit = targ)#only take 2nd row not interested in initial value
#     
#     cbind(col1[2,-1], col2[2,-1])
#     
#   }
#   
# 
#   # Set up grid -------------------------------------------------------------
#   
# 
#   
#   
#   
#   
#   # 2) $\Sigma (s, t)$ the multivariate error covariance ------------
#   Sigma_st <- function(s, t) {
#     matrix(data = c(0, 0,
#                     0, K_st(s = s, t = t)),
#            nrow = 2, ncol = 2, 
#            byrow = TRUE)
#   }
#   
#   
#   # 3)The integrand $exp(- \bar{B} (s)) \Sigma (s, t) ------------
#   integrand_st <- function(s, t) {
#     transition_matrix_sarg_to_targ(sarg = s, targ = t) %*% Sigma_st(s = s, t = t)
#   }
#   
#   
#   # 4) Split integrand matrix into entries so we can integrate separ --------
#   # THESE MUST BE VECTORISED IN S SO THAT THEY CAN BE INTEGRATED ALONG S
#   integrand_st_11 <- function(s, t) {
#     return_vector <- vector(mode = "numeric", length = length(s))
#     for(i in seq_along(s)) {
#       return_vector[i] <- as.matrix(integrand_st(s = s[i], t = t))[1, 1]
#     }
#     return_vector
#   }
#   
#   integrand_st_12 <- function(s, t) {
#     return_vector <- vector(mode = "numeric", length = length(s))
#     for(i in seq_along(s)) {
#       return_vector[i] <- as.matrix(integrand_st(s = s[i], t = t))[1, 2]
#     }
#     return_vector
#   }
#   
#   integrand_st_21 <- function(s, t) {
#     return_vector <- vector(mode = "numeric", length = length(s))
#     for(i in seq_along(s)) {
#       return_vector[i] <- as.matrix(integrand_st(s = s[i], t = t))[2, 1]
#     }
#     return_vector
#   }
#   
#   integrand_st_22 <- function(s, t) {
#     return_vector <- vector(mode = "numeric", length = length(s))
#     for(i in seq_along(s)) {
#       return_vector[i] <- as.matrix(integrand_st(s = s[i], t = t))[2, 2]
#     }
#     return_vector
#   }
#   
#   
#   # 5) Integrate the entries
#   integral_t_11 <- function(t) {
#     integrate(f = integrand_st_11, lower = 0, upper = t, t = t)$value
#   }
#   
#   integral_t_12 <- function(t) {
#     integrate(f = integrand_st_12, lower = 0, upper = t, t = t)$value
#   }
#   
#   integral_t_21 <- function(t) {
#     integrate(f = integrand_st_21, lower = 0, upper = t, t = t)$value
#   }
#   
#   integral_t_22 <- function(t) {
#     integrate(f = integrand_st_22, lower = 0, upper = t, t = t)$value
#   }
#   
#   
#   
#   # Function to compute integral matrix -------------------------------------
#   integral_matrix <- function(t) {
#     matrix(data = c(
#       integral_t_11(t = t), integral_t_12(t = t),
#       integral_t_21(t = t), integral_t_22(t = t)
#     ), nrow = 2, ncol = 2, byrow = TRUE)
#   }
#   
#   
#   
#   # And compute full bias ---------------------------------------------------
#   
#   bias_matrix_return_array <- function(t_grid) {
#     return_array <- array(NA, dim = c(2, 2, length(t_grid)))
#     for(tind in seq_along(t_grid)) {
#       return_array[,,tind] <- integral_matrix(t = t_grid[tind])
#     }
#     return_array
#   }
#   
#   
#   # final return: -----------------------------------------------------------
#   bias_matrix_return_array(t_grid = grid_points)
# }
# 
