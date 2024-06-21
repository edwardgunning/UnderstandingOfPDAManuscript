# True covariance function for SHM model.
# Re;lies on properties of time-invariant ODE solutions
calculate_covariance_shm_true <- function(K_st_est, s, t) {
  
  B_bar_q <- function(q) {
    matrix(data = c(0,  
                    -q,
                    q,
                    0), 
           nrow = 2, 
           ncol = 2, 
           byrow = FALSE)
  }
  
  integrand_uv_11 <- function(u, v) {
    (Matrix::expm(x = - B_bar_q(u)) %*%
       matrix(c(0, 0, 0, K_st_est(u, v)), ncol = 2, nrow = 2) %*%
       t(as.matrix(Matrix::expm(x = - B_bar_q(v)))))[1, 1]
  }
  integrand_uv_12 <- function(u, v) {
    (Matrix::expm(x = - B_bar_q(u)) %*%
       matrix(c(0, 0, 0, K_st_est(u, v)), ncol = 2, nrow = 2) %*%
       t(as.matrix(Matrix::expm(x = - B_bar_q(v)))))[1, 2]
  }
  integrand_uv_21 <- function(u, v) {
    (Matrix::expm(x = - B_bar_q(u)) %*%
       matrix(c(0, 0, 0, K_st_est(u, v)), ncol = 2, nrow = 2) %*%
       t(as.matrix(Matrix::expm(x = - B_bar_q(v)))))[2, 1]
  }
  integrand_uv_22 <- function(u, v) {
    (Matrix::expm(x = - B_bar_q(u)) %*%
       matrix(c(0, 0, 0, K_st_est(u, v)), ncol = 2, nrow = 2) %*%
       t(as.matrix(Matrix::expm(x = - B_bar_q(v)))))[2, 2]
  }
  
  # VECTORISE INTEGRANDS ----------------------------------------------------
  # Over u:
  integrand_uv_11_vectorise_u <- Vectorize(integrand_uv_11, vectorize.args = "u")
  integrand_uv_12_vectorise_u <- Vectorize(integrand_uv_12, vectorize.args = "u")
  integrand_uv_21_vectorise_u <- Vectorize(integrand_uv_21, vectorize.args = "u")
  integrand_uv_22_vectorise_u <- Vectorize(integrand_uv_22, vectorize.args = "u")
  
  
  # PART INTEGRATE integrands -----------------------------------------------
  # over u
  part_integral_vsarg_11 <- function(v, sarg) {
    integrate(f = integrand_uv_11_vectorise_u, lower = 0, upper = sarg, v = v)$value
  }
  part_integral_vsarg_12 <- function(v, sarg) {
    integrate(f = integrand_uv_12_vectorise_u, lower = 0, upper = sarg, v = v)$value
  }
  part_integral_vsarg_21 <- function(v, sarg) {
    integrate(f = integrand_uv_21_vectorise_u, lower = 0, upper = sarg, v = v)$value
  }
  part_integral_vsarg_22 <- function(v, sarg) {
    integrate(f = integrand_uv_22_vectorise_u, lower = 0, upper = sarg, v = v)$value
  }
  
  # VECTORISE PART INTEGRALS ------------------------------------------------
  # over v
  part_integral_vsarg_11_vectorise_v <- Vectorize(part_integral_vsarg_11, vectorize.args = "v")
  part_integral_vsarg_12_vectorise_v <- Vectorize(part_integral_vsarg_12, vectorize.args = "v")
  part_integral_vsarg_21_vectorise_v <- Vectorize(part_integral_vsarg_21, vectorize.args = "v")
  part_integral_vsarg_22_vectorise_v <- Vectorize(part_integral_vsarg_22, vectorize.args = "v")
  
  
  full_integral_sargtarg_matrix <- function(sarg, targ) {
    matrix(data = c(
      integrate(part_integral_vsarg_11_vectorise_v, lower = 0, upper = targ, sarg = sarg)$value,
      integrate(part_integral_vsarg_12_vectorise_v, lower = 0, upper = targ, sarg = sarg)$value,
      integrate(part_integral_vsarg_21_vectorise_v, lower = 0, upper = targ, sarg = sarg)$value,
      integrate(part_integral_vsarg_22_vectorise_v, lower = 0, upper = targ, sarg = sarg)$value
    ), nrow = 2, ncol = 2, byrow = TRUE)
  }
  
  as.matrix(Matrix::expm(x = B_bar_q(q = s))) %*%
    full_integral_sargtarg_matrix(sarg = s, targ = t) %*%
    t(as.matrix(Matrix::expm(x = B_bar_q(q = t))))
}


# -------------------------------------------------------------------------

calculate_covariance_shm_est <- function(K_st_est, beta_u_est, s, t) {
  
  # 1) TRANSIOTION MATRIX
  transition_matrix_sarg_to_targ <- function(sarg, targ) {
    
    dynamics_equations <- function(t, y, ...) {
      with(as.list(c(y)),{
        # rate of change
        dX <- Y
        dY <-  beta_u_est(t) * X  
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
  
  integrand_function_uvst <- function(u, v, t, s) {
    transition_matrix_sarg_to_targ(sarg = u, targ = s) %*%
      matrix(c(0, 0, 0, K_st_est(u, v)), ncol = 2, nrow = 2) %*%
      t(transition_matrix_sarg_to_targ(sarg = v, targ = t))
  }
  
  integrand_uvst_11 <- function(u, v, s, t) {
    integrand_function_uvst(u = u, v = v, s = s, t = t)[1, 1]
  }
  integrand_uvst_12 <- function(u, v, s, t) {
    integrand_function_uvst(u = u, v = v, s = s, t = t)[1, 2]
  }
  integrand_uvst_21 <- function(u, v, s, t) {
    integrand_function_uvst(u = u, v = v, s = s, t = t)[2, 1]
  }
  integrand_uvst_22 <- function(u, v, s, t) {
    integrand_function_uvst(u = u, v = v, s = s, t = t)[2, 2]
  }
  
  # VECTORISE INTEGRANDS ----------------------------------------------------
  # Over u:
  integrand_uvst_11_vectorise_u <- Vectorize(integrand_uvst_11, vectorize.args = "u")
  integrand_uvst_12_vectorise_u <- Vectorize(integrand_uvst_12, vectorize.args = "u")
  integrand_uvst_21_vectorise_u <- Vectorize(integrand_uvst_21, vectorize.args = "u")
  integrand_uvst_22_vectorise_u <- Vectorize(integrand_uvst_22, vectorize.args = "u")
  
  
  # PART INTEGRATE integrands -----------------------------------------------
  # over u
  part_integral_vsargtarg_11 <- function(v, sarg, targ) {
    integrate(f = integrand_uvst_11_vectorise_u, lower = 0, upper = sarg, v = v, s = sarg, t = targ)$value
  }
  part_integral_vsargtarg_12 <- function(v, sarg, targ) {
    integrate(f = integrand_uvst_12_vectorise_u, lower = 0, upper = sarg, v = v, s = sarg, t = targ)$value
  }
  part_integral_vsargtarg_21 <- function(v, sarg, targ) {
    integrate(f = integrand_uvst_21_vectorise_u, lower = 0, upper = sarg, v = v, s = sarg, t = targ)$value
  }
  part_integral_vsargtarg_22 <- function(v, sarg, targ) {
    integrate(f = integrand_uvst_22_vectorise_u, lower = 0, upper = sarg, v = v, s = sarg, t = targ)$value
  }
  
  # VECTORISE PART INTEGRALS ------------------------------------------------
  # over v
  part_integral_vsargtarg_11_vectorise_v <- Vectorize(part_integral_vsargtarg_11, vectorize.args = "v")
  part_integral_vsargtarg_12_vectorise_v <- Vectorize(part_integral_vsargtarg_12, vectorize.args = "v")
  part_integral_vsargtarg_21_vectorise_v <- Vectorize(part_integral_vsargtarg_21, vectorize.args = "v")
  part_integral_vsargtarg_22_vectorise_v <- Vectorize(part_integral_vsargtarg_22, vectorize.args = "v")
  
  
  full_integral_sargtarg_matrix <- function(sarg, targ) {
    matrix(data = c(
      integrate(part_integral_vsargtarg_11_vectorise_v, lower = 0, upper = targ, sarg = sarg, targ = targ)$value,
      integrate(part_integral_vsargtarg_12_vectorise_v, lower = 0, upper = targ, sarg = sarg, targ = targ)$value,
      integrate(part_integral_vsargtarg_21_vectorise_v, lower = 0, upper = targ, sarg = sarg, targ = targ)$value,
      integrate(part_integral_vsargtarg_22_vectorise_v, lower = 0, upper = targ, sarg = sarg, targ = targ)$value
    ), nrow = 2,
    ncol = 2, 
    byrow = TRUE)
  }

    full_integral_sargtarg_matrix(sarg = s, targ = t)
}

