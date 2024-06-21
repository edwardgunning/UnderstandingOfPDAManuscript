source(here::here("code", "SHM-bias-calculation-function.R"))

# Test with SHM Model: -----------------------------------------------------
test_grid_points <- seq(0, 2 * pi, length.out = 100)
test_K_st <-  function(s, t) {
  0.1 ^ 2 * dnorm(x = 2 * abs(t - s))
}
test_beta_u <- function(u) {
    - 1 + (u - u)
}

system.time(test_result <- bias_calculation(beta_u_est = test_beta_u,
                                K_st_est = test_K_st,
                                grid_points = test_grid_points))

system.time(test_result_tm_fast <- bias_calculation_with_transition_matrix_faster(beta_u_est = test_beta_u,
                                            K_st_est = test_K_st,
                                            grid_points = test_grid_points))

system.time(test_result_tm <- bias_calculation_with_transition_matrix(beta_u_est = test_beta_u,
                                K_st_est = test_K_st,
                                grid_points = test_grid_points))

# VisualCheck of results
par(mfrow = c(2, 2))
plot(test_result[1, 1, ], type = "l")
lines(test_result_tm[1, 1, ], col = 2)
lines(test_result_tm_fast[1, 1, ], col = 3)

plot(test_result[1, 2, ], type = "l")
lines(test_result_tm[1, 2, ], col = 2)
lines(test_result_tm_fast[1, 2, ], col = 3)

plot(test_result[2, 1, ], type = "l")
lines(test_result_tm[2, 1, ], col = 2)
lines(test_result_tm_fast[2, 1, ], col = 3)

plot(test_result[2, 2, ], type = "l")
lines(test_result_tm[2, 2, ], col = 2)
lines(test_result_tm_fast[2, 2, ], col = 3)

# Test with known functions for K(s, t) -----------------------------------
# Have known solutions here.

# K(s, t) = t*s
test_2_K_st <- function(s, t) {
 t * s
}

integral_12 <- function(t) {
  - t * (sin(t) - t * cos(t))
}
integral_22 <- function(t) {
  t * ( t * sin(t) + cos(t) - 1)
}
bias_12 <- function(t) {
  cos(t) * integral_12(t) + sin(t) * integral_22(t)
}
bias_22 <- function(t) {
  - sin (t) * integral_12(t) + cos(t) * integral_22(t)
}

system.time(test_2_result <- bias_calculation(beta_u_est = test_beta_u,
                                K_st_est = test_2_K_st,
                                grid_points = test_grid_points))

system.time(test_2_result_tm_fast <- bias_calculation_with_transition_matrix_faster(beta_u_est = test_beta_u,
                                              K_st_est = test_2_K_st,
                                              grid_points = test_grid_points)) # GIVES ERROR

system.time(test_2_result_tm <- bias_calculation_with_transition_matrix(beta_u_est = test_beta_u,
                                  K_st_est = test_2_K_st,
                                  grid_points = test_grid_points))

par(mfrow = c(2, 2))
plot(test_2_result[1, 1, ], type = "l")
lines(test_2_result_tm[1, 1, ], col = 2)
plot(test_2_result[1, 2, ], type = "l")
lines(bias_12(test_grid_points), col = 2, type = "p", cex = 0.5)
lines(test_2_result_tm[1, 2,], col = 4, type = "l")

plot(test_2_result[2, 1, ], type = "l")
plot(test_2_result[2, 2, ], type = "l")
lines(bias_22(test_grid_points), col = 2, type = "p", cex = 0.5)
lines(test_2_result_tm[2, 2,], col = 4, type = "l")



# Test 3 ------------------------------------------------------------------

test_3_K_st <- function(s, t) {
  t - s
}

integral_12 <- function(t) {
  - t + sin(t)
}
integral_22 <- function(t) {
  - cos(t) + 1
}
bias_12 <- function(t) {
  cos(t) * integral_12(t) + sin(t) * integral_22(t)
}
bias_22 <- function(t) {
  - sin (t) * integral_12(t) + cos(t) * integral_22(t)
}

system.time(test_3_result <- bias_calculation(beta_u_est = test_beta_u,
                                  K_st_est = test_3_K_st,
                                  grid_points = test_grid_points))
system.time(test_3_result_tm_fast <- bias_calculation_with_transition_matrix_faster(beta_u_est = test_beta_u,
                                                                        K_st_est = test_3_K_st,
                                                                        grid_points = test_grid_points))
test_3_result_tm <- bias_calculation_with_transition_matrix(beta_u_est = test_beta_u,
                                  K_st_est = test_3_K_st,
                                  grid_points = test_grid_points)
par(mfrow = c(2, 2))
plot(test_3_result[1, 1, ], type = "l")
plot(test_3_result[1, 2, ], type = "l")
lines(bias_12(test_grid_points), col = 2, type = "p", cex = 0.5)
lines(test_3_result_tm[1, 2,], col = 4)
plot(test_3_result[2, 1, ], type = "l")
lines(test_3_result_tm[2, 1,], col = 4)
plot(test_3_result[2, 2, ], type = "l")
lines(bias_22(test_grid_points), col = 2, type = "p", cex = 0.5)
lines(test_3_result_tm[2, 2,], col = 4)



# and just for timing -----------------------------------------------------

test_beta_u <- function(u) {
  u + 0.001 * u^2
}


system.time(test_result <- bias_calculation(beta_u_est = test_beta_u,
                                            K_st_est = test_K_st,
                                            grid_points = test_grid_points))

system.time(test_result_tm <- bias_calculation_with_transition_matrix(beta_u_est = test_beta_u,
                                                                      K_st_est = test_K_st,
                                                                      grid_points = test_grid_points))
  