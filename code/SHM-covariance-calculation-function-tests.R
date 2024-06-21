source(here::here("code", "SHM-covariance-calculation-function.R"))
library(parallel) # Support for Parallel computation in R 

# Test function with known functions plugged in for K(s, t) ---------------
# these may not be valid covariances...
# but will try to catch any obvious errors in function.

# -------------------------------------------------------------------------

# K(s, t) = 1
test_grid_points <- seq(0, 2 * pi, length.out = 25)
test_grid <- expand.grid(test_grid_points, test_grid_points)

# https://www.wolframalpha.com/input?i2d=true&i=Integrate%5B%7B%7Bcos%5C%2840%29s%5C%2841%29%2Csin%5C%2840%29s%5C%2841%29%7D%2C%7B-sin%5C%2840%29s%5C%2841%29%2Ccos%5C%2840%29s%5C%2841%29%7D%7D%7B%7Bcos%5C%2840%29u%5C%2841%29%2C-sin%5C%2840%29u%5C%2841%29%7D%2C%7Bsin%5C%2840%29u%5C%2841%29%2Ccos%5C%2840%29u%5C%2841%29%7D%7D%7B%7B0%2C0%7D%2C%7B0%2C1%7D%7D%7B%7Bcos%5C%2840%29v%5C%2841%29%2Csin%5C%2840%29v%5C%2841%29%7D%2C%7B-sin%5C%2840%29v%5C%2841%29%2Ccos%5C%2840%29v%5C%2841%29%7D%7D%7B%7Bcos%5C%2840%29t%5C%2841%29%2C-sin%5C%2840%29t%5C%2841%29%7D%2C%7Bsin%5C%2840%29t%5C%2841%29%2Ccos%5C%2840%29t%5C%2841%29%7D%7D%2C%7Bu%2C0%2Cs%7D%2C%7Bv%2C0%2Ct%7D%5D+
test_K_st <- function(s, t) {
  (t * s) - (t * s) + 1
}

beta_u_true <- function(u) {
  - 1 + (u - u)
}


ncores <- detectCores()

system.time(test <- mclapply(seq_len(nrow(test_grid)), function (i) {
  calculate_covariance_shm_true(K_st_est = test_K_st,
                           s = test_grid[i, 1],
                           t = test_grid[i, 2])},
  mc.cores = ncores))

system.time(test_est <- mclapply(seq_len(nrow(test_grid)), function (i) {
  calculate_covariance_shm_est(K_st_est = test_K_st,
                                         beta_u_est = beta_u_true,
                                         s = test_grid[i, 1],
                                         test_grid[i, 2])},
                     mc.cores = ncores))


truth_11 <- function(s, t) {
  (cos(s) - 1) * (cos(t)-1)
}
truth_12 <- function(s, t) {
  - (cos(s) - 1) * sin(t)
}

truth_21 <- function(s, t) {
  sin(s) * (-(cos(t)-1))
} 

truth_22 <- function(s, t) {
  sin(s) * sin(t)
} 


# Compare tests and truths ------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
test_11 <- sapply(test, function(x) {
  x[1, 1]
})

test_est_11 <- sapply(test_est, function(x) {
  x[1, 1]
})

truth_11_eval <- sapply(seq_len(nrow(test_grid)), function(i) {
  truth_11(s = test_grid[i, 1],
           t = test_grid[i, 2])})

plot(test_11, truth_11_eval)
abline(a=0, b=1)
plot(test_est_11, truth_11_eval)
abline(a=0, b=1)

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
test_12 <- sapply(test, function(x) {
  x[1, 2]
})
test_est_12 <- sapply(test_est, function(x) {
  x[1, 2]
})
truth_12_eval <- sapply(seq_len(nrow(test_grid)), function(i) {
  truth_12(s = test_grid[i, 1],
           t = test_grid[i, 2])})
plot(test_12, truth_12_eval)
abline(a=0, b=1)
plot(test_est_12, truth_12_eval)
abline(a=0, b=1)

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
test_21 <- sapply(test, function(x) {
  x[2, 1]
})
test_est_21 <- sapply(test_est, function(x) {
  x[2, 1]
})
truth_21_eval <- sapply(seq_len(nrow(test_grid)), function(i) {
  truth_21(s = test_grid[i, 1],
           t = test_grid[i, 2])})
plot(test_21, truth_21_eval)
abline(a = 0, b = 1)
plot(test_est_21, truth_21_eval)
abline(a = 0, b = 1)

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
test_22 <- sapply(test, function(x) {
  x[2, 2]
})
test_est_22 <- sapply(test_est, function(x) {
  x[2, 2]
})
truth_22_eval <- sapply(seq_len(nrow(test_grid)), function(i) {
  truth_22(s = test_grid[i, 1],
           t = test_grid[i, 2])})
plot(test_22, truth_22_eval)
abline(a=0, b =1)
plot(test_est_22, truth_22_eval)
abline(a=0, b =1)

