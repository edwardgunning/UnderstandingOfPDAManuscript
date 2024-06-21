generate_random_noise_curves <- function(grid_points, n, sigma) {
  n <- as.integer(n)
  stopifnot(sigma > 0 & n > 0 & is.numeric(grid_points))
  
  # Matrix of outer product |s-t|
  s_minus_t_mat <- outer(X = grid_points, Y = grid_points, FUN = "-")
  abs_s_minus_t_mat <- abs(s_minus_t_mat)
  
  # $K(|t-s|) = \sigma^2\phi(2*|t-s|)$:
  K <- sigma ^ 2 * dnorm(x = 2 * abs_s_minus_t_mat)
  
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

mu_init <- c(1, 0)
sigma_init <- diag(c(0.1, 0.1))




