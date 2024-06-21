library(ggplot2)    # CRAN v3.3.5
library(plot3D)     # CRAN v1.4
library(tikzDevice) # CRAN v0.12.3.1

cov_xs_xt <- readRDS(file = here::here("outputs",
                                       "SHM",
                                       "simulation-results",
                                       "cov_xs_xt_new.rds"))
n_grid <- 100
grid_range <- c(0, 2 * pi)
grid_points <- seq(from = grid_range[1],
                   to = grid_range[2],
                   length.out = n_grid)
M <- mesh(grid_points, grid_points)
test_grid <- expand.grid(grid_points, grid_points)


# -------------------------------------------------------------------------
zero_state_covariance_df <- cbind(test_grid, sapply(cov_xs_xt, function(x) {x[1,1]}))
names(zero_state_covariance_df) <- c("s", "t", "cov")
zero_state_covariance_df_wide <- reshape2::dcast(data = zero_state_covariance_df,
                                                 formula = s ~ t,
                                                 value.var = "cov")
stopifnot(zero_state_covariance_df_wide$s == grid_points)
stopifnot(colnames(zero_state_covariance_df_wide)[-1] ==grid_points)
zero_state_covariance_matrix_wide <- as.matrix(zero_state_covariance_df_wide[, -1])

# check symmetry:
stopifnot(max(zero_state_covariance_matrix_wide - t(zero_state_covariance_matrix_wide)) < 10^-8)
zero_state_covariance_matrix_wide  <- (zero_state_covariance_matrix_wide + t(zero_state_covariance_matrix_wide))/2

# -------------------------------------------------------------------------
phi_t0 <- function(t) {
  matrix(data = c(
    cos(t),
    sin(t),
    -sin(t),
    cos(t)
  ), nrow = 2, ncol = 2, byrow = TRUE)
} 

zero_input_covariance <- lapply(seq_len(nrow(test_grid)), function(i) {
  phi_t0(t = test_grid[i, 1]) %*% diag(c(0.05, 0.05)) %*% t(phi_t0(t = test_grid[i, 2]))
})

zero_input_covariance_11 <- sapply(zero_input_covariance, function(x) {
  x[1, 1]
})

zero_input_covariance_11_df <- cbind(test_grid, zero_input_covariance_11)
names(zero_input_covariance_11_df) <- c("s", "t", "cov")

zero_input_covariance_11_df_wide <- reshape2::dcast(data = zero_input_covariance_11_df,
                           formula = s ~ t,
                           value.var = "cov")
stopifnot(zero_input_covariance_11_df_wide[, 1]  == grid_points)
stopifnot(colnames(zero_input_covariance_11_df_wide)[-1] ==grid_points)
zero_input_covariance_11_matrix_wide <- as.matrix(zero_input_covariance_11_df_wide [, -1])

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937

# -------------------------------------------------------------------------
tikz(here::here("outputs", "SHM", "paper-plots", "cov_surface.tex"), 
     width = (1 * doc_width_inches),
     height = (1.3 * doc_width_inches)/2)
par(mfrow = c(1, 2), cex.main = 0.9, font.main = 2)
surf3D(x = M$x,
       y = M$y, 
       z = zero_input_covariance_11_matrix_wide,
       box = TRUE, 
       phi = 20,
       cex = 0.5,
       xlab = "$s$",
       ylab = "$t$",
       zlab = "Cov",
       main = "(a) Zero-input covariance",
       bty = "f",
       resfac = 3,
       cex = 0.1,
       colkey = list(side = 1, cex.clab = 0.5),
       theta = 100)
surf3D(x = M$x,
       y = M$y, 
       z = zero_state_covariance_matrix_wide,
       box = TRUE, 
       phi = 20,
       resfac = 3,
       cex = 0.9,
       bty = "f",
       ylab = "$t$",
       zlab = "Cov",
       xlab = "$s$",
       main = "(b) Zero-state covariance",
       colkey = list(side = 1),
       theta = 100)
dev.off()

