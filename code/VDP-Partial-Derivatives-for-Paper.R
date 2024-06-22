# Load packages and required functions: -----------------------------------
source(here::here("code", "VDP-helper-functions.R"))
library(fda)        # CRAN v5.5.1
library(refund)     # CRAN v0.1-26
library(data.table) # CRAN v1.14.2

# Set up settings: --------------------------------------------------------
Nrep <- 50

# Vary different parameters one at a time: --------------------------------
settings_vary_sigma_forcing <- expand.grid(rep = seq_len(Nrep),
                                           mu = 1,
                                           intensity_forcing = 2,
                                           sigma_forcing = 0.1)
settings <- settings_vary_sigma_forcing
N_sim_total <- nrow(settings)
# Parameters that will be fixed: ------------------------------------------
mu_ic <- c(1.9922213675594, -0.910974076470711)
Sigma_ic <- matrix(c(0.025, 0, 0, 0.025),
                   nrow = 2, 
                   ncol = 2) # Not too much IC noise
grid_points <- seq(0, 13, length.out = 200) 
n <- 200 # number of curves
results_list <- vector(mode = "list", length = nrow(settings))
n_iter <- 10L
simulation_seeds <- vector(mode = "list", length = nrow(settings))

partial_derivatives_array <- array(NA, dim = c(length(grid_points), 2, 2, N_sim_total))
beta_array <- array(NA, dim = c(length(grid_points), 3L, 2L, N_sim_total))
# Loop through simulation replicates: -------------------------------------
for(i in seq_len(N_sim_total)) {
  print('----------------------------------------')
  print(paste0("Simulation rep ", i))
  simulation_seeds[[i]] <- .Random.seed
  mu <- settings[i, "mu"]
  dataset_i <- generate_vdp_curves_stochastic_ic_and_forcing(n = n,
                                                             mu = mu,
                                                             mu_ic = mu_ic,
                                                             grid_points = grid_points,
                                                             Sigma_ic = Sigma_ic,
                                                             sigma_forcing = settings$sigma_forcing[i],
                                                             intensity = settings$intensity_forcing[i])
  x_i <- dataset_i$data[,1,]
  y_i <- dataset_i$data[,2,]
  eps_x_i <- dataset_i$forcing[,1,]
  eps_y_i <- dataset_i$forcing[,2,]
  dx_i <- mu * (x_i - (1/3) * (x_i^3) - y_i) + eps_x_i
  dy_i <- (1 / mu) * x_i + eps_y_i
  x_mean <- apply(x_i, 1, mean)
  y_mean <- apply(y_i, 1, mean)
  
  for(tind in seq_along(grid_points)) {
    # PDA/ regression part:
    lm_xtind <- lm(dx_i[tind, ] ~ x_i[tind, ] + y_i[tind, ])
    lm_ytind <- lm(dy_i[tind, ] ~ x_i[tind, ] + y_i[tind, ])
    
    beta_array[tind,,1,i] <- coef(lm_xtind)
    beta_array[tind,,2,i] <- coef(lm_ytind)
  }  

  # Non-parametric fit of Dx: -----------------------------------------------
  # Do pffr fit: ------------------------------------------------------------
  pffr_fit <- pffr(formula = dx_i ~ - 1 + c(s(x_i, y_i, k = 120)),
                   data = list(
                     dx_i = t(dx_i),
                     x_i = t(x_i),
                     y_i = t(y_i)),
                   yind = grid_points,
                   algorithm = "bam",
                   method = "fREML",
                   bs.yindex = list(bs="ps", k = 120, m=c(2, 1)))
  
  F_grid <- coefficients(object = pffr_fit,
                         se = FALSE, 
                         n1 = 200,
                         n2 = 200)$smterms$`c(s(x_i,y_i,120))`$coef
  
  surface_mat <- as.matrix(dcast.data.table(
    data = as.data.table(F_grid),
    value.var = "value", 
    formula = y_i ~ x_i))[, -1]
  rownames(surface_mat) <- unique(F_grid$y_i)
  # Calculate Partial Derivatives via Finite Differencing:
  # Get values and time-steps in xh and xk directions.
  # These will be the same for surface 
  (x_inds <- as.numeric(colnames(surface_mat)))
  (y_inds <- as.numeric(rownames(surface_mat)))
  
  stopifnot(length(unique(round(diff(x_inds), 8))) == 1)
  stopifnot(length(unique(round(diff(y_inds), 8))) == 1)
  
  # Time Steps (needed for denominator in differencing:)
  (delta_x <- unique(round(diff(x_inds), 8)))
  (delta_y <- unique(round(diff(y_inds), 8)))
  
  # Call F_H(x^h, x^k) "FH" i.e., our fitted/ predicted 
  # surface for dependent variable D2xH
  F_surf <- surface_mat
  # Set up matrices to store partial derivative surfaces:
  dFdx <- dFdy <- matrix(NA, nrow = nrow(F_surf), ncol = ncol(F_surf))
  
  # OK calculate forward finite-difference approximation for dFHdxh
  # Loop through rows (hold xk constant)
  # Difference columns (differencing wrt xh)
  for(i2 in seq_len(nrow(F_surf))) {
    dFdx[i2, ] <- c(NA, diff(F_surf[i2, ])) / delta_x
  }
  
  # Now do the opposite for xk
  for(i2 in seq_len(ncol(F_surf))) {
    dFdy[, i2] <- c(NA, diff(F_surf[, i2])) / delta_y
  }
  
  # get the indices of the surface matrix that these corrsespond to (i and j are rows and columns)
  i_location <- j_location <- inds <- vector(mode = "numeric", length  = length(x_mean))
  for(i2 in seq_len(length(x_mean))) {
    inds[i2] <- which.min((F_grid[, 1] - x_mean[i2]) ^ 2 + (F_grid[, 2] - y_mean[i2]) ^ 2) 
    i_location[i2] <- which.min(abs(round(F_grid$y[inds[i2]] - y_inds, 12)))
    j_location[i2] <- which.min(abs((round(F_grid$x[inds[i2]] - x_inds, 12))))
  }
  
  

  # K's mean different things... confusing
  for(k in seq_along(i_location)) {
    partial_derivatives_array[k,1,1,i] <- dFdx[i_location[k], j_location[k]]
    partial_derivatives_array[k,2,1,i] <- dFdy[i_location[k], j_location[k]]
  }
  
  

  # -------------------------------------------------------------------------
  # Do pffr fit: ------------------------------------------------------------
  pffr_fit <- pffr(formula = dy_i ~ - 1 + c(s(x_i, y_i, k = 120)),
                   data = list(
                     dy_i = t(dy_i),
                     x_i = t(x_i),
                     y_i = t(y_i)),
                   yind = grid_points,
                   algorithm = "bam",
                   method = "fREML",
                   bs.yindex = list(bs="ps", k = 120, m=c(2, 1)))
  
  F_grid <- coefficients(object = pffr_fit,
                         se = FALSE, 
                         n1 = 200,
                         n2 = 200)$smterms$`c(s(x_i,y_i,120))`$coef
  
  surface_mat <- as.matrix(dcast.data.table(
    data = as.data.table(F_grid),
    value.var = "value", 
    formula = y_i ~ x_i))[, -1]
  rownames(surface_mat) <- unique(F_grid$y_i)
  # Calculate Partial Derivatives via Finite Differencing:
  # Get values and time-steps in xh and xk directions.
  # These will be the same for surface 
  (x_inds <- as.numeric(colnames(surface_mat)))
  (y_inds <- as.numeric(rownames(surface_mat)))
  
  stopifnot(length(unique(round(diff(x_inds), 8))) == 1)
  stopifnot(length(unique(round(diff(y_inds), 8))) == 1)
  
  # Time Steps (needed for denominator in differencing:)
  (delta_x <- unique(round(diff(x_inds), 8)))
  (delta_y <- unique(round(diff(y_inds), 8)))
  
  # Call F_H(x^h, x^k) "FH" i.e., our fitted/ predicted 
  # surface for dependent variable D2xH
  F_surf <- surface_mat
  # Set up matrices to store partial derivative surfaces:
  dFdx <- dFdy <- matrix(NA, nrow = nrow(F_surf), ncol = ncol(F_surf))
  
  # OK calculate forward finite-difference approximation for dFHdxh
  # Loop through rows (hold xk constant)
  # Difference columns (differencing wrt xh)
  for(i2 in seq_len(nrow(F_surf))) {
    dFdx[i2, ] <- c(NA, diff(F_surf[i2, ])) / delta_x
  }
  
  # Now do the opposite for xk
  for(i2 in seq_len(ncol(F_surf))) {
    dFdy[, i2] <- c(NA, diff(F_surf[, i2])) / delta_y
  }
  
  # get the indices of the surface matrix that these corrsespond to (i and j are rows and columns)
  i_location <- j_location <- inds <- vector(mode = "numeric", length  = length(x_mean))
  for(i2 in seq_len(length(x_mean))) {
    inds[i2] <- which.min((F_grid[, 1] - x_mean[i2]) ^ 2 + (F_grid[, 2] - y_mean[i2]) ^ 2) 
    i_location[i2] <- which.min(abs(round(F_grid$y[inds[i2]] - y_inds, 12)))
    j_location[i2] <- which.min(abs((round(F_grid$x[inds[i2]] - x_inds, 12))))
  }
  

  # K's mean different things... confusing
  for(k in seq_along(i_location)) {
    partial_derivatives_array[k,1,2,i] <- dFdx[i_location[k], j_location[k]]
    partial_derivatives_array[k,2,2,i] <- dFdy[i_location[k], j_location[k]]
  }
  }


saveRDS(object = 
          list(partial_derivatives_array = partial_derivatives_array,
               beta_array = beta_array,
               simulation_seeds = simulation_seeds,
               session_info = sessionInfo()),
        file = here::here(
          "outputs",
          "VDP",
          "vdp-simulation-partial_derivatives_note.rds"))


