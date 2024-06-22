pda_results <- readRDS(file = here::here("outputs", "real-data", "pda-result-01.rds"))
prepared_data <- readRDS(here::here("outputs", "real-data", "prepared-data-01.rds"))
x <- prepared_data$x
Dx <- prepared_data$Dx
D2x <- prepared_data$D2x
N <- prepared_data$N
resid <- pda_results$pda_result$resid
milliseconds <- prepared_data$milliseconds
beta_hat <- pda_results$pda_result$beta
source("code/simulate_from_model_funs.R")
source("PDA-2-pw-general.R")

set.seed(1)
simulated_dataset_final <- generate_com_dataset_general(n = N, 
                             b0_fun = approxfun(milliseconds, beta_hat[,1,11]),
                             beta_0_fun = approxfun(milliseconds, beta_hat[,2,11]),
                             beta_1_fun = approxfun(milliseconds, beta_hat[,3,11]), 
                             mu_init_vec = c(mean(x[1,]), mean(Dx[1,])),
                             sigma_init_mat = cov(cbind(x[1,], Dx[1,])), 
                             C_grid = var(t(resid$final)),
                             grid_points = milliseconds)

simulated_dataset_init <- generate_com_dataset_general(n = N, 
                                                        b0_fun = approxfun(milliseconds, beta_hat[,1,1]),
                                                        beta_0_fun = approxfun(milliseconds, beta_hat[,2,1]),
                                                        beta_1_fun = approxfun(milliseconds, beta_hat[,3,1]), 
                                                        mu_init_vec = c(mean(x[1,]), mean(Dx[1,])),
                                                        sigma_init_mat = cov(cbind(x[1,], Dx[1,])), 
                                                        C_grid = var(t(resid$initial)),
                                                        grid_points = milliseconds)

par(mfrow = c(3, 2))
matplot(milliseconds, x, type = "l", col = scales::alpha(1, 0.5), lty = 1)
matplot(milliseconds, Dx, type = "l", col = scales::alpha(1, 0.5), lty = 1)

matplot(milliseconds, simulated_dataset_final$x, type = "l", col = scales::alpha(1, 0.5), lty = 1)
matplot(milliseconds, simulated_dataset_final$dx, type = "l", col = scales::alpha(1, 0.5), lty = 1)

matplot(milliseconds, simulated_dataset_init$x, type = "l", col = scales::alpha(1, 0.5), lty = 1)
matplot(milliseconds, simulated_dataset_init$dx, type = "l", col = scales::alpha(1, 0.5), lty = 1)

par(mfrow = c(1,1))
plot(apply(x, 1, sd),type = "l", ylim = c(0.002, 0.01))
lines(apply(simulated_dataset_final$x, 1, sd), col=2)
lines(apply(simulated_dataset_init$x, 1, sd), col=4)

# Do simulation: ----------------------------------------------------------
set.seed(1)
# some settings
num_iter <- 15L
k <- 50
resid_pve <- 0.99
n_cores <- parallel::detectCores()
N_sim <- 15
pda_simulation_results_list <- seeds_list <- vector("list", length = N_sim)
x_array <- Dx_array <- D2x_array <- array(data = NA, dim = c(length(milliseconds), N, N_sim))

for(i in seq_len(N_sim)) {
  cat("---------------------------------------------------------------------\n")
  print(paste("Iteration", i))
  cat("---------------------------------------------------------------------\n")
  
  seeds_list[[i]] <- .Random.seed
  
  print("Generating data")
  dataset_i <- generate_com_dataset_general(
    n = N, 
    b0_fun = approxfun(milliseconds, beta_hat[,1,11]),
    beta_0_fun = approxfun(milliseconds, beta_hat[,2,11]),
    beta_1_fun = approxfun(milliseconds, beta_hat[,3,11]), 
    mu_init_vec = c(mean(x[1,]), mean(Dx[1,])),
    sigma_init_mat = cov(cbind(x[1,], Dx[1,])), 
    C_grid = var(t(resid$final)),
    grid_points = milliseconds)
  
  x_array[,,i] <- x_i <- dataset_i$x
  fd_obj_i <- fda::Data2fd(argvals = milliseconds, 
                           y = x_i, 
                           basisobj = prepared_data$quintic_Bspline_basis,
                           lambda = 10^-12)
  Dx_array[,,i] <- Dx_i <- fda::eval.fd(milliseconds, fd_obj_i, Lfdobj = 1)
  D2x_array[,,i] <- D2x_i <- fda::eval.fd(milliseconds, fd_obj_i, Lfdobj = 2)
  
  # Do PDA: -----------------------------------------------------------------
  
    pda_simulation_results_list[[i]] <- try(
      do_pda_iterative_br_parallel_post_smooth(
        grid_points = milliseconds,
        x = x_i,
        Dx = Dx_i,
        D2x = D2x_i,
        n = N,
        num_iterations = num_iter, 
        verbose = TRUE,
        n_cores = n_cores,
        k = k,
        method = "REML", 
        resid_smooth = TRUE, 
        resid_k = k, 
        resid_pve = resid_pve))
  
  }

saveRDS(object = list(pda_simulation_results_list = pda_simulation_results_list,
                      x_array = x_array,
                      Dx_array = Dx_array,
                      D2x_array = D2x_array),
        file = here::here("outputs", "real-data", "pda-simulation-result-01.rds"))
