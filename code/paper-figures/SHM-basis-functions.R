# Packages: ---------------------------------------------------------------
library(data.table) # CRAN v1.14.2
library(tikzDevice) # CRAN v0.12.3.1
library(ggplot2)    # CRAN v3.3.5

# Results of estimated covariance functions: ------------------------------
results_list <- readRDS(file = here::here("outputs", "SHM",
                                          "simulation-results",
                                          "simulation-covariance-estimation.rds"))


covariance_est_list <- results_list[["covariance_est_list"]]
grid_range <- c(0, 2 * pi)
n_grid <- 50
grid_points <- seq(from = grid_range[1],
                   to = grid_range[2],
                   length.out = n_grid)


# Results of true covariance function -------------------------------------
# note that this was computed at 100 instead of 50 points:
cov_xs_xt <- readRDS(file = here::here("outputs",
                                       "SHM",
                                       "simulation-results",
                                       "cov_xs_xt_new.rds"))
grid_points_longer <- seq(from = grid_range[1],
                          to = grid_range[2],
                          length.out = 100)
test_grid_longer <- expand.grid(grid_points_longer, grid_points_longer)





# Extracting zero0state covariance estimates: -----------------------------
# We only computed upper-diagonal + diagonal elements to save time because 
# the covariance function is symmetric.
test_grid <- rbind(
  cbind(grid_points, grid_points),
  t(combn(grid_points, 2))
)
n_sim <- 10
covariance_results <- vector("list", n_sim) # to store covariance matrices initially

for(k in seq_len(n_sim)) {
  res1_dt <- data.table(cbind(test_grid, sapply(covariance_est_list[[k]], function(x){x[1,1]})))
  names(res1_dt) <- c("s", "t", "cov_xs_xt")
  cov_est_k_mat <- matrix(data = NA, ncol = n_grid, nrow = n_grid)
  
  for(ind in seq_len(nrow(res1_dt))) {
    sind <- res1_dt[ind, s]
    tind <- res1_dt[ind, t]
    
    i <- which(grid_points == sind)
    j <- which(grid_points == tind)
    
    cov_est_k_mat[i, j] <- cov_est_k_mat[j, i] <- res1_dt[ind, cov_xs_xt]
  }
  
  covariance_results[[k]] <- cov_est_k_mat
  
}


# Compute Eigenfunctions of estimated zero-state covariance ---------------
names(covariance_results) <- seq_len(10)
zero_state_efun_dt <- 
  data.table(
    purrr::map_dfr(.x = covariance_results, .f = function(x) {
      df <- data.frame(grid_points, eigen(x)$vectors[, 1:4])
      names(df) <- c("t", paste0("efun_", 1:4))
      df}, .id = "sim_rep"))




# Compute estimated canonical PDA basis functions: ------------------------
beta_est_mat <- results_list$results_array[,4,]
zero_input_efuns_list <- vector("list", length = n_sim)
names(zero_input_efuns_list) <- seq_len(n_sim)
for(k in seq_len(n_sim)) {
  print(paste0("Iteration", k))
  beta_k_t <- approxfun(grid_points, beta_est_mat[, k])
  dynamics_equations <- function(t, y, ...) {
    with(as.list(c(y)),{
      # rate of change
      dX <- Y
      dY <-  beta_k_t(t) * X  
      list(c(dX, dY))
    }
    )
    
  }
  
  tmp_df <- data.frame(
    cbind(deSolve::lsoda(y = c(X = 1, Y = 0),
                                 times = grid_points,
                                 func = dynamics_equations,
                                 tcrit = grid_range[2])[, c("time", "X")],
                  deSolve::lsoda(y = c(X = 0, Y = 1),
                                 times = grid_points,
                                 func = dynamics_equations,
                                 tcrit = grid_range[2])[, c("X")]))
  
  names(tmp_df) <- c("t", "cos_t", "sin_t")
  zero_input_efuns_list[[k]] <- tmp_df
}

zero_input_efun_dt <- data.table(
  purrr::map_dfr(zero_input_efuns_list, .f = function(x) {x}, .id = "sim_rep")
)

efun_dt <- merge.data.table(x = zero_state_efun_dt,
                            y =  zero_input_efun_dt,
                            by = c("t", "sim_rep"))
efun_dt_lng <- melt.data.table(data = efun_dt,
                               id.vars = c("t", "sim_rep"), 
                               variable.name = "efun_name", 
                               value.name = "phi",
                               verbose = TRUE,
                               variable.factor = FALSE,
                               value.factor = FALSE, measure.vars = c(
                                 "cos_t", "sin_t", paste0("efun_", 1:4)
                               ))

efun_dt_lng[, covariance_type := fifelse(efun_name %in% c("cos_t", "sin_t"), yes = "(c) Zero-input basis functions",
                                         no = "(d) Zero-state basis functions")]

efun_dt_lng[, efun_name := factor(efun_name,
                                  levels = c(
                                    "cos_t",
                                    "sin_t",
                                    paste0("efun_", 1:4)
                                  ), 
                                  labels = c(
                                    "$\\cos(t)$",
                                    "$\\sin(t)$",
                                    paste0("Zero-state Basis Function $", 1:4, "$")
                                  ))]








# Compute TRUE eiegenfunctions of zewro-state covariance ------------------
zero_state_covariance_df <- cbind(test_grid_longer, sapply(cov_xs_xt, function(x) {x[1,1]}))
names(zero_state_covariance_df) <- c("s", "t", "cov")
zero_state_covariance_df_wide <- reshape2::dcast(data = zero_state_covariance_df,
                                                 formula = s ~ t,
                                                 value.var = "cov")
stopifnot(zero_state_covariance_df_wide$s == grid_points_longer)
stopifnot(colnames(zero_state_covariance_df_wide)[-1] ==grid_points_longer)
zero_state_covariance_matrix_wide <- as.matrix(zero_state_covariance_df_wide[, -1])

zero_state_covariance_matrix_sym <- (zero_state_covariance_matrix_wide + t(zero_state_covariance_matrix_wide)) / 2
zero_state_eigen <- eigen(zero_state_covariance_matrix_sym)
# and the true canonical basis functions are just sin(t) and cosine(t).
efun_dt_true <- data.table(
  t = grid_points_longer, 
  sim_rep = "truth",
  cos_t = cos(grid_points_longer),
  sin_t = sin(grid_points_longer),
  zero_state_eigen$vectors[, 1:4]  * sqrt(2) # re-scaled because these are evuated at twice as many points as estimates
)
names(efun_dt_true)[-c(1:4)] <- paste0("efun_", 1:4)
efun_dt_true_lng <-  melt.data.table(data = efun_dt_true,
                                     id.vars = c("t", "sim_rep"), 
                                     variable.name = "efun_name", 
                                     value.name = "phi",
                                     verbose = TRUE,
                                     variable.factor = FALSE,
                                     value.factor = FALSE, measure.vars = c(
                                       "cos_t", "sin_t", paste0("efun_", 1:4)
                                     ))
efun_dt_true_lng[, covariance_type := fifelse(efun_name %in% c("cos_t", "sin_t"), yes = "(c) Zero-input basis functions",
                                         no = "(d) Zero-state basis functions")]

efun_dt_true_lng[, efun_name := factor(efun_name,
                                  levels = c(
                                    "cos_t",
                                    "sin_t",
                                    paste0("efun_", 1:4)
                                  ), 
                                  labels = c(
                                    "$\\cos(t)$",
                                    "$\\sin(t)$",
                                    paste0("Zero-state Basis Function $", 1:4, "$")
                                  ))]

# ------------------------------------------------------------------------#
# ------------------------------------------------------------------------#
# ------------------------------------------------------------------------#
# ------------------------------------------------------------------------#
# ------------------------------------------------------------------------#
# ------------------------------------------------------------------------#




                    # Create plot: ------------------------------------------------------------
theme_set(theme_bw())
p <- ggplot(data = efun_dt_lng) +
  facet_wrap(~ covariance_type, scales = "free_y") +
  aes(x = t, y = phi, group = interaction(efun_name, sim_rep), colour = efun_name) +
  geom_line(alpha = 0.5) +
  # geom_line(data = efun_dt_true_lng, col=1)+
  geom_line(data = efun_dt_true_lng, lwd = 0.75, lty =2) +
  scale_x_continuous(breaks = c(0, pi, 2 * pi), 
                     labels = c("$0$", "$\\pi$", "$2 \\pi$"),
                     limits = c(0, 2 * pi), expand = rep(0.025, 2)) +
  labs(colour = "Function", x = "$t$",
       y = "Basis Function Value") +
  theme(axis.title = element_text(size = 9),
        legend.position = "bottom",
        strip.text = element_text(face="bold"),
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        panel.spacing = unit(2, "lines"),
        legend.text = element_text(size = 8.5),
        plot.margin = margin(t = 10, r = 15, b = 10, l = 10),
        legend.box.margin = margin(t = -10),
        plot.title = element_text(size = 10, hjust = 0.5, face = "bold"))
p


# And save: ---------------------------------------------------------------
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937
tikz(here::here("outputs", "SHM", "paper-plots", "basis-funs.tex"), 
     width = 0.9 * (0.9 * doc_width_inches),
     height = 0.9 * (1.05 * doc_width_inches)/2,
     standAlone = TRUE)
print(p)              
dev.off()

tinytex::pdflatex(here::here("outputs", "SHM", "paper-plots", "basis-funs.tex"))


