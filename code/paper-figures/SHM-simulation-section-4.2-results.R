source(here::here("code", "SHM-bias-calculation-function.R"))
source(here::here("code", "SHM-helper-functions.R"))
library(data.table) # Extension of `data.frame`
library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics
library(tikzDevice) # R Graphics Output in LaTeX Format
options(tikzLatexPackages = c(getOption( "tikzLatexPackages" ),"\\usepackage{amsfonts}"))
options(tikzLatexPackages = c(getOption( "tikzLatexPackages" ),"\\usepackage{amsmath}"))
options(tikzLatexPackages = c(getOption( "tikzLatexPackages" ),"\\usepackage{amssymb}"))


# simulation_result_list <- readRDS(file = here::here("outputs",
#                                                     "SHM",
#                                                     "simulation-results",
#                                                     "simulation-section-4.2.rds"))
# results_array <- simulation_result_list$results_array
simulation_result_list <- readRDS(here::here("outputs",
                                             "SHM",
                                             "simulation-results",
                                             "simulation-section-4.2-updated.rds"))
n_grid <- 100
n_sim <- 500
results_array <- array(unlist(simulation_result_list$results_list), dim = c(n_grid, 4, n_sim))


# -------------------------------------------------------------------------

grid_range <- c(0, 2 * pi)
n_grid <- 100
grid_points <- seq(from = grid_range[1],
                   to = grid_range[2],
                   length.out = n_grid)


K_st_true <- function(s, t) {
  0.25 ^ 2 * dnorm(x = 2 * abs(t - s))
}

beta_u_true <- function(u) {
  - 1 * (u / u)
}

# true_bias_numerator <- bias_calculation(beta_u_est = beta_u_true, K_st_est = K_st_true, grid_points = grid_points)

B_bar_s_true_fun <- function(s) {
  matrix(data = c(0,  -s, s, 0), 
         nrow = 2, 
         ncol = 2, 
         byrow = FALSE)
}


# INTEGRAND ---------------------------------------------------------------
integrand_uv_11 <- function(u, v) {
  (Matrix::expm(x = - B_bar_s_true_fun(u)) %*%
     matrix(c(0, 0, 0, K_st_true(u, v)), ncol = 2, nrow = 2) %*%
     t(as.matrix(Matrix::expm(x = - B_bar_s_true_fun(v)))))[1, 1]
}
integrand_uv_12 <- function(u, v) {
  (Matrix::expm(x = - B_bar_s_true_fun(u)) %*%
     matrix(c(0, 0, 0, K_st_true(u, v)), ncol = 2, nrow = 2) %*%
     t(as.matrix(Matrix::expm(x = - B_bar_s_true_fun(v)))))[1, 2]
}
integrand_uv_21 <- function(u, v) {
  (Matrix::expm(x = - B_bar_s_true_fun(u)) %*%
     matrix(c(0, 0, 0, K_st_true(u, v)), ncol = 2, nrow = 2) %*%
     t(as.matrix(Matrix::expm(x = - B_bar_s_true_fun(v)))))[2, 1]
}
integrand_uv_22 <- function(u, v) {
  (Matrix::expm(x = - B_bar_s_true_fun(u)) %*%
     matrix(c(0, 0, 0, K_st_true(u, v)), ncol = 2, nrow = 2) %*%
     t(as.matrix(Matrix::expm(x = - B_bar_s_true_fun(v)))))[2, 2]
}



# VECTORISE INTEGRANDS ----------------------------------------------------
# Over u:
integrand_uv_11_vectorise_u <- Vectorize(integrand_uv_11, vectorize.args = "u")
integrand_uv_12_vectorise_u <- Vectorize(integrand_uv_12, vectorize.args = "u")
integrand_uv_21_vectorise_u <- Vectorize(integrand_uv_21, vectorize.args = "u")
integrand_uv_22_vectorise_u <- Vectorize(integrand_uv_22, vectorize.args = "u")


# PART INTEGRATE integrands -----------------------------------------------
# over u
part_integral_vt_11 <- function(v, t) {
  integrate(f = integrand_uv_11_vectorise_u, lower = 0, upper = t, v = v)$value
}
part_integral_vt_12 <- function(v, t) {
  integrate(f = integrand_uv_12_vectorise_u, lower = 0, upper = t, v = v)$value
}
part_integral_vt_21 <- function(v, t) {
  integrate(f = integrand_uv_21_vectorise_u, lower = 0, upper = t, v = v)$value
}
part_integral_vt_22 <- function(v, t) {
  integrate(f = integrand_uv_22_vectorise_u, lower = 0, upper = t, v = v)$value
}

# VECTORISE PART INTEGRALS ------------------------------------------------
# over v
part_integral_vt_11_vectorise_v <- Vectorize(part_integral_vt_11, vectorize.args = "v")
part_integral_vt_12_vectorise_v <- Vectorize(part_integral_vt_12, vectorize.args = "v")
part_integral_vt_21_vectorise_v <- Vectorize(part_integral_vt_21, vectorize.args = "v")
part_integral_vt_22_vectorise_v <- Vectorize(part_integral_vt_22, vectorize.args = "v")


full_integral_st_matrix <- function(s, t) {
  matrix(data = c(
    integrate(part_integral_vt_11_vectorise_v, lower = 0, upper = s, t = t)$value,
    integrate(part_integral_vt_12_vectorise_v, lower = 0, upper = s, t = t)$value,
    integrate(part_integral_vt_21_vectorise_v, lower = 0, upper = s, t = t)$value,
    integrate(part_integral_vt_22_vectorise_v, lower = 0, upper = s, t = t)$value
  ), nrow = 2, ncol = 2, byrow = TRUE)
}


full_covariance_st <- function(s, t) {
  as.matrix(Matrix::expm(x = B_bar_s_true_fun(s = s))) %*%
    full_integral_st_matrix(s = s, t = t) %*%
    t(as.matrix(Matrix::expm(x = B_bar_s_true_fun(s = t))))
}

true_bias_denom_term_2 <- sapply(grid_points, function(ti) {
  full_covariance_st(s = ti, t = ti)
})


# Calculate the Bias (Truth) ------------------------------------------------------
# K_st_true <- function(s, t) {
#   0.25 ^ 2 * dnorm(x = 2 * abs(t - s))
# }
# 
# beta_u_true <- function(u) {
#   - 1 * (u / u)
# }
# 
# B_bar_s_true_fun <- function(s) {
#   matrix(data = c(0,  -s, s, 0), 
#          nrow = 2, 
#          ncol = 2, 
#          byrow = FALSE)
# }
# integrand_uv <- function(u, v) {
#   (Matrix::expm(x = - B_bar_s_true_fun(u)) %*%
#      matrix(c(0, 0, 0, K_st_true(u, v)), ncol = 2, nrow = 2) %*%
#      t(as.matrix(Matrix::expm(x = - B_bar_s_true_fun(v)))))[1, 1]
# }
# integrand_uv_vectorise_u <- Vectorize(integrand_uv, vectorize.args = "u")
# part_integral_vt <- function(v, t) {
#   integrate(f = integrand_uv_vectorise_u, lower = 0, upper = t, v = v)$value
# }
# part_integral_vt_vectorise_v <- Vectorize(part_integral_vt, vectorize.args = "v")
# full_integral_t <- function(t) {
#   integrate(part_integral_vt_vectorise_v, lower = 0, upper = t, t = t)$value
# }
# full_integral_t_vectorise <- Vectorize(FUN = full_integral_t, vectorize.args = "t")
true_bias_denominator <- 0.05 + true_bias_denom_term_2[1, ]
true_bias_numerator <- bias_calculation(beta_u_est = beta_u_true,
                                        K_st_est = K_st_true, 
                                        grid_points = grid_points)
true_bias <- true_bias_numerator[1,2,] / true_bias_denominator

# par(mfrow = c(1, 2))
# plot(grid_points, true_bias_numerator[1,2,], type = "l", 
#      xlab = "t", ylab = expression("E[x(t)"~epsilon~"(t)]"))
# title("Numerator")
# plot(grid_points, true_bias_denominator, type = "l", 
#      xlab = "t", ylab = expression("E["~x(t)^2~"]"))
# title("Denominator")
# Calculated estimated bias from simulation -------------------------------

bias_result <- apply(results_array + 1, c(1, 2), mean)

# Create Datasets for Plotting --------------------------------------------

bias_beta_dt <- data.table(t_grid = grid_points, results_array[,1,])
bias_beta_dt_long <- melt.data.table(data = bias_beta_dt,
                                     id.vars = "t_grid",
                                     variable.name = "rep_fun",
                                     variable.factor = TRUE,
                                     value.name = "beta", 
                                     value.factor = FALSE,
                                     verbose = TRUE)
bias_dt <- data.table(t_grid = grid_points, bias_est = bias_result[,1], bias_true = true_bias)
bias_dt_long <- melt.data.table(data = bias_dt,
                                id.vars = "t_grid",
                                variable.name = "calc",
                                variable.factor = TRUE,
                                value.name = "bias", 
                                value.factor = FALSE,
                                verbose = TRUE)
bias_dt_long[, calc := fcase(calc ==  "bias_true", "Truth", calc == "bias_est", "Simulated")]


# Create Plots ------------------------------------------------------------

beta_plot <- ggplot(data = bias_beta_dt_long) +
  aes(x = t_grid, y = beta, group = rep_fun) +
  geom_path(colour= "darkgrey", alpha=0.5) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = c(0, pi, 2 * pi), labels = c("$0$", "$\\pi$", "$2 \\pi$"),
                     limits = c(0, 2 * pi), expand = rep(0.025, 2)) +
  labs(x = "$t$", 
       y = "$\\widehat{\\beta}_0 (t)$", 
       title = "Estimates")+#: $\\widehat{\\beta}_0 (t)$") +
  geom_hline(yintercept = -1, linetype = "dotted", col = "black", size = 0.75) +
  theme(axis.title = element_text(size = 8),
        plot.title = element_text(size = 9, hjust = 0.5, face = "bold"))


bias_plot <- ggplot(data = bias_dt_long) +
  aes(x = t_grid, y = bias, group = calc, colour = calc, linetype = calc) +
  geom_path() +
  theme_bw() +
  theme(legend.position = c(0.5, 0.3),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(size = 8)) +
  scale_x_continuous(breaks = c(0, pi, 2 * pi), labels = c("$0$", "$\\pi$", "$2 \\pi$"),
                     limits = c(0, 2 * pi), expand = rep(0.025, 2)) +
  labs(x = "$t$",
       y = "$\\mathbb{E}[\\widehat{\\beta}_0 (t)] - \\beta_0 (t)$",
       title = "Bias")+ # : $E[\\widehat{\\beta}_0 (t)] - \\beta_0 (t)$") +
  # geom_hline(yintercept = 0, linetype = "dotted", col = "black", size = 0.75) +
  theme(axis.title = element_text(size = 8),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.title = element_text(size = 9, hjust = 0.5, face = "bold"))

full_plot <- ggpubr::ggarrange(plotlist = list(beta_plot, bias_plot),
                               nrow = 1,
                               ncol = 2, 
                               labels = c("(a)", "(b)"), 
                               align = "hv", 
                               font.label = list(size = 8))
# full_plot <- ggpubr::ggarrange(
#   plotlist = list(beta_plot, NULL,  bias_plot),
#   nrow = 1,
#   ncol = 3,
#   widths = c(1, 0.05, 1),
#   labels = c("(a)", "", "(b)"))


full_plot
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937
tikz(here::here("outputs", "SHM", "paper-plots", "bias-beta.tex"), 
     width = (0.9 * doc_width_inches),
     height = (0.8 * doc_width_inches)/2)
print(full_plot)
dev.off()

