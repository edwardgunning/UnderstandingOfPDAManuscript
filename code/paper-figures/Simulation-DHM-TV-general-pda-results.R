# ------------------------------------------------------------------------#
# Simple simulation to demonstrate general proposed bias correction on
# SHM model.
# ------------------------------------------------------------------------#

# Load necessary functions and packages: ----------------------------------
source(here::here("code", "SHM-helper-functions.R"))
source(here::here("code", "DHM-TV-general-data-generation.R"))
library(fda) # Functional Data Analysis, CRAN v5.5.1
library(data.table) # CRAN v1.14.2
library(ggplot2)    # CRAN v3.4.0
library(tikzDevice) # CRAN v0.12.3.1


# Plot settings: ----------------------------------------------------------

# Plot settings -----------------------------------------------------------
source(here::here("code", "theme_gunning.R"))
theme_gunning()
theme_update(legend.position = "bottom", 
             legend.title = element_text(face = "bold", size = 13),
             strip.text = element_text(size = 12),
             axis.text = element_text(size = 11),
             axis.title = element_text(size = 12),
             legend.text = element_text(size = 13))
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937


# Path to store results: --------------------------------------------------
results_path <- here::here("outputs", "SHM", "simulation-results", "general-DHM-TV-simulation.rds")
results <- readRDS(file = results_path)

grid_points <- results$grid_points
sigma <- results$sigma
mu_init <- results$mu_init
intensity <- results$intensity
sigma_init <- results$sigma_init
n <- results$n_i
betas <- results$beta
num_simulations <- results$num_simulations
simulation_seeds <- results$simulation_seeds
damping_fun_test <- function(t) {
  0.01 *  (t - 2 * pi) ^2
}
.Random.seed <- simulation_seeds[[1]]
dataset_1 <- generate_DHM_TV_dataset_general(n = n,
                                          grid_points = grid_points,
                                          sigma = 0.4,
                                          mu_init = c(1, 0),
                                          intensity = 1,
                                          sigma_init = diag(rep(0.05, 2)),
                                          damping_fun = damping_fun_test)

x_1 <- dataset_1$x
Dx_1 <- dataset_1$dx
eps_1 <- dataset_1$noise
D2x_1 <- - x_1 + eps_1

# Convert objects to data.table for plotting:

# -------------------------------------------------------------------------


## x ----------------------------------------------------------------------
x_1_dt <- data.table(
  grid_points, x_1
)
names(x_1_dt) <- c("t", paste0("obs_", seq_len(n)))
x_1_dt_lng <- melt.data.table(data = x_1_dt,
                              id.vars = "t",
                              value.factor = FALSE,
                              value.name = "x",
                              variable.name = "obs",
                              variable.factor = FALSE,
                              verbose = TRUE, 
                              measure.vars = paste0("obs_", seq_len(n)))

## Dx ---------------------------------------------------------------------
Dx_1_dt <- data.table(
  grid_points, Dx_1
)
names(Dx_1_dt) <- c("t", paste0("obs_", seq_len(n)))
Dx_1_dt_lng <- melt.data.table(data = Dx_1_dt,
                               id.vars = "t",
                               value.factor = FALSE,
                               value.name = "Dx",
                               variable.name = "obs",
                               variable.factor = FALSE,
                               verbose = TRUE, 
                               measure.vars = paste0("obs_", seq_len(n)))

## D2x ---------------------------------------------------------------------
D2x_1_dt <- data.table(
  grid_points, D2x_1
)
names(D2x_1_dt) <- c("t", paste0("obs_", seq_len(n)))
D2x_1_dt_lng <- melt.data.table(data = D2x_1_dt,
                                id.vars = "t",
                                value.factor = FALSE,
                                value.name = "D2x",
                                variable.name = "obs",
                                variable.factor = FALSE,
                                verbose = TRUE, 
                                measure.vars = paste0("obs_", seq_len(n)))

## eps ---------------------------------------------------------------------
eps_1_dt <- data.table(
  grid_points, eps_1
)
names(eps_1_dt) <- c("t", paste0("obs_", seq_len(n)))
eps_1_dt_lng <- melt.data.table(data = eps_1_dt,
                                id.vars = "t",
                                value.factor = FALSE,
                                value.name = "eps",
                                variable.name = "obs",
                                variable.factor = FALSE,
                                verbose = TRUE, 
                                measure.vars = paste0("obs_", seq_len(n)))



# Create data frame for plotting: -----------------------------------------

x_dx_dt <- merge.data.table(x_1_dt_lng, y = Dx_1_dt_lng, by = c("t", "obs"), all = TRUE)
stopifnot(nrow(x_dx_dt) == nrow(Dx_1_dt_lng))
x_dx_d2x_dt <-  merge.data.table(x_dx_dt, y = D2x_1_dt_lng, by = c("t", "obs"), all = TRUE)
stopifnot(nrow(x_dx_dt) == nrow(x_dx_d2x_dt))
plot_dt_wide <- merge.data.table(x_dx_d2x_dt, eps_1_dt_lng, by = c("t", "obs"), all = TRUE)
stopifnot(nrow(plot_dt_wide) == nrow(x_dx_d2x_dt))


# Reshape: ----------------------------------------------------------------

plot_dt_lng <- melt.data.table(plot_dt_wide,
                               id.vars = c("t", "obs"),
                               variable.name = "fun", 
                               variable.factor = FALSE,
                               value.name = "value", 
                               value.factor = FALSE,
                               measure.vars = c("x", "Dx", "D2x", "eps"),
                               verbose = TRUE)
plot_dt_lng[, fun := factor(fun, 
                            levels = c("x", "Dx", "D2x", "eps"),
                            labels = c("$x(t)$","$Dx(t)$", "$D^2x(t)$", "$\\epsilon(t)$"))]


# Plot: -------------------------------------------------------------------

dataset_plot <- ggplot(plot_dt_lng) +
  aes(x = t, y = value, group = obs) +
  geom_line(alpha = 0.9) +
  geom_line(data = plot_dt_lng[obs == "obs_1"], color = "orange") +
  facet_wrap(~ fun) +
  labs(x = "$t$") +
  theme(legend.position = "none",
        axis.title.y = element_blank()) +
  scale_x_continuous(breaks = c(0, pi, 2 * pi, 3 * pi, 4 * pi), labels = c("$0$", "$\\pi$", "$2 \\pi$", "$3 \\pi$", "$4 \\pi$"),
                     limits = c(0, 4 * pi), expand = rep(0.025, 2))
dataset_plot
tikz(here::here("outputs", "SHM", "paper-plots", "general-dhm-tv-datasets.tex"),
     width = doc_width_inches, 
     height =  doc_width_inches)
dataset_plot
dev.off()  


# Now plot of results: ----------------------------------------------------

# Compare initial and bias-reduced estimates
matplot(betas[,1,1,], type = "l", col = 1)
matlines(betas[,1,11,], col = 2)

matplot(betas[,2,1,], type = "l", col = 1)
matlines(betas[,2,11,], col = 2)

matplot(betas[,3,1,], type = "l", col = 1)
matlines(betas[,3,11,], col = 2)




# Intercept Function: -----------------------------------------------------
b_0_ols_dt <- data.table(t = grid_points, 
                         betas[,1,1,])
names(b_0_ols_dt)[-1] <- paste0("rep_", seq_len(num_simulations))
b_0_br_dt <- data.table(t = grid_points, 
                        betas[,1,11,])
names(b_0_br_dt)[-1] <- paste0("rep_", seq_len(num_simulations))

b_0_ols_dt_lng <- melt.data.table(data = b_0_ols_dt, 
                                  id.vars = "t", 
                                  measure.vars = paste0("rep_", seq_len(num_simulations)), 
                                  variable.name = "rep",
                                  value.name = "ols", 
                                  variable.factor = FALSE, 
                                  value.factor = FALSE)
b_0_br_dt_lng <- melt.data.table(data = b_0_br_dt, 
                                 id.vars = "t", 
                                 measure.vars = paste0("rep_", seq_len(num_simulations)), 
                                 variable.name = "rep",
                                 value.name = "br", 
                                 variable.factor = FALSE, 
                                 value.factor = FALSE)

b_0_dt <- merge.data.table(x = b_0_ols_dt_lng,
                           y = b_0_br_dt_lng,
                           by = c("t", "rep"), 
                           all = TRUE)
stopifnot(nrow(b_0_ols_dt_lng) == nrow(b_0_dt))
b_0_dt_lng <- melt.data.table(data = b_0_dt, 
                              id.vars = c("t", "rep"), 
                              measure.vars = c("ols", "br"),
                              variable.name = "estimator", 
                              variable.factor = FALSE, 
                              value.factor = FALSE,
                              value.name = "b_0")



# Beta_0 ------------------------------------------------------------------

beta_0_ols_dt <- data.table(t = grid_points, 
                            betas[,2,1,])
names(beta_0_ols_dt)[-1] <- paste0("rep_", seq_len(num_simulations))
beta_0_br_dt <- data.table(t = grid_points, 
                           betas[,2,11,])
names(beta_0_br_dt)[-1] <- paste0("rep_", seq_len(num_simulations))

beta_0_ols_dt_lng <- melt.data.table(data = beta_0_ols_dt, 
                                     id.vars = "t", 
                                     measure.vars = paste0("rep_", seq_len(num_simulations)), 
                                     variable.name = "rep",
                                     value.name = "ols", 
                                     variable.factor = FALSE, 
                                     value.factor = FALSE)
beta_0_br_dt_lng <- melt.data.table(data = beta_0_br_dt, 
                                    id.vars = "t", 
                                    measure.vars = paste0("rep_", seq_len(num_simulations)), 
                                    variable.name = "rep",
                                    value.name = "br", 
                                    variable.factor = FALSE, 
                                    value.factor = FALSE)

beta_0_dt <- merge.data.table(x = beta_0_ols_dt_lng,
                              y = beta_0_br_dt_lng,
                              by = c("t", "rep"), 
                              all = TRUE)
stopifnot(nrow(beta_0_ols_dt_lng) == nrow(beta_0_dt))
beta_0_dt_lng <- melt.data.table(data = beta_0_dt, 
                                 id.vars = c("t", "rep"), 
                                 measure.vars = c("ols", "br"),
                                 variable.name = "estimator", 
                                 variable.factor = FALSE, 
                                 value.factor = FALSE,
                                 value.name = "beta_0")

# Beta_1 ------------------------------------------------------------------
beta_1_ols_dt <- data.table(t = grid_points, 
                            betas[,3,1,])
names(beta_1_ols_dt)[-1] <- paste0("rep_", seq_len(num_simulations))
beta_1_br_dt <- data.table(t = grid_points, 
                           betas[,3,11,])
names(beta_1_br_dt)[-1] <- paste0("rep_", seq_len(num_simulations))

beta_1_ols_dt_lng <- melt.data.table(data = beta_1_ols_dt, 
                                     id.vars = "t", 
                                     measure.vars = paste0("rep_", seq_len(num_simulations)), 
                                     variable.name = "rep",
                                     value.name = "ols", 
                                     variable.factor = FALSE, 
                                     value.factor = FALSE)
beta_1_br_dt_lng <- melt.data.table(data = beta_1_br_dt, 
                                    id.vars = "t", 
                                    measure.vars = paste0("rep_", seq_len(num_simulations)), 
                                    variable.name = "rep",
                                    value.name = "br", 
                                    variable.factor = FALSE, 
                                    value.factor = FALSE)

beta_1_dt <- merge.data.table(x = beta_1_ols_dt_lng,
                              y = beta_1_br_dt_lng,
                              by = c("t", "rep"), 
                              all = TRUE)
stopifnot(nrow(beta_1_ols_dt_lng) == nrow(beta_1_dt))
beta_1_dt_lng <- melt.data.table(data = beta_1_dt, 
                                 id.vars = c("t", "rep"), 
                                 measure.vars = c("ols", "br"),
                                 variable.name = "estimator", 
                                 variable.factor = FALSE, 
                                 value.factor = FALSE,
                                 value.name = "beta_1")


b_0_beta_0_dt <- merge.data.table(x = b_0_dt_lng, y = beta_0_dt_lng,
                                  by = c("t", "rep", "estimator"),
                                  all = TRUE)
stopifnot(nrow(b_0_beta_0_dt) == nrow(b_0_dt_lng))
stopifnot(nrow(b_0_beta_0_dt) == nrow(beta_0_dt_lng))

stopifnot(nrow(b_0_beta_0_dt) == nrow(beta_1_dt_lng))
plot_params_dt_wide <- merge.data.table(x = b_0_beta_0_dt, 
                                        y = beta_1_dt_lng,
                                        all = TRUE,
                                        by = c("t", "rep", "estimator"))
stopifnot(nrow(plot_params_dt_wide) == nrow(b_0_beta_0_dt))
plot_params_dt_lng <- melt.data.table(data = plot_params_dt_wide,
                                      id.vars = c("t", "rep", "estimator"),
                                      measure.vars = c("b_0", "beta_0", "beta_1"), 
                                      variable.name = "parameter", variable.factor = FALSE,
                                      value.name = "value", value.factor = FALSE)
plot_params_dt_lng[,
                   parameter := factor(parameter, 
                                       levels = c("b_0", "beta_0", "beta_1"),
                                       labels = c("$\\alpha (t)$", "$\\beta_0 (t)$", "$\\beta_1 (t)$"))]
plot_params_dt_lng[, estimator := factor(estimator, levels = c("ols", "br"), labels = c("OLS", "Bias-Corrected"))]


plot_params_dt_lng_average <- plot_params_dt_lng[, 
                                                 .(value = mean(value)),
                                                 by = .(
                                                   t, estimator, parameter
                                                 )
]


truth_dataset <- data.table(
  t = rep(grid_points, times = 3),
  parameter = rep(c("$\\alpha (t)$", "$\\beta_0 (t)$", "$\\beta_1 (t)$"), each = length(grid_points)),
  value = c(rep(0, length(grid_points)), rep(-1, length(grid_points)), damping_fun_test(grid_points))
)

truth_dataset[,
              parameter := factor(parameter, 
                                  levels = c("$\\alpha (t)$", "$\\beta_0 (t)$", "$\\beta_1 (t)$"),
                                  labels = c("$\\alpha (t)$", "$\\beta_0 (t)$", "$\\beta_1 (t)$"))]

results_plot <- ggplot(data = plot_params_dt_lng) +
  aes(x = t, y = value) +
  facet_wrap(~ parameter, scales = "free_y", ncol = 3) +
  geom_line(alpha = 0.15, aes(group = interaction(rep, estimator), colour = estimator)) +
  geom_line(data = plot_params_dt_lng_average, aes(group = estimator, colour = estimator), lwd = 1.1, lty = 1) +
  geom_line(data = truth_dataset, col = 1, lty = 2) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.title.y = element_blank()) +
  labs(x = "$t$") +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 0.75))) +
  scale_color_manual(values = c("#CC79A7", "#009E73")) +
  scale_x_continuous(breaks = c(0, pi, 2 * pi, 3 * pi, 4 * pi), labels = c("$0$", "$\\pi$", "$2 \\pi$", "$3 \\pi$", "$4 \\pi$"),
                     limits = c(0, 4 * pi), expand = rep(0.025, 2))
results_plot


tikz(here::here("outputs", "SHM", "paper-plots", "general-dhm-tv-result.tex"),
     width = 1.5 * doc_width_inches, 
     height =  0.65 * doc_width_inches,
     standAlone = TRUE)
results_plot
dev.off()  

tinytex::lualatex(here::here("outputs", "SHM", "paper-plots", "general-dhm-tv-result.tex"))

