# Packages: ---------------------------------------------------------------
library(data.table) # CRAN v1.14.2
library(ggplot2)    # CRAN v3.4.0
library(tikzDevice)

# Scripts: ----------------------------------------------------------------
source(here::here("code", "VDP-helper-functions.R"))


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


# Import and unpack results: ----------------------------------------------
results_list <- readRDS(file = here::here(
  "outputs",
  "VDP",
  "vdp-simulation-partial_derivatives_note.rds"))
grid_points <- seq(0, 13, length.out = 200) 
partial_derivatives_array <- results_list$partial_derivatives_array
beta_array <- results_list$beta_array

par(mfrow = c(2, 2))
matplot(beta_array[,2,1,], type = "l", col = 2)
matlines(partial_derivatives_array[,1,1,], col = 3)

matplot(beta_array[,3,1,], type = "l", col = 2)
matlines(partial_derivatives_array[,2,1,], col = 3)

matplot(beta_array[,2,2,], type = "l", col = 2)
matlines(partial_derivatives_array[,1,2,], col = 3)

matplot(beta_array[,3,2,], type = "l", col = 2)
matlines(partial_derivatives_array[,2,2,], col = 3)

ols_dt <- data.table(
  t_grid = rep(grid_points, 4),
  parameter = rep(c( "beta_xx", "beta_xy", "beta_yx", "beta_yy"), each = length(grid_points)),
  rbind(beta_array[,2,1,],
        beta_array[,3,1,],
        beta_array[,2,2,],
        beta_array[,3,2,]))
names(ols_dt)[-c(1:2)] <- paste0("rep_", 1:50)


jacobian_dt <- data.table(
  t_grid = rep(grid_points, 4),
  parameter = rep(c( "beta_xx", "beta_xy", "beta_yx", "beta_yy"), each = length(grid_points)),
  rbind(partial_derivatives_array[,1,1,],
        partial_derivatives_array[,2,1,],
        partial_derivatives_array[,1,2,],
        partial_derivatives_array[,2,2,]))
names(jacobian_dt)[-c(1:2)] <- paste0("rep_", 1:50)

ols_dt_lng <- melt.data.table(ols_dt, id.vars = c("t_grid", "parameter"),
                              measure.vars = paste0("rep_", 1:50),
                              variable.name = "rep",
                              value.name = "param_est", 
                              variable.factor = FALSE,
                              value.factor = FALSE)

jacobian_dt_lng <- melt.data.table(jacobian_dt, id.vars = c("t_grid", "parameter"),
                              measure.vars = paste0("rep_", 1:50),
                              variable.name = "rep",
                              value.name = "param_est", 
                              variable.factor = FALSE,
                              value.factor = FALSE)

ols_dt_lng$method  <- "ols"
jacobian_dt_lng$method <- "np"

plot_dt_lng <- rbind(ols_dt_lng, jacobian_dt_lng)



# Generate true values: ---------------------------------------------------
mu_ic <- c(1.9922213675594, -0.910974076470711)
mu<-1
limit_cycle_data <- generate_vdp_deterministic(mu = 1, 
                                               grid_points = grid_points,
                                               y0 = mu_ic[2],
                                               x0 = mu_ic[1])
limit_x <- limit_cycle_data[,2]
limit_y <- limit_cycle_data[,3]

b_0x_limit <- mu * (limit_x - (1/3)*limit_x^3 - limit_y)  - mu * (1 - limit_x^2) * limit_x + mu * limit_y
beta_xx_limit <- mu * (1-limit_x^2)

truth_dt <- data.table(
  t_grid = grid_points,
  b0_x = b_0x_limit,
  b0_y = 0,
  beta_xx = beta_xx_limit,
  beta_yx = 1/mu,
  beta_xy = - mu,
  beta_yy = 0)
truth_dt_lng <- melt.data.table(truth_dt, 
                                id.vars = "t_grid",
                                variable.name = "parameter",
                                variable.factor = FALSE, 
                                value.factor = FALSE,
                                value.name = "param_est")
truth_dt_lng <- truth_dt_lng[parameter %in% c( "beta_xx", "beta_xy", "beta_yx", "beta_yy")]


# -------------------------------------------------------------------------

truth_dt_lng[, parameter_label := factor(parameter,
                                         levels = c( "beta_xx", "beta_xy", "beta_yx", "beta_yy"),
                                         labels = c( "$\\beta_{xx}(t)$", "$\\beta_{xy}(t)$", "$\\beta_{yx}(t)$", "$\\beta_{yy}(t)$"))]
plot_dt_lng[, parameter_label := factor(parameter,
                                         levels = c( "beta_xx", "beta_xy", "beta_yx", "beta_yy"),
                                         labels = c( "$\\beta_{xx}(t)$", "$\\beta_{xy}(t)$", "$\\beta_{yx}(t)$", "$\\beta_{yy}(t)$"))]

plot_dt_lng[, method_label := factor(method,
                                     levels = c("ols", "np"),
                                     labels = c("OLS", "Non-Parametric"))]
truth_dt_lng$method_label = "True"


p <- ggplot(data = plot_dt_lng) +
  aes(x = t_grid, y = param_est, colour = method_label) +
  facet_wrap(~ parameter_label, scales = "free_y") +
  geom_line(aes(group = interaction(rep, method)), alpha = 0.25) +
  geom_line(data = truth_dt_lng, linewidth = 0.75) +
  xlab("$t$") +
  ylab("$\\beta(t)$") +
  labs(colour = "Method:") +
  theme(legend.position = "bottom") +
  scale_colour_manual(values = c("red3", "cornflowerblue", "black")) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1, alpha = 1)))

p
tikz(here::here("outputs", "VDP", "paper-plots", "partial-derivative-plot.tex"),
     width = 1 * doc_width_inches, 
     height = 0.9 *  doc_width_inches,
     standAlone = TRUE)
p
dev.off()

tinytex::lualatex(file = here::here("outputs", "VDP", "paper-plots", "partial-derivative-plot.tex"))

