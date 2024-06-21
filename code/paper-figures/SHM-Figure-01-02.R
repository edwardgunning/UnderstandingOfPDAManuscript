

# Packages ----------------------------------------------------------------
library(tikzDevice)  # CRAN v0.12.3.1 
library(data.table)  # CRAN v1.14.2
library(ggplot2)     # CRAN v3.3.5
library(randomcoloR) # CRAN v1.1.0.1


# Functions to generate data ----------------------------------------------
source(here::here("code", "SHM-helper-functions.R"))


# Preliminaries -----------------------------------------------------------
N <- 20 # Number of replicate curves
grid_range <- c(0, 2 * pi)
n_grid <- 100
grid_points <- seq(from = grid_range[1],
                   to = grid_range[2],
                   length.out = n_grid)

# Simulate functional observations ----------------------------------------
# Random Palette was gicing different colours every time, even with set seed
# https://github.com/ronammar/randomcoloR/pull/14
# set.seed(1996)
# palette <- randomColor(N) # simulate random colours to plot each curve
palette <- c("#2c6f89", "#82c62f", "#b9f293", "#8d78db", "#f21fc1", "#59ed1e", 
             "#6fe2db", "#f7c3bb", "#6996db", "#aefcd9", "#cc3079", "#d86c97", 
             "#7b9b06", "#266993", "#c4a235", "#e8f77b", "#562aa3", "#d85db9", 
             "#9fd4e0", "#8eb8ed")
set.seed(1996)
dataset_fig <- generate_SHM_dataset(n = N, 
                                    grid_points = grid_points,
                                    sigma = 0.25,
                                    mu_init = c(0, 0), 
                                    intensity = 2,
                                    sigma_init = diag(rep(0.05, 2)))

saveRDS(object = list(datset_fig = dataset_fig, seed = 1996, time = timestamp()),
        file = here::here("outputs", "SHM", "simulation-results", "dataset-fig.rds"))

# Reshape data for plotting -----------------------------------------------
xfig <- dataset_fig$x
D2xfig <- calculate_deriv(grid_points = grid_points, x = xfig, norder = 2)
Dxfig <- calculate_deriv(grid_points = grid_points, x = xfig, norder = 1)
colnames(xfig) <- colnames(D2xfig) <- colnames(Dxfig) <- paste0("rep", seq_len(N))

xfig_dt <- data.table(t_grid = grid_points, xfig)
D2xfig_dt <- data.table(t_grid = grid_points, D2xfig)
Dxfig_dt <- data.table(t_grid = grid_points, Dxfig)

xfig_dt_lng <- melt.data.table(data = xfig_dt, id.vars = "t_grid",
                               variable.name = "rep_fun",
                               variable.factor = TRUE,
                               value.name = "xt", 
                               value.factor = FALSE,
                               verbose = TRUE)
D2xfig_dt_lng <- melt.data.table(data = D2xfig_dt, id.vars = "t_grid",
                               variable.name = "rep_fun",
                               variable.factor = TRUE,
                               value.name = "D2xt", 
                               value.factor = FALSE,
                               verbose = TRUE)

Dxfig_dt_lng <- melt.data.table(data = Dxfig_dt, id.vars = "t_grid",
                                 variable.name = "rep_fun",
                                 variable.factor = TRUE,
                                 value.name = "Dxt", 
                                 value.factor = FALSE,
                                 verbose = TRUE)

plot_dt <- merge(merge(xfig_dt_lng,
                 D2xfig_dt_lng, by = c("t_grid", "rep_fun")),
                 Dxfig_dt_lng, by = c("t_grid", "rep_fun"))


# Plot and save -----------------------------------------------------------
y_lims_square <- c(-max(xfig, D2xfig, Dxfig), max(xfig, D2xfig, Dxfig))
p1 <- ggplot(data = plot_dt) +
  aes(x = xt, y = D2xt, group = rep_fun, colour = rep_fun) +
  geom_path() +
  theme_bw() +
  scale_color_manual(values = palette) +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = c(-1, 0, 1), limits = y_lims_square) +
  scale_y_continuous(breaks = c(-1, 0, 1), limits = y_lims_square) +
  labs(x = "$x(t)$", y = "$D^2 x(t)$", title = "$D^2 x(t) \\sim x(t)$") +
  geom_abline(slope = - 1, linetype = "dotted", col = "black", size = 1) +
  theme(axis.title = element_text(size = 8),
        plot.title = element_text(size = 9, hjust = 0.5))

p2 <- ggplot(data = plot_dt) +
  aes(x = xt, y = Dxt, group = rep_fun, colour = rep_fun) +
  geom_path() +
  theme_bw() +
  scale_color_manual(values = palette) +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = c(-1, 0, 1), limits = y_lims_square) +
  scale_y_continuous(breaks = c(-1, 0, 1), limits = y_lims_square) +
  labs(x = "$x(t)$", y = "$D x(t)$", title = "$Dx(t) \\sim x(t)$") +
  theme(axis.title = element_text(size = 8),
        plot.title = element_text(size = 9, hjust = 0.5))


p3 <- ggplot(data = plot_dt) +
  aes(x = Dxt, y = D2xt, group = rep_fun, colour = rep_fun) +
  geom_path() +
  theme_bw() +
  scale_color_manual(values = palette) +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = c(-1, 0, 1), limits = y_lims_square) +
  scale_y_continuous(breaks = c(-1, 0, 1), limits = y_lims_square) +
  labs(x = "$Dx(t)$", y = "$D^2 x(t)$", title =  "$D^2x(t) \\sim Dx(t)$") +
  theme(axis.title = element_text(size = 8),
        plot.title = element_text(size = 9, hjust = 0.5))

full_plot <- ggpubr::ggarrange(plotlist = list(p2, p1, p3),
                               nrow = 1,
                               ncol = 3, 
                               labels = c("(a)", "(b)", "(c)"), 
                               align = "hv", font.label = list(size = 8))
full_plot
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937
tikz(here::here("outputs", "SHM", "paper-plots", "shm-phase-plane.tex"),
     width = doc_width_inches, 
     height = doc_width_inches/3)
print(full_plot)
dev.off()





# Second plot -------------------------------------------------------------

epsfig_dt <- data.table(t_grid = grid_points, dataset_fig$noise)

epsfig_dt_lng <- melt.data.table(data = epsfig_dt, id.vars = "t_grid",
                               variable.name = "rep_fun",
                               variable.factor = TRUE,
                               value.name = "epst", 
                               value.factor = FALSE,
                               verbose = TRUE)



p4 <- ggplot(data = plot_dt) +
  aes(x = t_grid, y = D2xt, group = rep_fun, colour = rep_fun) +
  geom_path() +
  theme_bw() +
  scale_color_manual(values = palette) +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = c(0, pi, 2 * pi), labels = c("$0$", "$\\pi$", "$2 \\pi$"),
                     limits = c(0, 2 * pi), expand = rep(0.025, 2)) +
  scale_y_continuous(breaks = c(-1, 0, 1), limits = y_lims_square) +
  labs(x = "$t$", y = "$x(t)$", title =  "Generated Functions") +
  theme(axis.title = element_text(size = 9),
        plot.title = element_text(size = 9, hjust = 0.5, face = "bold"))
p5 <- ggplot(data = epsfig_dt_lng) +
  aes(x = t_grid, y = epst, group = rep_fun, colour = rep_fun) +
  geom_path() +
  theme_bw() +
  scale_color_manual(values = palette) +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = c(0, pi, 2 * pi), labels = c("$0$", "$\\pi$", "$2 \\pi$"),
                     limits = c(0, 2 * pi), expand = rep(0.025, 2)) +
  scale_y_continuous(breaks = c(-1, 0, 1), limits = y_lims_square) +
  labs(x = "$t$", y = "$\\epsilon(t)$", title =  "Stochastic Disturbances") +
  theme(axis.title = element_text(size = 9),
        plot.title = element_text(size = 9, hjust = 0.5, face = "bold"))


tikz(here::here("outputs", "SHM", "paper-plots", "shm-dataset.tex"), 
     width = (0.75 * doc_width_inches),
     height = (0.75 * doc_width_inches)/2)
full_plot_2 <- ggpubr::ggarrange(plotlist = list(p4, p5),
                               nrow = 1,
                               ncol = 2, 
                               labels = c("(a)", "(b)"), 
                               align = "hv",
                               font.label = list(size = 8))
full_plot_2
dev.off()
