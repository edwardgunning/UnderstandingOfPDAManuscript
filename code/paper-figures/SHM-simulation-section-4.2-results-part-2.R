source(here::here("code", "SHM-bias-calculation-function.R"))
source(here::here("code", "SHM-helper-functions.R"))
library(data.table) # Extension of `data.frame`
library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics
library(tikzDevice) # R Graphics Output in LaTeX Format
options(tikzLatexPackages = c(getOption( "tikzLatexPackages" ),"\\usepackage{amsfonts}"))
options(tikzLatexPackages = c(getOption( "tikzLatexPackages" ),"\\usepackage{amsmath}"))
options(tikzLatexPackages = c(getOption( "tikzLatexPackages" ),"\\usepackage{amssymb}"))
#                                                     "simulation-section-4.2.rds"))
# results_array <- simulation_result_list$results_array
simulation_result_list <- readRDS(here::here("outputs",
                                             "SHM",
                                             "simulation-results",
                                             "simulation-section-4.2-updated.rds"))
grid_range <- c(0, 2 * pi)
n_grid <- 100
grid_points <- seq(from = grid_range[1],
                   to = grid_range[2],
                   length.out = n_grid)
n_sim <- 500
results_array <- array(unlist(simulation_result_list$results_list), dim = c(n_grid, 4, n_sim))



bias_result <- apply(results_array + 1, c(1, 2), mean)
bias_SE <- apply(results_array + 1, c(1, 2), function(x) {
  sd(x) / sqrt(500)
})

bias_confint_lower <- bias_result - 2 * bias_SE
bias_confint_upper <- bias_result + 2 * bias_SE

bias_result_dt <- data.table(grid_points, bias_result, "estimate")
bias_lower_dt <- data.table(grid_points, bias_confint_lower, "lower")
bias_upper_dt <- data.table(grid_points, bias_confint_upper, "upper")

names(bias_result_dt) <- names(bias_lower_dt) <- names(bias_upper_dt) <- c("t_grid",
                                                                           "No correction",
                                                                           "1 iteration",
                                                                           paste0(2:3, " iterations"),
                                                                           "quantity")
bias_plot_dt <- rbind(bias_result_dt, bias_lower_dt, bias_upper_dt)
bias_plot_dt_lng <- melt.data.table(data = bias_plot_dt,
                                    id.vars = c("t_grid", "quantity"), 
                                    measure.vars = c("No correction",
                                                     "1 iteration",
                                                     paste0(2:3, " iterations")), 
                                    variable.name = "method",
                                    variable.factor = TRUE,
                                    value.name = "bias_t", value.factor = FALSE,
                                    verbose = TRUE)
bias_plot_dt_lng
p <- ggplot(bias_plot_dt_lng) +
  aes(x = t_grid, y = bias_t,
      group = interaction(method, quantity),
      linetype = quantity,
      color = method) +
  geom_line() +
  scale_linetype_manual(values = c(1, 2, 2), guide = NULL) +
  theme_bw() +
  theme(legend.position = c(0.75, 0.831),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(size = 8.75)) +
  scale_x_continuous(breaks = c(0, pi, 2 * pi), labels = c("$0$", "$\\pi$", "$2 \\pi$"),
                     limits = c(0, 2 * pi), expand = rep(0.025, 2)) +
  labs(x = "$t$",
       y = "$\\mathbb{E}[\\widehat{\\beta}_0 (t)] - \\beta_0 (t)$",
       title = "Bias") +  #: $E[\\widehat{\\beta}_0 (t)] - \\beta_0 (t)$") +
  # geom_hline(yintercept = 0, linetype = "dotted", col = "black", size = 0.75) +
  theme(axis.title = element_text(size = 9),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        plot.title = element_text(size = 9, hjust = 0.5, face = "bold")) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.1)))
p
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937
tikz(here::here("outputs", "SHM", "paper-plots", "bias-corrected-beta.tex"), 
     width = (0.51 * doc_width_inches),
     height = (0.46 * doc_width_inches))
print(p)
dev.off()

