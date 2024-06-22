# ------------------------------------------------------------------------#
# Plot results of short simulation
# ------------------------------------------------------------------------#

# Packages: ---------------------------------------------------------------
library(data.table) # CRAN v1.14.2
library(ggplot2)    # CRAN v3.4.2
library(tikzDevice) # CRAN v0.12.3.1

# Some Graphics Settings: -------------------------------------------------
source(here::here("code", "theme_gunning.R"))
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937
theme_gunning()
theme_update(axis.text = element_text(size = 9))
options(tikzLatexPackages = c(getOption( "tikzLatexPackages" ),"\\usepackage{amsfonts}"))
options(tikzLatexPackages = c(getOption( "tikzLatexPackages" ),"\\usepackage{amsmath}"))
options(tikzLatexPackages = c(getOption( "tikzLatexPackages" ),"\\usepackage{amssymb}"))

# Import Results and Extract: ---------------------------------------------
pda_results <- readRDS(file = here::here("outputs", "real-data", "pda-result-01.rds"))
pda_simulation_results <- readRDS(file = here::here("outputs", "real-data", "pda-simulation-result-01.rds"))
prepared_data <- readRDS(here::here("outputs", "real-data", "prepared-data-01.rds"))
N <- prepared_data$N
x <- prepared_data$x
Dx <- prepared_data$Dx
D2x <- prepared_data$D2x
milliseconds <- prepared_data$milliseconds
beta_hat <- pda_results$pda_result$beta




initial_estimates <- sapply(pda_simulation_results$pda_simulation_results_list, function(x) {
  x$beta[,,1]
})

dim(initial_estimates)


initial_dt <- data.table(
  t = rep(milliseconds, times = 3),
  parameter = rep(c("$\\alpha(t)$", "$\\beta_0(t)$", "$\\beta_1(t)$"), each = length(milliseconds)),
  initial_estimates
)

matplot(initial_estimates[1:142,], type = "l")

names(initial_dt)[-c(1:2)] <- paste0("sim_rep_", seq_len(15))

initial_dt_lng <- melt.data.table(initial_dt, 
                                  id.vars = c("t", "parameter"),
                                  variable.name = "sim_rep", 
                                  variable.factor = FALSE, 
                                  value.name = "init_est")

ggplot(initial_dt_lng) +
  aes(x = t, y = init_est, group = sim_rep, colour = sim_rep) +
  facet_wrap(~ parameter) +
  geom_line(alpha = 0.5) +
  theme(legend.position = "none") +
  labs(x = "$t$")


final_estimates <- sapply(pda_simulation_results$pda_simulation_results_list, function(x) {
  x$beta[,,11]
})
final_dt <- data.table(
  t = rep(milliseconds, times = 3),
  parameter = rep(c("$\\alpha(t)$", "$\\beta_0(t)$", "$\\beta_1(t)$"), each = length(milliseconds)),
  final_estimates
)

matplot(final_estimates[1:142,], type = "l")

names(final_dt)[-c(1:2)] <- paste0("sim_rep_", seq_len(15))

final_dt_lng <- melt.data.table(final_dt, 
                                  id.vars = c("t", "parameter"),
                                  variable.name = "sim_rep", 
                                  variable.factor = FALSE, 
                                  value.name = "final_est")

ggplot(final_dt_lng) +
  aes(x = t, y = final_est, group = sim_rep, colour = sim_rep) +
  facet_wrap(~ parameter, scales = "free_y") +
  geom_line(alpha = 0.5) +
  theme(legend.position = "none") +
  labs(x = "$t$")



# -------------------------------------------------------------------------

plot_dt <- merge.data.table(initial_dt_lng, final_dt_lng,
                            all = TRUE, by = c("t", "parameter", "sim_rep"))

plot_dt_lng <- melt.data.table(plot_dt, id.vars = c("t", "parameter", "sim_rep"))


truth_dt <- data.table(
  t = rep(milliseconds, times = 6),
  parameter = rep(rep(c("$\\alpha(t)$", "$\\beta_0(t)$", "$\\beta_1(t)$"), each = length(milliseconds)), times = 2),
  variable = rep(c("init_est", "final_est"), each = 3 * length(milliseconds)),
  value = c(c(pda_results$pda_result$beta[,,1]), c(pda_results$pda_result$beta[,,11]))
)

truth_dt[, variable := factor(variable, levels = c("init_est", "final_est"),
                              labels = c("Initial Estimates", "Final Estimates"))]
plot_dt_lng[, variable := factor(variable, levels = c("init_est", "final_est"),
                              labels = c("Initial Estimates", "Final Estimates"))]

p <- ggplot(plot_dt_lng) +
  facet_grid(parameter ~ variable, scales = "free_y") +
  aes(x = t, y = value) +
  geom_line(alpha = 0.5, aes(group = sim_rep, colour = sim_rep)) +
  geom_line(data = truth_dt, linewidth = 0.75) +
  theme(legend.position = "none") +
  labs(x = "$t$", y = "Parameter Estimates")

tikz(here::here("outputs", "real-data", "paper-plots", "com-pda-simulation-results.tex"),
     width = doc_width_inches,
     height = doc_width_inches,
     standAlone = TRUE)
p
dev.off()

tinytex::lualatex(here::here("outputs", "real-data", "paper-plots", "com-pda-simulation-results.tex"))
