# Packages: ---------------------------------------------------------------
library(data.table) # CRAN v1.14.2
library(tikzDevice) # CRAN v0.12.3.1
library(ggplot2)    # CRAN v3.4.0


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


# Import results: ---------------------------------------------------------
results_list <- readRDS(file = here::here("outputs", "VDP", "LNA-demo-02.rds"))
settings <- results_list$settings
grid_points <- results_list$grid_points



# Reshape data: -----------------------------------------------------------
x_true_dt <- data.table(t(results_list$results_x_true))
names(x_true_dt) <- paste(grid_points)
x_true_dt <- data.table(settings, x_true_dt)
x_true_dt_lng <- melt.data.table(
  data = x_true_dt, 
  id.vars = c("rep", "mu"),
  measure.vars =  paste(grid_points),
  variable.name = "t", 
  value.name = "fun_val", 
  variable.factor = FALSE,
  value.factor = FALSE,
  verbose = TRUE
)
x_true_dt_lng$var <- "x"
x_true_dt_lng$method <- "true"

y_true_dt <- data.table(t(results_list$results_y_true))
names(y_true_dt) <- paste(grid_points)
y_true_dt <- data.table(settings, y_true_dt)
y_true_dt_lng <- melt.data.table(
  data = y_true_dt, 
  id.vars = c("rep", "mu"),
  measure.vars =  paste(grid_points),
  variable.name = "t", 
  value.name = "fun_val", 
  variable.factor = FALSE,
  value.factor = FALSE,
  verbose = TRUE
)
y_true_dt_lng$var <- "y"
y_true_dt_lng$method <- "true"



x_lin_dt <- data.table(t(results_list$results_x_lin))
names(x_lin_dt) <- paste(grid_points)
x_lin_dt <- data.table(settings, x_lin_dt)
x_lin_dt_lng <- melt.data.table(
  data = x_lin_dt, 
  id.vars = c("rep", "mu"),
  measure.vars =  paste(grid_points),
  variable.name = "t", 
  value.name = "fun_val", 
  variable.factor = FALSE,
  value.factor = FALSE,
  verbose = TRUE
)
x_lin_dt_lng$var <- "x"
x_lin_dt_lng$method <- "lin"

y_lin_dt <- data.table(t(results_list$results_y_lin))
names(y_lin_dt) <- paste(grid_points)
y_lin_dt <- data.table(settings, y_lin_dt)
y_lin_dt_lng <- melt.data.table(
  data = y_lin_dt, 
  id.vars = c("rep", "mu"),
  measure.vars =  paste(grid_points),
  variable.name = "t", 
  value.name = "fun_val", 
  variable.factor = FALSE,
  value.factor = FALSE,
  verbose = TRUE
)
y_lin_dt_lng$var <- "y"
y_lin_dt_lng$method <- "lin"



# Get data ready for plotting: --------------------------------------------
plot_dt <- rbind(
  x_lin_dt_lng,
  x_true_dt_lng, 
  y_lin_dt_lng,
  y_true_dt_lng
)

plot_dt[, t := as.numeric(t)]
plot_dt[, rep := factor(rep)]
plot_dt[, method := factor(method, levels= c("true", "lin"), labels = c("True Non-Linear Model",
                                                                        "Linearised Approximation"))]
plot_dt[, mu := factor(mu, 
                       levels = c(0.5, 1, 2),
                       labels = paste0("$\\mu = ", c(0.5, 1, 2), "$"))]

plot_dt[, var := factor(var, levels = c("x", "y"), labels = c("$x(t)$", "$y(t)$"))]


# Plot --------------------------------------------------------------------

p <- ggplot(data = plot_dt) +
  aes(x = t, y = fun_val, group = interaction(rep, method), colour = rep, lty = method) +
  facet_wrap(var ~ mu, nrow = 2, ncol = 3, scales = "free_y") +
  geom_line() +
  labs(x = "$t$", y = "Function Value") +
  scale_color_discrete(guide = "none") +
  scale_linetype_discrete("Generative Model: ") 



# Save: -------------------------------------------------------------------

tikz(here::here("outputs", "VDP", "paper-plots", "LNA-demo-02.tex"),
     width = (3/2) * doc_width_inches, 
     height =  doc_width_inches)
p
dev.off()


