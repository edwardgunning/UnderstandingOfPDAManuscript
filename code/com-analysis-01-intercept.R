

# Load packages: ----------------------------------------------------------
library(data.table) # CRAN v1.14.2
library(ggplot2)    # CRAN v3.4.2
library(tikzDevice) # CRAN v0.12.3.1

# Functions: --------------------------------------------------------------
source(here::here("code", "theme_gunning.R"))
source(here::here("code", "calculate_non_parametric_fit.R"))

# Some graphics settings: -------------------------------------------------
theme_gunning()
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937
options(tikzLatexPackages = c(getOption( "tikzLatexPackages" ),"\\usepackage{amsfonts}"))
options(tikzLatexPackages = c(getOption( "tikzLatexPackages" ),"\\usepackage{amsmath}"))
options(tikzLatexPackages = c(getOption( "tikzLatexPackages" ),"\\usepackage{amssymb}"))


# Import and Unpack Saved Results: ----------------------------------------
pda_results <- readRDS(file = here::here("outputs", "real-data", "pda-result-01.rds"))
prepared_data <- readRDS(here::here("outputs", "real-data", "prepared-data-01.rds"))
x <- prepared_data$x
Dx <- prepared_data$Dx
D2x <- prepared_data$D2x
milliseconds <- prepared_data$milliseconds

# Do a simple PDA with NO intercept: --------------------------------------
beta_mat <- matrix(NA, nrow = length(milliseconds), ncol = 2)
beta_hat <- pda_results$pda_result$beta
for(tind in seq_along(milliseconds)) {
  beta_mat[tind,] <- coefficients(lm(D2x[tind,]~0 + x[tind, ] + Dx[tind,]))
}



# Do a non-parametric fit and get gradients -------------------------------
non_parametric_fit <- calculate_non_parametric_fit(x = x,
                                                   Dx = Dx,
                                                   D2x = D2x,
                                                   grid_points = milliseconds)


# Quick and dirty plots for comparison: -----------------------------------
plot(beta_mat[,1], type = "l", ylim = range(beta_mat[,1], beta_hat[,2,1]))
lines(beta_hat[,2,1], col=2)
lines(non_parametric_fit$B0_approx, col=4)

plot(beta_mat[,2], type = "l", ylim = range(beta_mat[,2], beta_hat[,3,1]))
lines(beta_hat[,3,1], col=2)
lines(non_parametric_fit$B1_approx, col = 4)




# Create publishable-squality plots w beamer font: ------------------------
intercept_dt <- data.table(t = milliseconds, 
                           beta_hat[,2:3,1])
no_intercept_dt <- data.table(t = milliseconds, 
                              beta_mat)

names(intercept_dt)[2:3] <- names(no_intercept_dt)[2:3] <- c("$\\beta_0(t)$", "$\\beta_1(t)$")

non_parametric_dt <- data.table(
  t = milliseconds,
  "$\\beta_0(t)$" = non_parametric_fit$B0_approx,
  "$\\beta_1(t)$" = non_parametric_fit$B1_approx,
  model = "Gradient of Non-linear Fit"
)

intercept_dt$model <- "PDA Intercept Included"
no_intercept_dt$model <- "PDA No Intercept Included"


plot_dt <- rbind(no_intercept_dt, intercept_dt, non_parametric_dt)
plot_dt_lng <- melt.data.table(data = plot_dt, id.vars = c("t", "model"))
p <- ggplot(plot_dt_lng) +
  aes(x = t, y = value, group = model, colour = model) +
  geom_hline(yintercept = 0, col = "grey") +
  geom_line() +
  facet_wrap(~ variable) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  labs(x = "$t$", y = "$\\beta(t)$") +
  guides(colour = guide_legend(override.aes = list(linewidth = 1)))


tikz(here::here("outputs", "real-data", "paper-plots", "pda-intercept-estimates.tex"), 
     width = doc_width_inches,
     height = 0.5 * doc_width_inches,
     standAlone = TRUE)
p
dev.off()
tinytex::lualatex(here::here("outputs", "real-data", "paper-plots", "pda-intercept-estimates.tex"))


options(tikzLatexPackages = c(getOption( "tikzLatexPackages" ),"\\renewcommand{\\familydefault}{\\sfdefault}"))

tikz(here::here("outputs", "real-data", "rough", "pda-intercept-estimates.tex"), 
     width = doc_width_inches,
     height = 0.5 * doc_width_inches,
     standAlone = TRUE)
p
dev.off()
tinytex::lualatex(here::here("outputs", "real-data", "rough", "pda-intercept-estimates.tex"))

