# ------------------------------------------------------------------------#
# Create Plot of the Decomposition of E[x(t)]
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
prepared_data <- readRDS(here::here("outputs", "real-data", "prepared-data-01.rds"))
x <- prepared_data$x
Dx <- prepared_data$Dx
D2x <- prepared_data$D2x
milliseconds <- prepared_data$milliseconds
beta_hat <- pda_results$pda_result$beta


# Using un-corrected parameter estimates: ---------------------------------
# (i.e., initial estimates)
alpha_fun_uc <- approxfun(x = milliseconds, beta_hat[,1,1])
beta_0_fun_uc <- approxfun(x = milliseconds, beta_hat[,2,1])
beta_1_fun_uc <- approxfun(x = milliseconds, beta_hat[,3,1])

dynamics_equations_hom_uc <- function(t, y, ...) {
  with(as.list(c(y)),{
    # rate of change
    dX <- Y
    dY <-  beta_0_fun_uc(t) * X   + beta_1_fun_uc(t) * Y
    list(c(dX, dY))
  }
  )
}


u1_mat_uc <- deSolve::lsoda(y = c(X = 1, Y = 0),
                            times = milliseconds,
                            func = dynamics_equations_hom_uc,
                            tcrit = max(milliseconds))
u2_mat_uc <- deSolve::lsoda(y = c(X = 0, Y = 1),
                            times = milliseconds,
                            func = dynamics_equations_hom_uc,
                            tcrit = max(milliseconds))

plot_dt_uc <- data.table(t=milliseconds,
                         "$u_1(t)$" = u1_mat_uc[,2],
                         "$u_2(t)$" = u2_mat_uc[,2])


# Using corrected Parameter Estimates: ------------------------------------
# (i.e., final estimates)
alpha_fun_c <- approxfun(x = milliseconds, beta_hat[,1,11])
beta_0_fun_c <- approxfun(x = milliseconds, beta_hat[,2,11])
beta_1_fun_c <- approxfun(x = milliseconds, beta_hat[,3,11])

dynamics_equations_hom_c <- function(t, y, ...) {
  with(as.list(c(y)),{
    # rate of change
    dX <- Y
    dY <-  beta_0_fun_c(t) * X   + beta_1_fun_c(t) * Y
    list(c(dX, dY))
  }
  )
}


u1_mat_c <- deSolve::lsoda(y = c(X = 1, Y = 0),
                           times = milliseconds,
                           func = dynamics_equations_hom_c,
                           tcrit = max(milliseconds))
u2_mat_c <- deSolve::lsoda(y = c(X = 0, Y = 1),
                           times = milliseconds,
                           func = dynamics_equations_hom_c,
                           tcrit = max(milliseconds))

plot_dt_c <- data.table(t=milliseconds,
                        "$u_1(t)$" = u1_mat_c[,2],
                        "$u_2(t)$" = u2_mat_c[,2])



# Combine uncorrected and corrected estimates: ----------------------------
plot_dt_uc$correction <- "Initial PDA Estimate"
plot_dt_c$correction <- "Bias-Corrected Estimate"
plot_dt <- rbind(plot_dt_c, plot_dt_uc)
plot_dt$correction <- factor(plot_dt$correction,
                             levels = c("Initial PDA Estimate",
                                        "Bias-Corrected Estimate"))
plot_dt_lng <- melt.data.table(data = plot_dt, id.vars = c("t", "correction"))

(p_both <- ggplot(data = plot_dt_lng) +
    aes(x = t, y = value, group = variable, colour = variable) +
    facet_wrap(~ correction) +
    geom_line() +
    labs(x = "$t$", y = "$u(t)$", title = "PDA Basis Functions") +
    theme(legend.position = "bottom", legend.title = element_blank()))


# Do Decomposition of $\\mathbb{E}[x(t)]$: --------------------------------------------

# compute first term:
mu_0 <- c(mean(x[1,]), mean(Dx[1,]))
first_term_c <- (matrix(mu_0, nrow = 1, ncol = 2) %*% t(cbind(u1_mat_c[,2], u2_mat_c[,2])))[1,]
first_term_uc <- (matrix(mu_0, nrow = 1, ncol = 2) %*% t(cbind(u1_mat_uc[,2], u2_mat_uc[,2])))[1,]

plot(milliseconds, first_term_uc, type = "l")
lines(milliseconds, first_term_c, col = 4)

# compute second term:
# (for both hom and inhom).
dynamics_equations_inhom_c <- function(t, y, ...) {
  with(as.list(c(y)),{
    # rate of change
    dX <- Y
    dY <-  beta_0_fun_c(t) * X   + beta_1_fun_c(t) * Y + alpha_fun_c(t)
    list(c(dX, dY))
  }
  )
}

dynamics_equations_inhom_uc <- function(t, y, ...) {
  with(as.list(c(y)),{
    # rate of change
    dX <- Y
    dY <-  beta_0_fun_uc(t) * X   + beta_1_fun_uc(t) * Y + alpha_fun_uc(t)
    list(c(dX, dY))
  }
  )
}


second_term_uc <- (deSolve::lsoda(y = c(X = 1, Y = 0), times = milliseconds, 
                                 func = dynamics_equations_inhom_uc, tcrit = max(milliseconds))[,2] -
  deSolve::lsoda(y = c(X = 1, Y = 0), times = milliseconds, 
                 func = dynamics_equations_hom_uc, tcrit = max(milliseconds))[,2])

second_term_c <- (deSolve::lsoda(y = c(X = 1, Y = 0), times = milliseconds, 
                                  func = dynamics_equations_inhom_c, tcrit = max(milliseconds))[,2] -
                     deSolve::lsoda(y = c(X = 1, Y = 0), times = milliseconds, 
                                    func = dynamics_equations_hom_c, tcrit = max(milliseconds))[,2])


# Plot of Mean Decomposition: ---------------------------------------------

mean_decomp_dt <- data.table(
  t = rep(milliseconds, times = 2),
  correction = rep(x = c("Initial", "Final"), each = length(milliseconds)),
  first_term = c(first_term_uc, first_term_c),
  second_term = c(second_term_uc, second_term_c)
)

mean_decomp_dt_lng <- melt.data.table(data = mean_decomp_dt, 
                                      id.vars = c("t", "correction"), 
                                      variable.name = "term", 
                                      variable.factor = FALSE,
                                      value.factor = FALSE)

mean_decomp_dt_lng[,
  term := factor(term, levels = c("$\\mathbb{E}[x(t)]$", "first_term", "second_term"),
                 labels = c("$\\mathbb{E}[x(t)]$",
                            "$\\Phi(t, 0) \\boldsymbol{\\mu}_0$",
                            "$\\int_0^t\\boldsymbol{\\Phi}(t, s) (0, \\alpha(s))^\\top \\mathrm{d}s$"))
]

mean_dt <- data.table(t = milliseconds, term = factor("$\\mathbb{E}[x(t)]$", 
                                                      levels = c("$\\mathbb{E}[x(t)]$",
                                                                 "$\\Phi(t, 0) \\boldsymbol{\\mu}_0$",
                                                                 "$\\int_0^t\\boldsymbol{\\Phi}(t, s) (0, \\alpha(s))^\\top \\mathrm{d}s$")),
                                                      value = apply(x,1,mean))



(p <- ggplot(data = mean_decomp_dt_lng) +
  aes(x = t, y = value, group = term, colour = term) +
  facet_wrap(~ correction) +
  geom_line(data = mean_dt) +
  geom_line() +
  ylim(c(-0.2, 0.2)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.subtitle = element_text(hjust = 0.5)) +
  scale_colour_manual(values = c("black", "cornflowerblue", "red3"),
                      breaks = c("$\\mathbb{E}[x(t)]$",
                                 "$\\Phi(t, 0) \\boldsymbol{\\mu}_0$",
                                 "$\\int_0^t\\boldsymbol{\\Phi}(t, s) (0, \\alpha(s))^\\top \\mathrm{d}s$")) +
  labs(x = "$t$", y = "Value", 
       title = "Decomposition of",
       subtitle = "$\\mathbb{E}[x(t)] =\\boldsymbol{\\Phi}(t, 0) \\boldsymbol{\\mu}_0 + \\int_0^t\\boldsymbol{\\Phi}(t, s) (0, \\alpha(s))^\\top \\mathrm{d}s$") +
  guides(colour = guide_legend(override.aes = list(linewidth = 1))))

tikz(here::here("outputs", "real-data", "paper-plots", "com-mean-decomposition.tex"),
     width = doc_width_inches,
     height = 0.6 *doc_width_inches,
     standAlone = TRUE)
print(p)  
dev.off()

tinytex::lualatex(here::here("outputs", "real-data", "paper-plots", "com-mean-decomposition.tex"))

