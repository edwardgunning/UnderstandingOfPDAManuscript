# ------------------------------------------------------------------------#
# Create Plot of the Decomposition of Cov[x(t)]
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
N <- prepared_data$N
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



# De-mean x ---------------------------------------------------------------
x_cent <- sweep(x, 1, apply(x, 1, mean), "-")
matplot(x_cent, type = "l")


# Calculate zero-input covariance: ----------------------------------------
Sigma_0 <- cov(cbind(x[1,], Dx[1,]))
U_mat_uc <- cbind(u1_mat_uc[,2], u2_mat_uc[,2])
U_mat_c <- cbind(u1_mat_c[,2], u2_mat_c[,2])

zero_input_covariance_uc <- U_mat_uc %*% Sigma_0 %*% t(U_mat_uc)
zero_input_covariance_c <- U_mat_c %*% Sigma_0 %*% t(U_mat_c)


# Use that total Cov = Zero-input + Zero-State to get Zero-State ----------
Cov_x <- cov(t(x))
zero_state_covariance_c <- Cov_x - zero_input_covariance_c
zero_state_covariance_uc <- Cov_x - zero_input_covariance_uc

stochastic_basis_obj_c <- eigen(zero_state_covariance_c)
stochastic_basis_obj_c$values <- stochastic_basis_obj_c$values[stochastic_basis_obj_c$values>1e-12]
stochastic_basis_obj_c$vectors <- stochastic_basis_obj_c$vectors[,stochastic_basis_obj_c$values>1e-12]
stochastic_basis_funs_c <- stochastic_basis_obj_c$vectors[, 1:4]

stochastic_basis_obj_uc <- eigen(zero_state_covariance_uc)
stochastic_basis_obj_uc$values <- stochastic_basis_obj_uc$values[stochastic_basis_obj_uc$values>1e-12]
stochastic_basis_obj_uc$vectors <- stochastic_basis_obj_uc$vectors[,stochastic_basis_obj_uc$values>1e-12]
stochastic_basis_funs_uc <- stochastic_basis_obj_uc$vectors[, 1:4]


# quick and dirty plots:
par(mfrow = c(1, 2))
matplot(milliseconds, stochastic_basis_funs_c, type = "l",lty = 1)
matplot(milliseconds, stochastic_basis_funs_uc, type = "l",lty = 1)




# Now do OLS regression and Decompose Variability: ------------------------
combined_basis_c <- cbind(U_mat_c, stochastic_basis_funs_c)
combined_basis_uc <- cbind(U_mat_uc, stochastic_basis_funs_uc)

ols_c <- lsfit(x = combined_basis_c, y = x_cent, intercept = FALSE)
ols_uc <- lsfit(x = combined_basis_uc, y = x_cent, intercept = FALSE)

coef_U_c <- ols_c$coefficients[1:2, ]
coef_stochastic_c <- ols_c$coefficients[-c(1:2), ]

fitted_U_c <- U_mat_c %*%  coef_U_c 
fitted_stochastic_c <- stochastic_basis_funs_c %*% coef_stochastic_c

coef_U_uc <- ols_uc$coefficients[1:2, ]
coef_stochastic_uc <- ols_uc$coefficients[-c(1:2), ]

fitted_U_uc <- U_mat_uc %*%  coef_U_uc 
fitted_stochastic_uc <- stochastic_basis_funs_uc %*% coef_stochastic_uc

par(mfrow = c(2, 2))
ylim = range(fitted_stochastic_uc, fitted_stochastic_c, fitted_U_c, fitted_U_uc)
matplot(fitted_U_c, type = "l", ylim = ylim)
matplot(fitted_stochastic_c, type = "l", ylim = ylim)

matplot(fitted_U_uc, type = "l", ylim = ylim)
matplot(fitted_stochastic_uc, type = "l", ylim =ylim)

# Create Publishable Quality Figures: -------------------------------------


## 1) PDA Basis Functions: ------------------------------------------------
basis_funs_dt <- data.table(
  t = rep(milliseconds, 2),
  correction = rep(x = c("Initial", "Final"), 
                   each = length(milliseconds)),
  rbind(combined_basis_uc, combined_basis_c)
)
names(basis_funs_dt)[-c(1:2)] <- c(
  paste0("zero_input_", 1:2),
  paste0("zero_state_", 1:4)
)

basis_funs_dt_lng <- melt.data.table(basis_funs_dt,
                                     id.vars = c("t", "correction"))

basis_funs_dt_lng[, covar_name := stringr::str_extract(variable, "zero_(input|state)")]
basis_funs_dt_lng[, basis_fun_num := stringr::str_extract(variable, "(?<=zero_(input|state)_)\\d")]
basis_funs_dt_lng[, covar_name := factor(covar_name,
                                         levels = c('zero_input', "zero_state"),
                                         labels = c("Zero-input", "Zero-state"))]
basis_funs_dt_lng[, basis_fun_name := paste(covar_name, "Basis Function", basis_fun_num)]

p <- ggplot(data = basis_funs_dt_lng) +
  facet_wrap(correction ~ covar_name, scales = "free_y") +
  aes(x = t, y = value, group = basis_fun_num, colour = basis_fun_name) +
  geom_line() +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  labs(x = "$t$", y = "Basis Function Value") +
  guides(colour = guide_legend(override.aes = list(linewidth = 1)))

tikz(here::here("outputs", "real-data", "paper-plots", "com-pda-basis.tex"),
     width = doc_width_inches,
     height = 1.01 *doc_width_inches,
     standAlone = TRUE)
print(p)  
dev.off()

tinytex::lualatex(here::here("outputs", "real-data", "paper-plots", "com-pda-basis.tex"))



# Variance Decomposition: -------------------------------------------------

fitted_dt_c <- data.table(
  t = milliseconds,
  correction = "Final",
  fitted_stochastic_c,
  fitted_U_c)

fitted_dt_uc <- data.table(
  t = milliseconds,
  correction = "Initial",
  fitted_stochastic_uc,
  fitted_U_uc)

fitted_dt <- rbind(fitted_dt_uc, fitted_dt_c)
names(fitted_dt)[-c(1:2)] <- paste(
  rep(c("stochastic", "ic"), each = N),
  rep(paste0("sub_", seq_len(N)), times = 2), 
  sep = "_")

fitted_dt_lng <- melt.data.table(data = fitted_dt, 
                                 id.vars = c("t", "correction"), 
                                 variable.factor = FALSE, 
                                 value.factor = FALSE)
fitted_dt_lng[, component := stringr::str_extract(variable, "(stochastic|ic)")]
fitted_dt_lng[, subject := stringr::str_extract(variable, "(?<=(stochastic|ic)_)sub_\\d{1,2}")]

subjects_sample <- sample(unique(fitted_dt_lng$subject), size = 20)

fitted_dt_lng[, component := factor(component,
                                    levels = c("ic", "stochastic"),
                                    labels = c("Zero-input (Initial Conditions)", "Zero-state (Stochastic Disturbance)")
                                    )]

fitted_dt_lng[, correction := factor(correction, 
                                     levels = c("Final", "Initial"),
                                     labels = c("Using Final (Bias-Reduced) Estimates", "Using Initial Estimates"))]





p <- ggplot(fitted_dt_lng[subject %in% subjects_sample]) +
  aes(x = t, y = value, group = subject, colour = subject) +
  facet_grid(correction ~ component) +
  geom_line() +
  theme(legend.position = "none") +
  labs(x = "$t$", y = "Contribution to $x(t) - \\mathbb{E}[x(t)]$")

tikz(here::here("outputs", "real-data", "paper-plots", "com-pda-decomposition.tex"),
     width = doc_width_inches,
     height = doc_width_inches * 0.9,
     standAlone = TRUE)
print(p)  
dev.off()

tinytex::lualatex(here::here("outputs", "real-data", "paper-plots", "com-pda-decomposition.tex"))
