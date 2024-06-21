# -------------------------------------------------------------------------
library(ggplot2)    # CRAN v3.3.5 
library(tikzDevice) # CRAN v0.12.3.1
library(data.table) # CRAN v1.14.2

# -------------------------------------------------------------------------
cov_xs_xt <- readRDS(file = here::here("outputs",
                                       "SHM",
                                       "simulation-results",
                                       "cov_xs_xt_new.rds"))
dataset_fig <- readRDS(here::here("outputs", 
                                  "SHM",
                                  "simulation-results", 
                                  "dataset-fig.rds"))

response <- dataset_fig$datset_fig$x
# -------------------------------------------------------------------------
n_grid <- 100
grid_range <- c(0, 2 * pi)
grid_points <- seq(from = grid_range[1],
                   to = grid_range[2],
                   length.out = n_grid)
M <- plot3D::mesh(grid_points, grid_points)
test_grid <- expand.grid(grid_points, grid_points)

# -------------------------------------------------------------------------
zero_state_covariance_df <- cbind(test_grid, sapply(cov_xs_xt, function(x) {x[1,1]}))
names(zero_state_covariance_df) <- c("s", "t", "cov")
zero_state_covariance_df_wide <- reshape2::dcast(data = zero_state_covariance_df,
                                                 formula = s ~ t,
                                                 value.var = "cov")
stopifnot(zero_state_covariance_df_wide$s == grid_points)
stopifnot(colnames(zero_state_covariance_df_wide)[-1] ==grid_points)
zero_state_covariance_matrix_wide <- as.matrix(zero_state_covariance_df_wide[, -1])

zero_state_covariance_matrix_sym <- (zero_state_covariance_matrix_wide + t(zero_state_covariance_matrix_wide)) / 2
filled.contour(zero_state_covariance_matrix_wide)
zero_state_eigen <- eigen(zero_state_covariance_matrix_sym)
plot(zero_state_eigen$values)


# -------------------------------------------------------------------------
matplot(zero_state_eigen$vectors[, 1:4], type = "l", lty = 1)

zero_state_basis <- zero_state_eigen$vectors[, 1:4]
zero_input_basis <- cbind(cos(grid_points), sin(grid_points))
combined_basis <- cbind(zero_input_basis, zero_state_basis)
matplot(combined_basis, type = "l")
bspline_basis <- fda::eval.basis(evalarg = grid_points, 
                                 basisobj = fda::create.bspline.basis(rangeval = grid_range,
                                                                      nbasis = 4, 
                                                                      norder = 4))
combined_basis_bspline <- cbind(zero_input_basis, bspline_basis)


lsfit_U <- lsfit(y = response, x = zero_input_basis, intercept = FALSE)
lsfit_UV <- lsfit(y = response, x = combined_basis, intercept = FALSE)
lsfit_Ubspline <- lsfit(y = response, x = combined_basis_bspline, intercept = FALSE)

pca_basis <- prcomp(x = t(lsfit_U$residuals), center = FALSE)$rotation[, 1:4]
combined_basis_pca <- cbind(zero_input_basis, pca_basis)

lsfit_Upca <- lsfit(y = response, x = combined_basis_pca, intercept = FALSE)

# -------------------------------------------------------------------------
# Calculate point-wise R-Squared for model without intercept.
# (denominator = E[x^2] rather than var[x] because no constant term there).
Rsq_t_U <- 1 - (apply(lsfit_U$residuals^2, 1, mean) / apply(response^2, 1, mean))
Rsq_t_UV <- 1 - (apply(lsfit_UV$residuals^2, 1, mean) / apply(response^2, 1, mean))
Rsq_t_Ubspline <- 1 - (apply(lsfit_Ubspline$residuals^2, 1, mean) / apply(response^2, 1, mean))
Rsq_t_Upca <- 1 - (apply(lsfit_Upca$residuals^2, 1, mean) / apply(response^2, 1, mean))


# Rsq_int_U <- (1 / (2 * pi)) * sum(0.063467 * Rsq_t_U)
# Rsq_int_UV <- (1 / (2 * pi)) * sum(0.063467 * Rsq_t_UV)
# Rsq_int_Upca <- (1 / (2 * pi)) * sum(0.063467 * Rsq_t_Upca)
# Rsq_int_U <- (1 / (2 * pi)) * sum(0.063467 * Rsq_t_UV)


Rsq_int_U <- 100 * (integrate(f = approxfun(x = grid_points, y = Rsq_t_U, method = "linear"),
          lower = 0, upper = 2 * pi)$value / (2*pi))
Rsq_int_UV <- 100 * (integrate(f = approxfun(x = grid_points, y = Rsq_t_UV, method = "linear"),
                              lower = 0, upper = 2 * pi)$value / (2*pi))

Rsq_int_Ubspline <- 100 * (integrate(f = approxfun(x = grid_points, y = Rsq_t_Ubspline, method = "linear"),
                              lower = 0, upper = 2 * pi)$value / (2*pi))
Rsq_int_Upca <- 100 * (integrate(f = approxfun(x = grid_points, y = Rsq_t_Upca, method = "linear"),
                               lower = 0, upper = 2 * pi)$value / (2*pi))






# -------------------------------------------------------------------------







resid_dt <- data.table(t_grid = rep(grid_points, times = 4),
                       method = factor(rep(c("$U(t)$", "$U(t)$ + Zero-state basis","$U(t)$ + B-spline basis", "$U(t)$ + PCA basis"), each = 100),
                                       levels = c("$U(t)$", "$U(t)$ + Zero-state basis","$U(t)$ + B-spline basis", "$U(t)$ + PCA basis")),
           rbind(resid(lsfit_U), resid(lsfit_UV), resid(lsfit_Ubspline), resid(lsfit_Upca)))
resid_dt_lng <- melt.data.table(data = resid_dt,
                                id.vars = c("t_grid", "method"),
                                measure.vars = paste0("Y", 1:20),
                                variable.name = "subject",
                                value.name = "resid",
                                variable.factor = FALSE,
                                value.factor = FALSE)

resid_dt_lng[, subject := stringr::str_remove(string = subject, pattern = "Y")]

labels_df <- data.frame(
  method = factor(c("$U(t)$", "$U(t)$ + Zero-state basis", "$U(t)$ + B-spline basis", 
    "$U(t)$ + PCA basis"), levels = levels(resid_dt_lng$method)),
  pc_var = paste0("$R^2_{int} = ", round(c(Rsq_int_U,
                                                       Rsq_int_UV,
                                                       Rsq_int_Ubspline,
                                                       Rsq_int_Upca), 2),
                  "\\%$")
)


(upper_panel <- ggplot(data = resid_dt_lng) +
  facet_wrap( ~ method, ncol = 4, nrow = 1) +
  aes(x = t_grid, y = resid, group = subject, color = subject) +
  geom_line() +
  theme_bw() +
  theme(legend.position =  "none",
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 9),
        plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
        strip.text = element_text(size = 9)) +
  ggtitle("Residual Functions") +
  xlab("$t$") +
  geom_label(data = labels_df, inherit.aes = FALSE,
             label.r = unit(0.4, units = "lines"),
             aes(x = pi, y =  0.5, label = pc_var),
             size = 2.15,
             label.padding = unit(0.4, units = "lines")) +
  ylab("$x(t) - \\hat{x}(t)$") +
  scale_x_continuous(breaks = c(0, pi, 2 * pi),
                     labels = c("$0$", "$\\pi$", "$2 \\pi$"),
                     limits = c(0, 2 * pi), expand = rep(0.025, 2)) +
  scale_y_continuous(breaks = seq(-0.2, 0.6, by = 0.2), labels = paste0(c(-0.2, 0, 0.2, 0.4, 0.6))))



ic_dt <- data.table(method = factor(rep(c("$U(t)$", "$U(t)$ + Zero-state basis","$U(t)$ + B-spline basis", "$U(t)$ + PCA basis"), each = 20),
                           levels = c("$U(t)$", "$U(t)$ + Zero-state basis","$U(t)$ + B-spline basis", "$U(t)$ + PCA basis")),
           x0 = rep(response[1, ], 4),
           coef_val = c(lsfit_U$coefficients[1, ],
             lsfit_UV$coefficients[1, ],
             lsfit_Ubspline$coefficients[1, ],
             lsfit_Upca$coefficients[1,]))

lower_panel <- ggplot(data = ic_dt) +
  facet_wrap( ~ method, ncol = 4, nrow = 1) +
  theme_bw() +
  theme(legend.position =  "none",
        axis.title = element_text(size = 9),
        plot.title = element_text(size = 10, hjust = 0.5, face = "bold"),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 9)) + 
  aes(x = x0, y = coef_val) + 
  geom_point(size = 0.75) +
  geom_abline(slope = 1, intercept = 0) +
  ylab("Coefficient of $\\cos(t)$") +
  xlab("Initial Condition $x(0)$") +
  ggtitle("Initial Conditions") +
  scale_y_continuous(breaks = c(0.75, 1, 1.25, 1.5),
                     labels = paste(c(0.75, 1, 1.25, 1.5))) +
  scale_x_continuous(breaks = c(0.8, 0.9, 1, 1.1, 1.2, 1.3), 
                     labels = paste0(c(0.8, 0.9, 1, 1.1, 1.2, 1.3)))




doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937
tikz(here::here("outputs", "SHM", "paper-plots", "shm-bases-fits.tex"),
     width = doc_width_inches, 
     height =  doc_width_inches/1.35) #,
     # standAlone = TRUE)
ggpubr::ggarrange(upper_panel,lower_panel, nrow = 2)
dev.off()

# tinytex::lualatex(here::here("outputs", "SHM", "paper-plots", "shm-bases-fits.tex"))


doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937
tikz(here::here("outputs", "SHM", "paper-plots", "shm-bases-fits.tex"),
     width = doc_width_inches, 
     height =  doc_width_inches/1.35,
standAlone = TRUE)
ggpubr::ggarrange(upper_panel,lower_panel, nrow = 2)
dev.off()

tinytex::lualatex(here::here("outputs", "SHM", "paper-plots", "shm-bases-fits.tex"))






