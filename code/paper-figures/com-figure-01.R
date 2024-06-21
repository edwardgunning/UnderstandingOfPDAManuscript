library(data.table)    # CRAN v1.14.2
library(fda)           # CRAN v5.5.1
library(rgl)           # CRAN v0.108.3
library(scatterplot3d) # CRAN v0.3-41
library(dichromat)     # CRAN v2.0-0
library(ggplot2)       # CRAN v3.4.0
library(tikzDevice)    # CRAN v0.12.3.1

source(here::here("code", "theme_gunning.R"))
theme_gunning()
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937

# -------------------------------------------------------------------------
participant_dataset <- readRDS(file = here::here("data", "real-dataset-4161.rds"))
participant_dataset_com_x <- participant_dataset[side == "right" & location == "CentreOfMass" & plane_of_motion == "z"] 
# participant_dataset_com_x <- participant_dataset[side == "right" & location == "Knee" & plane_of_motion == "fle"] 
participant_dataset_com_x[, .N]

subset_dt <- participant_dataset_com_x[, c("stride_num", paste0("data_", 0:197))]
y <- as.matrix((subset_dt[, `data_0`:data_141]))
y <- (y - mean(y))/1000
sample_hz <- 200
second_per_frame <- 1/sample_hz
frame <- 0:141
seconds <- second_per_frame * frame
milliseconds <- 10 * seconds
quintic_Bspline_basis <- create.bspline.basis(
  rangeval = range(milliseconds), 
  norder = 6,
  nbasis = length(milliseconds))
fd_par <- fdPar(fdobj = quintic_Bspline_basis, Lfdobj = 4, lambda = 10^-2)

lambda_grid <- seq(-12, 10, by = 1)
gcv_vec <- vector(mode = "numeric", length = length(lambda_grid))
for(i in seq_along(lambda_grid)) {
  fd_par <- fdPar(fdobj = quintic_Bspline_basis, Lfdobj = 4, lambda = 10^lambda_grid[i])
  smooth_i <- smooth.basis(argvals = milliseconds, y = t(y), fdParobj = fd_par)
  gcv_vec[i] <- mean(smooth_i$gcv)
}

plot(lambda_grid, gcv_vec)
abline(v = lambda_grid[which.min(gcv_vec)])
(log10_lambda_best <- lambda_grid[which.min(gcv_vec)])



fd_par <- fdPar(fdobj = quintic_Bspline_basis, Lfdobj = 4, lambda = 10^log10_lambda_best)
fd_obj <- smooth.basis(argvals = milliseconds, y = t(y), fdParobj = fd_par)$fd

par(mfrow = c(1, 3))
plot(fd_obj, Lfdobj = 0)
plot(fd_obj, Lfdobj = 1)
plot(fd_obj, Lfdobj = 2)


# -------------------------------------------------------------------------
x <- eval.fd(milliseconds, fdobj = fd_obj, Lfdobj = 0)
Dx <- eval.fd(milliseconds, fdobj = fd_obj, Lfdobj = 1)
D2x <- eval.fd(milliseconds, fdobj = fd_obj, Lfdobj = 2)




# -------------------------------------------------------------------------

xfig_dt <- data.table(t_grid = milliseconds, x)
D2xfig_dt <- data.table(t_grid = milliseconds, D2x)
Dxfig_dt <- data.table(t_grid = milliseconds, Dx)

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



# -------------------------------------------------------------------------
theme_update(plot.title = element_text(size = 8))
(p1 <- ggplot(data = plot_dt[rep_fun %in% paste0("V", 1:20)]) +
   aes(x = t_grid, y = xt, group = rep_fun, colour = rep_fun) +
   geom_path() +
   scale_color_manual(values = seq_len(ncol(x))) +
   theme(legend.position = "none") +
   scale_x_continuous(expand = c(0, 0)) +
   # scale_x_continuous(breaks = c(-1, 0, 1), limits = y_lims_square) +
   # scale_y_continuous(breaks = c(-1, 0, 1), limits = y_lims_square) +
   labs(x = "$t$ (seconds $\\times 10$)", y = "$x(t)$", title = "Displacement"))

(p2 <- ggplot(data = plot_dt[rep_fun %in% paste0("V", 1:20)]) +
    aes(x = t_grid, y = Dxt, group = rep_fun, colour = rep_fun) +
    geom_path() +
    scale_color_manual(values = seq_len(ncol(x))) +
    theme(legend.position = "none") +
    scale_x_continuous(expand = c(0, 0)) +
    # scale_x_continuous(breaks = c(-1, 0, 1), limits = y_lims_square) +
    # scale_y_continuous(breaks = c(-1, 0, 1), limits = y_lims_square) +
    labs(x = "$t$ (seconds $\\times 10$)", y = "$Dx(t)$", title = "Velocity"))


(p3 <- ggplot(data = plot_dt[rep_fun %in% paste0("V", 1:20)]) +
  aes(x = t_grid, y = D2xt, group = rep_fun, colour = rep_fun) +
  geom_path() +
  scale_color_manual(values = seq_len(ncol(x))) +
  theme(legend.position = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  # scale_x_continuous(breaks = c(-1, 0, 1), limits = y_lims_square) +
  # scale_y_continuous(breaks = c(-1, 0, 1), limits = y_lims_square) +
  labs(x = "$t$ (seconds $\\times 10$)", y = "$D^2x(t)$", title = "Acceleration"))


full_plot <- ggpubr::ggarrange(plotlist = list(p1, p2, p3),
                               nrow = 1,
                               ncol = 3, 
                               labels = c("(a)", "(b)", "(c)"), 
                               align = "hv", font.label = list(size = 8))

full_plot


tikz(here::here("outputs", "real-data", "paper-plots", "curves.plot.tex"), 
     width = (1 * doc_width_inches),
     height = (0.33 * doc_width_inches))#,
     #standAlone = TRUE)
full_plot
dev.off()

# tinytex::lualatex(here::here("outputs", "real-data", "paper-plots", "curves.plot.tex"))


