library(data.table)    # CRAN v1.14.2
library(fda)           # CRAN v5.5.1
library(rgl)           # CRAN v0.108.3
library(scatterplot3d) # CRAN v0.3-41
library(dichromat)     # CRAN v2.0-0
library(tikzDevice)  # CRAN v0.12.3.1 
library(data.table)  # CRAN v1.14.2
library(ggplot2)     # CRAN v3.3.5
library(randomcoloR) # CRAN v1.1.0.1

doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937

participant_dataset <- readRDS(file = here::here("data", "real-dataset-4161.rds"))
participant_dataset_com_x <- participant_dataset[side == "right" & location == "CentreOfMass" & plane_of_motion == "z"] 
# participant_dataset_com_x <- participant_dataset[side == "right" & location == "Knee" & plane_of_motion == "fle"] 
participant_dataset_com_x[, .N]

subset_dt <- participant_dataset_com_x[, c("stride_num", paste0("data_", 0:197))]
subset_dt_lng <- melt.data.table(data = subset_dt,
                                 id.vars = c("stride_num"), 
                                 measure.vars = paste0("data_", 0:197),
                                 variable.name = "frame",
                                 value.name = "angle", 
                                 variable.factor = FALSE,
                                 value.factor = FALSE,
                                 na.rm = TRUE,
                                 verbose = TRUE)
subset_dt_lng[, frame := as.numeric(stringr::str_remove(frame, "data_"))]
subset_dt_lng
setorderv(subset_dt_lng, c("stride_num", "frame"))
head(subset_dt_lng)
subset_dt_lng[, frame_new := 0:(.N-1)]


sample_hz <- 200
second_per_frame <- 1/sample_hz
# millisecond_per_frame <- 1000 * second_per_frame
subset_dt_lng[, seconds := second_per_frame * frame_new]

segment_times <- unique(subset_dt_lng[frame == 0, seconds])


plot(subset_dt_lng$seconds[1:2000], subset_dt_lng$angle[1:2000], type = "l")





# -------------------------------------------------------------------------

subset_dt_lng <- subset_dt_lng[frame <= 141]
subset_dt_lng

y <- as.matrix((subset_dt[, `data_0`:data_141]))
y <- (y - mean(y))/1000


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


# Now SHM -----------------------------------------------------------------

# Packages ----------------------------------------------------------------



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
                                    mu_init = c(2, 0), 
                                    intensity = 2,
                                    sigma_init = diag(rep(0.05, 2)))

# saveRDS(object = list(datset_fig = dataset_fig, seed = 1996, time = timestamp()),
#         file = here::here("outputs", "SHM", "simulation-results", "dataset-fig.rds"))

# Reshape data for plotting -----------------------------------------------
xfig <- dataset_fig$x
D2xfig <- calculate_deriv(grid_points = grid_points, x = xfig, norder = 2)
Dxfig <- calculate_deriv(grid_points = grid_points, x = xfig, norder = 1)




tikz(here::here("outputs", "SHM", "paper-plots", "3d-phase-plane.tex"), 
     width = (1 * doc_width_inches),
     height = (0.5 * doc_width_inches)) #,
     # standAlone = TRUE)
par(mfrow = c(1, 2))
xlim <- range(xfig)
ylim <- range(Dxfig)
zlim <- range(D2xfig)
sp <- scatterplot3d(x = xfig[, 1],
                    y = Dxfig[,1],
                    z = D2xfig[,1],
                    type = "l",
                    cex.lab = 1,
                    cex.axis = 0.75,
                    xlim = xlim,
                    ylim = ylim,
                    zlim = zlim,
                    ylab = "",
                    zlab = "$D^2 x(t)$",
                    xlab = "$x(t)$",
                    angle = 20,
                    pch = 20,
                    cex.main = 0.9,
                    main = "(a) Simulated SHM Data")
# par("usr")
text(y = -3.3,  x = 7, "$Dx(t)$", cex = 1, xpd=TRUE, srt = 30)

for(j in seq_len(N)) {
  sp$points3d(x = xfig[, j],
              y = Dxfig[,j],
              z = D2xfig[,j],
              type = "l",
              col = j)
}

xlim <- range(x)
ylim <- range(Dx)
zlim <- range(D2x)
sp <- scatterplot3d(x = x[, 1],
                    y = Dx[,1],
                    z = D2x[,1],
                    type = "l",
                    cex.lab = 1,
                    cex.axis = 0.75,
                    xlim = xlim,
                    ylim = ylim,
                    zlim = zlim,
                    ylab = "",
                    zlab = "$D^2 x(t)$",
                    xlab = "$x(t)$",
                    angle = 20,
                    x.ticklabs = c("", "0.04", "", 0, "", "-0.04", ""),
                    y.ticklabs =  c("-0.1", "", "0", "0.05", ""),
                    z.ticklabs = c('-0.2', '', '-0.1', '', '0', '', "0.1", "", "0.2", "0"),
                    pch = 20,
                    cex.main = 0.9,
                    main = "(b) Running Data")
# par("usr")
text(y = -5,  x = 6, "$Dx(t)$", cex = 1, xpd=TRUE, srt = 30)

for(j in seq_len(10)) {
  sp$points3d(x = x[,j],
              y = Dx[,j],
              z = D2x[,j],
              type = "l",
              col = j)
}

dev.off()

# tinytex::lualatex(here::here("outputs", "SHM", "paper-plots", "3d-phase-plane.tex"))


par(mfrow = c(1, 2))
xlim <- range(xfig)
ylim <- range(Dxfig)
zlim <- range(D2xfig)
sp <- scatterplot3d(x = xfig[, 1],
                    y = Dxfig[,1],
                    z = D2xfig[,1],
                    type = "l",
                    cex.lab = 1,
                    cex.axis = 0.75,
                    xlim = xlim,
                    ylim = ylim,
                    zlim = zlim,
                    ylab = "",
                    zlab = "$D^2 x(t)$",
                    xlab = "$x(t)$",
                    angle = 30,
                    pch = 20,
                    main = "Simulated SHM Data")
# par("usr")
text(y = -3.37,  x = 6, "$Dx(t)$", srt = 35, cex = 1)

for(j in seq_len(N)) {
  sp$points3d(x = xfig[, j],
              y = Dxfig[,j],
              z = D2xfig[,j],
              type = "l",
              col = j)
}

xlim <- range(x)
ylim <- range(Dx)
zlim <- range(D2x)
sp <- scatterplot3d(x = x[, 1],
                    y = Dx[,1],
                    z = D2x[,1],
                    type = "l",
                    cex.lab = 1,
                    cex.axis = 0.75,
                    xlim = xlim,
                    ylim = ylim,
                    zlim = zlim,
                    ylab = "",
                    zlab = "$D^2 x(t)$",
                    xlab = "$x(t)$",
                    angle = 30,
                    pch = 20,
                    main = "Running Data")
# par("usr")
text(y = -3.37,  x = 6, "$Dx(t)$", srt = 35, cex = 1)

for(j in seq_len(10)) {
  sp$points3d(x = x[,j],
              y = Dx[,j],
              z = D2x[,j],
              type = "l",
              col = j)
}

