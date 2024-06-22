# Script to do the data analysis for PDA paper: ---------------------------
library(data.table)    # CRAN v1.14.2
library(fda)           # CRAN v5.5.1

# Read in data: -----------------------------------------------------------
participant_dataset <- readRDS(file = here::here("data", "real-dataset-4161.rds"))
participant_dataset_com_x <- participant_dataset[side == "right" &
                                                   location == "CentreOfMass" &
                                                   plane_of_motion == "z" &
                                                   !(stride_num %in% c(5, 26, 60, 80))]
N <- participant_dataset_com_x[, .N]
N # [1] 82
# Reshape data for a short plot: ------------------------------------------
subset_dt <- participant_dataset_com_x[, c("stride_num", paste0("data_", 0:197))]
subset_dt_lng <- melt.data.table(data = subset_dt,
                                 id.vars = c("stride_num"), 
                                 measure.vars = paste0("data_", 0:197),
                                 variable.name = "frame",
                                 value.name = "com", 
                                 variable.factor = FALSE,
                                 value.factor = FALSE,
                                 na.rm = TRUE,
                                 verbose = TRUE)
subset_dt_lng[, frame := as.numeric(stringr::str_remove(frame, "data_"))]
subset_dt_lng
setorderv(subset_dt_lng, c("stride_num", "frame"))
head(subset_dt_lng)
subset_dt_lng[, max(frame), by = stride_num][, min(V1)]
subset_dt_lng[, frame_new := 0:(.N-1)]
sample_hz <- 200
second_per_frame <- 1/sample_hz
subset_dt_lng[, seconds := second_per_frame * frame_new]
segment_times <- unique(subset_dt_lng[frame == 0, seconds])
plot(x = subset_dt_lng$seconds[1:2000],
     y = subset_dt_lng$com[1:2000], 
     type = "l",
     xlab = "time (secs)",
     ylab = "vertical position")
start_point <- subset_dt_lng$com[1]

subset_dt_lng[, .(max_frame = max(frame)), by = stride_num][, min(max_frame)]
subset_dt_lng[, .(max_frame = max(frame)), by = stride_num][, unique(max_frame)]


# Set up basis for representing the data: ---------------------------------
y <- as.matrix((subset_dt[, `data_0`:data_141])) / 1000
y <- y - mean(y) # center around average point
frame <- 0:141
seconds <- second_per_frame * frame
milliseconds <- 10 * seconds
quintic_Bspline_basis <- create.bspline.basis(
  rangeval = range(milliseconds), 
  norder = 6,
  nbasis = length(milliseconds) + 6 - 2)
par(mfrow = c(1, 2))
matplot(milliseconds, t(y), type = "l")
plot(quintic_Bspline_basis)



# Grid search for sp: -----------------------------------------------------
lambda_grid <- seq(-14, -2, by = 1)
gcv_vec <- vector(mode = "numeric", length = length(lambda_grid))
for(i in seq_along(lambda_grid)) {
  fd_par <- fdPar(fdobj = quintic_Bspline_basis, Lfdobj = 4, lambda = 10^lambda_grid[i])
  smooth_i <- smooth.basis(argvals = milliseconds, y = t(y), fdParobj = fd_par)
  gcv_vec[i] <- mean(smooth_i$gcv)
}
log10_lambda_best <- lambda_grid[which.min(gcv_vec)]
fd_par <- fdPar(fdobj = quintic_Bspline_basis, Lfdobj = 4, lambda = 10^log10_lambda_best)
fd_obj <- smooth.basis(argvals = milliseconds, y = t(y), fdParobj = fd_par)$fd


# Unpack objects on grid: -------------------------------------------------
x <- eval.fd(milliseconds, fdobj = fd_obj, Lfdobj = 0)
Dx <- eval.fd(milliseconds, fdobj = fd_obj, Lfdobj = 1)
D2x <- eval.fd(milliseconds, fdobj = fd_obj, Lfdobj = 2)

saveRDS(object = list(x = x,
                      Dx = Dx,
                      D2x = D2x,
                      lambda = 10^log10_lambda_best,
                      fd_obj = fd_obj,
                      quintic_Bspline_basis = quintic_Bspline_basis,
                      milliseconds = milliseconds,
                      N = N),
        file = here::here("outputs", "real-data", "prepared-data-01.rds"))




