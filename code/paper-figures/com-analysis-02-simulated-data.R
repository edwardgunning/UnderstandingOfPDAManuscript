library(scatterplot3d) # CRAN v0.3-41
library(tikzDevice)    # CRAN v0.12.3.1
doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937

pda_simulation_result <- readRDS(here::here("outputs", "real-data", "pda-simulation-result-02.rds"))
pda_results <- readRDS(file = here::here("outputs", "real-data", "pda-result-02.rds"))
prepared_data <- readRDS(here::here("outputs", "real-data", "prepared-data-02.rds"))

x <- prepared_data$x
Dx <- prepared_data$Dx
D2x <- prepared_data$D2x

milliseconds <- prepared_data$milliseconds

x_i <- pda_simulation_result$x_array[,,1]
Dx_i <- pda_simulation_result$Dx_array[,,1]
D2x_i <- pda_simulation_result$D2x_array[,,1]

(N <- prepared_data$N)


tikz(here::here("outputs", "real-data", "paper-plots", "real-vs-simulated-3d-phase-plane-02.tex"), 
     width = (1 * doc_width_inches),
     height = (0.5 * doc_width_inches),
     standAlone = TRUE)
par("usr")
par(mfrow = c(1, 2))
xlim <- range(c(x, x_i))
ylim <- range(c(Dx, Dx_i))
zlim <- range(c(D2x, D2x_i))
sp <- scatterplot3d(x = x[, 1],
                    y = Dx[,1],
                    z = D2x[,1],
                    type = "l",
                    cex.lab = 1,
                    cex.axis = 0.75,
                    cex.main = 0.9,
                    xlim = xlim,
                    ylim = ylim,
                    zlim = zlim,
                    ylab = "",
                    zlab = "$D^2 x(t)$",
                    xlab = "$x(t)$",
                    angle = 30,
                    pch = 20,
                    label.tick.marks = FALSE,
                    tick.marks = FALSE,
                    main = "(a) Real Running Data")


# par("usr")
# par("usr")
text(y = -4,  x = 6, "$Dx(t)$", srt = 35, cex = 1, xpd = TRUE)

for(j in seq_len(N)) {
  sp$points3d(x = x[,j],
              y = Dx[,j],
              z = D2x[,j],
              type = "l",
              col = j)
}



sp <- scatterplot3d(x = x_i[, 1],
                    y = Dx_i[,1],
                    z = D2x_i[,1],
                    type = "l",
                    cex.lab = 1,
                    cex.main = 0.9,
                    cex.axis = 0.75,
                    xlim = xlim,
                    ylim = ylim,
                    zlim = zlim,
                    ylab = "",
                    zlab = "$D^2 x(t)$",
                    xlab = "$x(t)$",
                    angle = 30,
                    pch = 20,
                    label.tick.marks = FALSE,
                    tick.marks = FALSE,
                    main = "(b) Simulated Running Data")
# par("usr")
text(y = -4,  x = 6, "$Dx(t)$", srt = 35, cex = 1, xpd = TRUE)

for(j in seq_len(N)) {
  sp$points3d(x = x_i[,j],
              y = Dx_i[,j],
              z = D2x_i[,j],
              type = "l",
              col = j)
}
dev.off()
tinytex::lualatex(
  here::here("outputs", "real-data", "paper-plots", "real-vs-simulated-3d-phase-plane-02.tex"))

