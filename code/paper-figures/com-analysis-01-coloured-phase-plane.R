library(tikzDevice)    # CRAN v0.12.3.1
library(scatterplot3d) # CRAN v0.3-41
prepared_data <- readRDS(here::here("outputs", "real-data", "prepared-data-01.rds"))
x <- prepared_data$x
Dx <- prepared_data$Dx
D2x <- prepared_data$D2x
milliseconds <- prepared_data$milliseconds
N <- prepared_data$N

x_mean <- apply(x, 1, mean)
Dx_mean <- apply(Dx, 1, mean)
D2x_mean <- apply(D2x, 1, mean)


xlim <- range(c(x_mean))
ylim <- range(c(Dx_mean))
zlim <- range(c(D2x_mean))

myColorRamp <- function(colors, values) {
  v <- (values - min(values))/diff(range(values))
  x <- colorRamp(colors)(v)
  rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
}
cols <- myColorRamp(c("red", "blue"), milliseconds) 

(0.75/0.85) * 0.75
tikz(here::here("outputs", "real-data", "paper-plots", "real-data-3d-phase-plane-coloured.tex"), 
     width = (0.75 * doc_width_inches),
     height = (0.66 * doc_width_inches),
     standAlone = TRUE)
par(mfrow = c(1,1), xpd = TRUE)
sp <- scatterplot3d(x = x_mean,
                    y = Dx_mean,
                    z = D2x_mean,
                    type = "p",
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
                    color = cols,
                    main = "Running Data",
                    xpd= TRUE)
sp$points3d(x = x_mean,
            y = Dx_mean,
            z = D2x_mean,
            type = "l",
            col = "grey")
text(y = -4,  x = 6, "$Dx(t)$", srt = 35, cex = 1, xpd = TRUE)
col_inds <- round(seq(1, length(milliseconds), length.out = 8))
legend(6, 4.2, legend = milliseconds[col_inds], pch = 20, col = cols[col_inds], title = "$t$:", xpd=TRUE)
dev.off()
tinytex::lualatex(here::here("outputs", "real-data", "paper-plots", "real-data-3d-phase-plane-coloured.tex"))
