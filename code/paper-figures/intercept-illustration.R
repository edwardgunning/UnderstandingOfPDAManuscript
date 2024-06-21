# -------------------------------------------------------------------------
# Make illustrative figure to show why we should include an intercept when
# approximating an unknown non-linear function.
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Packages:
library(ggplot2)    # CRAN v3.3.5
library(tikzDevice) # CRAN v0.12.3.1


# -------------------------------------------------------------------------
# Simulate some data:
set.seed(1996)
intercept1 <- slope1 <- slope2 <- vector(mode = "numeric", length = 1000)
for(i in 1:1000) {
  x <- rnorm(n = 10000, mean = 0.6, sd = 0.1)
  y <- x^2 + rnorm(n = 1000, mean = 0, sd = 0.2)
  
  lm1 <- lm(y ~ x)
  slope1[i] <- coef(lm1)[2]
  intercept1[i] <- coef(lm1)[1]
  
  lm2 <- lm(y ~ 0 + x)
  slope2[i] <- coef(lm2)[1]
}

x_seq <- seq(0, 1, by = 0.005)
density_df <- data.frame(x = x_seq, y = 0.05 * dnorm(x_seq, 0.6, 0.1))
plot_df <- data.frame(x = x_seq, 
                      fx_unknown = x_seq^2,
                      taylor = -(0.6^2) + 1.2 * x_seq,
                      linear_int =mean(intercept1) + mean(slope1) * x_seq,
                      linear_noint = mean(slope2) * x_seq)
plot_df_lng <- reshape2::melt(data = plot_df, 
               id.vars = "x",
               variable.name = "fun", variable.factor = FALSE,
               measure.vars= c("fx_unknown", "taylor", "linear_int", "linear_noint"))
plot_df_lng$fun <- factor(plot_df_lng$fun, levels = c("fx_unknown", "taylor", "linear_int", "linear_noint"),
                          labels = c("Non-linear function $g(x)$",
                                     "Taylor approximation",
                                     "Least-squares approximation",
                                     "Least-squares approximation \\emph{without} intercept"))

p <- ggplot(data = plot_df_lng) +
  #xlim(0, 1) +
  geom_line(aes(x=x, y = value, group = fun, colour = fun), size = 0.85) +
  geom_area(data = density_df, mapping = aes(x  = x, y = y), alpha = 0.5, fill = "darkgrey", colour = "darkgrey") +
  theme_bw() +
  geom_point(pch = 1, aes(x=0,y=0), size = 2) +
  # ylim(c(0, 1)) +
  scale_colour_manual(values = c("black","#1B9E77", "#D95F02", "#7570B3")) +
  theme(legend.position = c(0.38, 0.81),
        text = element_text(size = 8),
        axis.text = element_blank(),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.title = element_blank(),
        legend.text = element_text(),
        axis.ticks = element_blank()) +
  scale_x_continuous(expand = c(0.05, 0.05), limits = c(0, 1))+
  scale_y_continuous(expand = c(0.05, 0.05), limits = c(0, 1))+
  annotate(geom = "text", x = 0.8, y = 0.3,# fill= "lightgrey", alpha = 0.1,
           label = "Distribution of covariate $x$", size = 7 * (5/14)) +
  geom_segment(mapping = aes(x = 0.75, y = 0.25, xend = 0.675, yend = 0.175),
               arrow = arrow(length = unit(0.3,"cm"), type = "closed"),
               arrow.fill = "lightgrey", colour = "darkgrey") +
  guides(colour = guide_legend(override.aes = list(size = 1.5)))
   p

doc_width_cm <- 16
doc_width_inches <- doc_width_cm *  0.3937
tikz(here::here("outputs", "VDP", "paper-plots", "intercept-demo.tex"), 
     width = (0.5 * doc_width_inches),
     height = (0.4 * doc_width_inches))
print(p)
dev.off()


