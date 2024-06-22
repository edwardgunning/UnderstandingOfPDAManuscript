# ------------------------------------------------------------------------#
# Create Plot of PDA Estimates.
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


# Load in Results and Unpack: ---------------------------------------------
pda_results <- readRDS(file = here::here("outputs", "real-data", "pda-result-01.rds"))
prepared_data <- readRDS(here::here("outputs", "real-data", "prepared-data-01.rds"))
beta_hat <- pda_results$pda_result$beta
milliseconds <- prepared_data$milliseconds


# Reshape Data: -----------------------------------------------------------
param_dt <- data.table(
  t = rep(milliseconds, times = 3),
  parameter = rep(c("$\\alpha(t)$", "$\\beta_0(t)$", "$\\beta_1(t)$"),
                  each = length(milliseconds)),
  rbind(beta_hat[,1,c(1, 10, 11)],
        beta_hat[,2,c(1, 10, 11)],
        beta_hat[,3,c(1, 10, 11)])
  )
names(param_dt)[-c(1:2)] <- c("est_1", "est_10", "est_11")

param_dt_lng <- melt.data.table(param_dt, id.vars = c("t", "parameter"))
param_dt_lng[, variable := factor(variable,
                                  levels = c("est_1", "est_10", "est_11"),
                                  labels = c("Initial", "$9$th Iteration", "$10$th Iteration (Final)"))]

(p <- ggplot(data = param_dt_lng) +
      aes(x = t, y = value, group = variable, colour = variable) +
      facet_wrap(~ parameter, scales = "free_y") +
      geom_line() +
      labs(x = "$t$", y = "$\\beta(t)$", colour = "Estimate:") +
      theme(legend.position = "bottom") +
      guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 0.75))) +
      # scale_color_manual(values = rev(c("#CC79A7", "#009E73","#4E84C4"))) +
      geom_hline(yintercept = 0, lty = 2, col = "darkgrey"))

tikz(here::here("outputs", "real-data", "paper-plots", "pda-estimates.tex"), 
     width = (1 * doc_width_inches),
     height = (0.4 * doc_width_inches),
     standAlone = TRUE)
p
dev.off()



tinytex::lualatex(here::here("outputs", "real-data", "paper-plots", "pda-estimates.tex"))


# Change font to Beamer for Giles' JSM Presentation: ----------------------
options(tikzLatexPackages = c(getOption( "tikzLatexPackages"),
                              "\\renewcommand{\\familydefault}{\\sfdefault}"))
tikz(here::here("outputs", "real-data", "rough", "pda-estimates.tex"), 
     width = (1 * doc_width_inches),
     height = (0.4 * doc_width_inches),
     standAlone = TRUE)
p
dev.off()
tinytex::lualatex(here::here("outputs", "real-data", "rough", "pda-estimates.tex"))

