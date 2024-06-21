# Packages needed: --------------------------------------------------------
library(ggplot2)    # CRAN v3.3.5
library(data.table) # CRAN v1.14.2
library(tikzDevice) # CRAN v0.12.3.1

# Some initial parameters: ------------------------------------------------
grid_range <- c(0, 2 * pi)
n_grid <- 100
grid_points <- seq(from = grid_range[1],
                   to = grid_range[2],
                   length.out = n_grid)
(delta_t <- unique(round(diff(grid_points), 8)))
# doc_width_cm <- 16
# doc_width_inches <- doc_width_cm *  0.3937

# Read in results of Simulation: ------------------------------------------
simulation_result_list_baseline <- readRDS(here::here("outputs",
                                                      "SHM",
                                                      "simulation-results",
                                                      "simulation-section-4.2-updated.rds"))


results_array_baseline <- array(unlist(simulation_result_list_baseline$results_list), dim = c(n_grid, 4, 500))
error_array_baseline <- results_array_baseline + 1

options(scipen = 999)
integrated_squared_error_baseline <- apply(error_array_baseline, c(2,3), function(x) {
  sum(x^2) * delta_t
})
(mise_baseline <- round(apply(integrated_squared_error_baseline, 1, mean), 4))
(mise_se_baseline <- round(sqrt(1/500) * apply(integrated_squared_error_baseline, 1, sd), 4))


# Extension 01 ------------------------------------------------------------


simulation_result_list_extension_01 <- readRDS(here::here("outputs",
                                                          "SHM",
                                                          "simulation-results",
                                                          "simulation-section-4.2-extension-ic.rds"))

settings_extension_01 <- simulation_result_list_extension_01$settings_dt
results_array_extension_01 <- simulation_result_list_extension_01$results_array
error_array_extension_01 <- results_array_extension_01 + 1

integrated_squared_error_extension_01_low  <- apply(error_array_extension_01[,,settings_extension_01$ic_setting_level==1], c(2,3), function(x) {
  sum(x^2) * delta_t
})
(mise_extension_01_low<- round(apply(integrated_squared_error_extension_01_low, 1, mean), 4))
(mise_extension_01_low <- round(sqrt(1/200) * apply(integrated_squared_error_extension_01_low, 1, sd), 4))

integrated_squared_error_extension_01_high  <- apply(error_array_extension_01[,,settings_extension_01$ic_setting_level==2], c(2,3), function(x) {
  sum(x^2) * delta_t
})
(mise_extension_01_high<- round(apply(integrated_squared_error_extension_01_high, 1, mean), 4))
(mise_extension_01_high <- round(sqrt(1/200) * apply(integrated_squared_error_extension_01_high, 1, sd), 4))




# Extension 02 ------------------------------------------------------------

simulation_result_list_extension_02 <- readRDS(here::here("outputs",
                                                          "SHM",
                                                          "simulation-results",
                                                          "simulation-section-4.2-extension-intensity.rds"))

settings_extension_02 <- simulation_result_list_extension_02$settings_dt
results_array_extension_02 <- simulation_result_list_extension_02$results_array
error_array_extension_02 <- results_array_extension_02 + 1

integrated_squared_error_extension_02_low  <- apply(error_array_extension_02[,,settings_extension_02$intensity_setting_level==1], c(2,3), function(x) {
  sum(x^2) * delta_t
})
(mise_extension_02_low<- round(apply(integrated_squared_error_extension_02_low, 1, mean), 4))
(mise_extension_02_low <- round(sqrt(1/200) * apply(integrated_squared_error_extension_02_low, 1, sd), 4))

integrated_squared_error_extension_02_high  <- apply(error_array_extension_02[,,settings_extension_02$intensity_setting_level==3], c(2,3), function(x) {
  sum(x^2) * delta_t
})
(mise_extension_02_high<- round(apply(integrated_squared_error_extension_02_high, 1, mean), 4))
(mise_extension_02_high <- round(sqrt(1/200) * apply(integrated_squared_error_extension_02_high, 1, sd), 4))



# -------------------------------------------------------------------------
# Extension 02 ------------------------------------------------------------

simulation_result_list_extension_03 <- readRDS(here::here("outputs",
                                                          "SHM",
                                                          "simulation-results",
                                                          "simulation-section-4.2-extension-forcing-variance.rds"))

settings_extension_03 <- simulation_result_list_extension_03$settings_dt
results_array_extension_03 <- simulation_result_list_extension_03$results_array
error_array_extension_03 <- results_array_extension_03 + 1

integrated_squared_error_extension_03_low  <- apply(error_array_extension_03[,,settings_extension_03$sigma_forcing_level==0.15], c(2,3), function(x) {
  sum(x^2) * delta_t
})
(mise_extension_03_low<- round(apply(integrated_squared_error_extension_03_low, 1, mean), 4))
(mise_extension_03_low <- round(sqrt(1/200) * apply(integrated_squared_error_extension_03_low, 1, sd), 4))

integrated_squared_error_extension_03_high  <- apply(error_array_extension_03[,,settings_extension_03$sigma_forcing_level==0.4], c(2,3), function(x) {
  sum(x^2) * delta_t
})
(mise_extension_03_high<- round(apply(integrated_squared_error_extension_03_high, 1, mean), 4))
(mise_extension_03_high <- round(sqrt(1/200) * apply(integrated_squared_error_extension_03_high, 1, sd), 4))

