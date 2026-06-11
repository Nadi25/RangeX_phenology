

# 03_temp_sens ------------------------------------------------------------

# Temperature sensitivity plotting -------------------------------------------------

# set theme for plots ------------------------------------
# base_size = 20 for presentations
theme_set(theme_bw())


# Source temperature sensitivity analysis script --------------------------
source("Temperature_sensitivity_analysis.R")


# plot temp sens predictions julian days ----------------------------------------------
# for both regions 
# but used separate models
p_sens <- plot_temp_sens_predictions(sens_all_gs, plot_predictions_gs,
                                     expression("Temperature sensitivity (days/"*degree*"C)"))
p_sens

# ggsave(filename = "Output/Sensitivity/Temperature_sensitivity_gs_bud_flower_fruit_seed_onset_NOR_CHE.png", 
#       plot = p_sens,
#        width = 15, height = 10, units = "in")


p_sens3 <- plot_temp_sens_predictions3(sens_all_gs, plot_predictions_gs3,
                                       expression("Temperature sensitivity (days/"*degree*"C)"))

p_sens3

# ggsave(filename = "Output/Sensitivity/Temperature_sensitivity_gs_bud_flower_fruit_seed_onset_NOR_CHE_violin.png", 
#       plot = p_sens3,
#        width = 15, height = 10, units = "in")


# plot temp sens predictions growing degree days (GDD) ----------------------------------------------
# for both regions 
# but used separate models
p_sens_gdd <- plot_temp_sens_predictions(sens_all_gs_gdd, plot_predictions_gs_gdd, 
                                         expression("Temperature sensitivity (GDD2/"*degree*"C)"))
p_sens_gdd

# ggsave(filename = "Output/Sensitivity/Temperature_sensitivity_GDD_gs_bud_flower_fruit_seed_onset_NOR_CHE.png", 
#       plot = p_sens_gdd,
#        width = 15, height = 10, units = "in")


p_sens_gdd3 <- plot_temp_sens_predictions3(sens_all_gs_gdd, plot_predictions_gs_gdd3, 
                                         expression("Temperature sensitivity (GDD2/"*degree*"C)"))
p_sens_gdd3

# ggsave(filename = "Output/Sensitivity/Temperature_sensitivity_GDD_gs_bud_flower_fruit_seed_onset_NOR_CHE_violin.png", 
#       plot = p_sens_gdd3,
#        width = 15, height = 10, units = "in")

















