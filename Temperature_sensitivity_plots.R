

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
p_sens <- plot_temp_sens_predictions(sens_all_gs, plot_predictions_gs)
p_sens

# ggsave(filename = "Output/Sensitivity/Temperature_sensitivity_gs_bud_flower_fruit_seed_onset_NOR_CHE.png", 
#       plot = p_sens,
#        width = 15, height = 10, units = "in")



























