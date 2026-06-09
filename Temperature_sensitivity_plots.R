

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



# violins for raw ---------------------------------------------------------


p_sens2 <- ggplot() +
  
  # raw species values
  geom_violin(
    data = sens_all_gs,
    aes(
      x = treat_competition,
      y = temp_sens,
      fill = treat_competition
    ),
    position = position_dodge(width = 0.9),
    alpha = 0.25,
    color = NA,
    trim = FALSE
  ) +
  
  # model predictions
  geom_point(
    data = plot_predictions_gs,
    aes(
      x = treat_competition,
      y = temp_sens,
      color = treat_competition
    ),
    size = 4, stroke = 1.2
  ) +
  
  geom_errorbar(
    data = plot_predictions_gs,
    aes(
      x = treat_competition,
      ymin = plo,
      ymax = phi,
      color = treat_competition
    ),
    width = 0.12
  ) +
  
  facet_grid(region ~ stage) +
  
  scale_color_manual(values = c(
    "with" = "#528B8B",
    "without" = "#CD950C"
  )) +
  scale_fill_manual(values = c(
    "with" = "#528B8B",
    "without" = "#CD950C"
  )) +
  labs(
    x = "Biotic interactions",
    y =  expression("Temperature sensitivity (days/"*degree*"C)")
  ) +
  theme(
    legend.position = "none"
  )+
  geom_hline(yintercept=0, linetype = "dashed")
p_sens2

# ggsave(filename = "Output/Sensitivity/Temperature_sensitivity_gs_bud_flower_fruit_seed_onset_NOR_CHE_violin.png", 
#       plot = p_sens2,
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


p_sens_gdd2 <- plot_temp_sens_predictions(sens_all_gs_gdd, plot_predictions_gs_gdd2, 
                                         expression("Temperature sensitivity (GDD2/"*degree*"C)"))
p_sens_gdd2



















