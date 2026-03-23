

# Combine all four plots --------------------------------------------------

# obviously need to run the other script first to get
p_nb
p_nb_warm
p_bio_nb_species
p_bio_nb_species_warm

library(patchwork)

final_plot <- 
  ((p_nb + labs(title = "Transplantation beyond current range")) | 
     (p_bio_nb_species + labs(title = "adjusted for biomass"))) /
  ((p_nb_warm + labs(title = "Warming at beyond range site")) | 
     (p_bio_nb_species_warm + labs(title = "adjusted for biomass"))) +
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    title = "Effects of temperature and competition on flower number NOR"
  )

final_plot


final_plot <- 
  final_plot +
  plot_layout(guides = "collect") & 
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

final_plot

# ggsave(filename = "Output/Biomass/Cooling_warming_with_without_biomass_NOR.png", 
#        plot = final_plot, 
#        width = 20, height = 18, units = "in")









