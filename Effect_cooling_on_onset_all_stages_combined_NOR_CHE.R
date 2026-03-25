
# Effect of transplantation on Budding, Flowering, Fruiting onset with and without competition ---------------------------------------

## Date:      25.03.26
## Author:    Nadine Arzt
## Purpose:   Effect of transplantation on Budding, Flowering, Fruiting onset NOR and CHE



library(patchwork)


source("Effect_cooling_on_budding_onset_NOR_CHE.R")

source("Effect_cooling_on_flowering_onset_NOR_CHE.R")

source("Effect_cooling_on_fruiting_onset_NOR_CHE.R")


# the three plots
#budding
nor_che_delta_raw_cool_bud <- nor_che_delta_raw_cool_bud +
  guides(color = "none", shape = "none") +
  labs(title = "Effect transplantation on:\nBudding onset")

# flowering
nor_che_delta_raw_cool <- nor_che_delta_raw_cool +
  guides(color = "none", shape = "none")+
  labs(title = "Flowering onset")

# fruiting
# dont remove legend here to have it in panel
nor_che_delta_raw_cool_fruit <- nor_che_delta_raw_cool_fruit +
  labs(title = "Fruiting onset")



# combine the three plots in a panel --------------------------------------
panel <- (nor_che_delta_raw_cool_bud + nor_che_delta_raw_cool + nor_che_delta_raw_cool_fruit) +
  plot_layout(ncol = 3, guides = "collect")&
  theme(legend.position = "bottom")
panel


# Add labels A, B, C
panel_labeled <- panel + 
  plot_annotation(
    tag_levels = "A",   # produces A, B, C
    tag_prefix = "",    # no extra characters
    tag_suffix = ""     # no punctuation
  )
panel_labeled



ggsave(filename = "Output/Onset/Effect_cooling_bud_flower_fruit_NOR_CHE.png", 
       plot = panel_labeled, width = 20, height = 12, units = "in")





