
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



# ggsave(filename = "Output/Onset/Effect_cooling_bud_flower_fruit_NOR_CHE.png", 
#        plot = panel_labeled, width = 20, height = 12, units = "in")



# One figure transplantation all stages -------------------------------------------
contrast_df_bud_cool$stage <- "Budding"
contrast_df_n_c_cool$stage <- "Flowering"
contrast_df_fruit_cool$stage <- "Fruiting"

contr_all_cool <- bind_rows(contrast_df_bud_cool, contrast_df_n_c_cool, contrast_df_fruit_cool)

contr_all_cool$stage <- factor(contr_all_cool$stage,
                               levels = c("Fruiting", "Flowering", "Budding"))

contr_all_cool <- contr_all_cool |>
  mutate(sig = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE ~ ""
  ))

theme_set(theme_bw(base_size = 20))


onset_all_stages_cooling <- ggplot(contr_all_cool,
                                   aes(x = estimate, y = stage, color = treat_competition, shape = region)) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  geom_point(position = position_dodge(width = 0.5), size = 4) +
  
  geom_errorbar(
    aes(xmin = lower.CL, xmax = upper.CL),
    position = position_dodge(width = 0.5),
    width = 0.2,
    orientation = "y" 
  ) +
  
  facet_wrap(~ region) +
  
  scale_color_manual(values = c("#528B8B", "#CD950C"))+
  
  scale_shape_manual(values = c("Norway" = 16, "Switzerland" = 17))+
  guides(shape = "none")+
  
  labs(title = "Effect of transplantation",
    x = "Δ days shifted onset (high − low)",
    y = "Phenological stage",
    color = "Biotic interactions"
  )+
  
  geom_text(aes(label = sig),
            vjust = -1, position = position_dodge(width = 0.7),
            show.legend = FALSE)
onset_all_stages_cooling


# ggsave(filename = "Output/Onset/Effect_cooling_bud_flower_fruit_NOR_CHE_all_in_one.png", 
#        plot = onset_all_stages_cooling, width = 15, height = 12, units = "in")



# Warming -----------------------------------------------------------------

source("Effect_warming_on_budding_onset_NOR_CHE.R")

source("Effect_warming_on_flower_onset_NOR_CHE.R")

source("Effect_warming_on_fruiting_onset_NOR_CHE.R")


# the three plots
#budding
nor_che_delta_raw_warm_bud <- nor_che_delta_raw_warm_bud +
  guides(color = "none", shape = "none") +
  labs(title = "Effect warming on:\nBudding onset")

# flowering
nor_che_delta_raw_warm_flower <- nor_che_delta_raw_warm_flower +
  guides(color = "none", shape = "none")+
  labs(title = "Flowering onset")

# fruiting
# dont remove legend here to have it in panel
nor_che_delta_raw_warm_fruit <- nor_che_delta_raw_warm_fruit +
  labs(title = "Fruiting onset")


# combine the three plots in a panel --------------------------------------
panel_warm <- (nor_che_delta_raw_warm_bud + nor_che_delta_raw_warm_flower + 
                 nor_che_delta_raw_warm_fruit) +
  plot_layout(ncol = 3, guides = "collect")&
  theme(legend.position = "bottom")
panel_warm


# Add labels D, E, F
panel_warm_labeled <- panel_warm + 
  plot_annotation(
    tag_levels = list(c("D", "E", "F")),
    tag_prefix = "",    # no extra characters
    tag_suffix = "",    # no punctuation
    tag = "D" 
  )
panel_warm_labeled



# ggsave(filename = "Output/Onset/Effect_warming_bud_flower_fruit_NOR_CHE.png", 
#        plot = panel_warm_labeled, width = 20, height = 12, units = "in")




# One figure warming all stages -------------------------------------------
contr_df_warm_bud$stage <- "Budding"
contr_df_warm_flo$stage <- "Flowering"
contr_df_warm_fruit$stage <- "Fruiting"

contr_all_warm <- bind_rows(contr_df_warm_bud, contr_df_warm_flo, contr_df_warm_fruit)

contr_all_warm$stage <- factor(contr_all_warm$stage,
                          levels = c("Fruiting", "Flowering", "Budding"))

contr_all_warm <- contr_all_warm |>
  mutate(sig = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE ~ ""
  ))

theme_set(theme_bw(base_size = 20))


onset_all_stages_warming <- ggplot(contr_all_warm,
       aes(x = estimate, y = stage, color = treat_competition, shape = region)) +
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  geom_point(position = position_dodge(width = 0.5), size = 4) +
  
  geom_errorbar(
    aes(xmin = lower.CL, xmax = upper.CL),
    position = position_dodge(width = 0.5),
    width = 0.2,
    orientation = "y" 
  ) +
  
  facet_wrap(~ region) +
  
  scale_color_manual(values = c("#528B8B", "#CD950C"))+
  
  scale_shape_manual(values = c("Norway" = 16, "Switzerland" = 17))+
  guides(shape = "none")+
  
  labs(title = "Effect of warming",
    x = "Δ days shifted onset (warmed − ambient)",
    y = "Phenological stage",
    color = "Biotic interactions"
  )+
  
  geom_text(aes(label = sig),
            vjust = -1, position = position_dodge(width = 0.7),
            show.legend = FALSE)
onset_all_stages_warming

# ggsave(filename = "Output/Onset/Effect_warming_bud_flower_fruit_NOR_CHE_all_in_one.png", 
#        plot = onset_all_stages_warming, width = 15, height = 12, units = "in")






