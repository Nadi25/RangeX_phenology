
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

  annotate("rect",
           xmin = -Inf, xmax = 0,
           ymin = -Inf, ymax = Inf,
           fill = "salmon1", alpha = 0.2)+
  
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
            show.legend = FALSE)+
  scale_x_continuous(
    limits = c(-18, 30))

onset_all_stages_cooling


# ggsave(filename = "Output/Onset/Effect_cooling_bud_flower_fruit_NOR_CHE_all_in_one.png", 
#        plot = onset_all_stages_cooling, width = 15, height = 10, units = "in")



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
  
  annotate("rect",
           xmin = -Inf, xmax = 0,
           ymin = -Inf, ymax = Inf,
           fill = "salmon1", alpha = 0.2)+
  
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
            show.legend = FALSE)+
  scale_x_continuous(
    limits = c(-18, 30))
onset_all_stages_warming

# ggsave(filename = "Output/Onset/Effect_warming_bud_flower_fruit_NOR_CHE_all_in_one.png", 
#        plot = onset_all_stages_warming, width = 15, height = 10, units = "in")




# Warming interaction -----------------------------------------------------

# Combined plot with effect of warming, biotic interaction --------
# source all the effect warming scripts

source("Effect_warming_on_budding_onset_NOR_CHE.R")

source("Effect_warming_on_flower_onset_NOR_CHE.R")

source("Effect_warming_on_fruiting_onset_NOR_CHE.R")


stage_colors <- c(
  "Budding"   = "#4F9EC9",   
  "Flowering" = "pink3",   
  "Fruiting"  =  "#F4A636"   
)


plot_df_warm_bud  <- plot_df_warm_bud   |> mutate(stage = "Budding")

plot_df_warm_flower <- plot_df_warm_flower |> mutate(stage = "Flowering")

plot_df_warm_fruit <- plot_df_warm_fruit |> mutate(stage = "Fruiting")


plot_df_all_warm <- bind_rows(
  plot_df_warm_bud,
  plot_df_warm_flower,
  plot_df_warm_fruit
) |>
  select(region, effect, group, estimate, lower.CL, upper.CL, stars, stage)
plot_df_all_warm



plot_df_all_warm <- plot_df_all_warm |>
  mutate(effect = factor(effect,
                         levels = c("Warming",
                                    "Biotic interactions",
                                    "Interaction")))

plot_df_all_warm <- plot_df_all_warm |>
  mutate(
    group = case_when(
      group == "with" ~ "with bi",
      group == "without" ~ "without bi",
      group == "warm" ~ "warmed",
      group == "ambi" ~ "ambient",
      group == "all interactions" ~ "warming × bi",
      TRUE ~ group
    )
  ) |>
  mutate(
    group = stringr::str_wrap(group, width = 12)
  ) |>
  mutate(
    group = factor(group,
                   levels = stringr::str_wrap(
                     c("warming × bi",
                       "with bi",
                       "without bi",
                       "warmed",
                       "ambient"),
                     width = 12
                   ))
  )
plot_df_all_warm

warming_all <- ggplot(plot_df_all_warm,
                      aes(x = group, y = estimate, color = stage, shape  = region)) +
  
  geom_point(size = 6,
             position = position_dodge(width = 0.6)) +
  
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0.1, linewidth = 1, 
                position = position_dodge(width = 0.6)) +
  
  geom_text(aes(y = upper.CL, label = stars),
            position = position_dodge(width = 0.6),
            vjust = -0.3,
            size = 7,
            show.legend = FALSE) +
  
  facet_grid(effect ~ region, scales = "free_y") +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  labs(x = NULL,
       y = "Shift in onset (days)",
       color = "Phenological stage")+
  
  scale_color_manual(values = stage_colors)+
  
  theme(
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20)
  ) +
  
  scale_shape_manual(values = c("Norway" = 16, "Switzerland" = 17))+
  guides(shape = "none")

warming_all


# flip axis looks better
warming_all <- warming_all + coord_flip()
warming_all

# ggsave(filename = "Output/Onset/Onset_warming_effect_NOR_CHE_all_interactions.png", 
#        plot = warming_all, width = 17, height = 16, units = "in")



# GDD NOR ---------------------------------------------------------------------

source("GDD_Cooling_budding_onset_NOR.R")

source("GDD_Cooling_flowering_onset_NOR.R")

source("GDD_Cooling_fruiting_onset_NOR.R")



# One figure all stages -------------------------------------------
contrast_df_nor_gdd_ambi_bud$stage <- "Budding"
contrast_df_nor_gdd_ambi$stage <- "Flowering"
contrast_df_nor_gdd_ambi_fruit$stage <- "Fruiting"

contr_all_gdd_cool <- bind_rows(contrast_df_nor_gdd_ambi_bud, contrast_df_nor_gdd_ambi, contrast_df_nor_gdd_ambi_fruit)

contr_all_gdd_cool$stage <- factor(contr_all_gdd_cool$stage,
                               levels = c("Fruiting", "Flowering", "Budding"))

contr_all_gdd_cool <- contr_all_gdd_cool |>
  mutate(sig = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE ~ ""
  ))
contr_all_gdd_cool


onset_all_stages_gdd_cool <- ggplot(contr_all_gdd_cool,
                                   aes(x = estimate, y = stage, color = treat_competition)) +
  
  annotate("rect",
           xmin = -Inf, xmax = 0,
           ymin = -Inf, ymax = Inf,
           fill = "salmon1", alpha = 0.2)+
  
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  geom_point(position = position_dodge(width = 0.5), size = 4) +
  
  geom_errorbar(
    aes(xmin = lower.CL, xmax = upper.CL),
    position = position_dodge(width = 0.5),
    width = 0.2,
    orientation = "y" 
  ) +
  
  scale_color_manual(values = c("#528B8B", "#CD950C"))+
  
  labs(title = "Effect of transplantation using GDD",
       x = "Δ days shifted onset (high − low)",
       y = "Phenological stage",
       color = "Biotic interactions"
  )+
  
  geom_text(aes(label = sig),
            vjust = -1, position = position_dodge(width = 0.7),
            show.legend = FALSE)
onset_all_stages_gdd_cool

# ggsave(filename = "Output/Onset/GDD_effect_transplantation_bud_flower_fruit_NOR_CHE_all_in_one.png", 
#        plot = onset_all_stages_gdd_cool, width = 15, height = 10, units = "in")






# Interactions ------------------------------------------------------------

# Combined plot with effect of transplantation, biotic interaction --------
# source all the effect cooling scripts

source("Effect_cooling_on_budding_onset_NOR_CHE.R")

source("Effect_cooling_on_flowering_onset_NOR_CHE.R")

source("Effect_cooling_on_fruiting_onset_NOR_CHE.R")


stage_colors <- c(
  "Budding"   = "#7FBF4F",   
  "Flowering" = "#4F9EC9",   
  "Fruiting"  =  "#F4A636"   
)

stage_colors <- c(
  "Budding"   = "#4F9EC9",   
  "Flowering" = "pink3",   
  "Fruiting"  =  "#F4A636"   
)


df_cooling_bud   <- df_cooling_bud   |> mutate(stage = "Budding")
df_comp_bud      <- df_comp_bud      |> mutate(stage = "Budding")
df_interaction_bud <- df_interaction_bud |> mutate(stage = "Budding")

df_cooling_flower <- df_cooling |> mutate(stage = "Flowering")
df_comp_flower    <- df_comp    |> mutate(stage = "Flowering")
df_interaction_flower <- df_interaction |> mutate(stage = "Flowering")

df_cooling_fruit <- df_cooling_fruit |> mutate(stage = "Fruiting")
df_comp_fruit    <- df_comp_fruit    |> mutate(stage = "Fruiting")
df_interaction_fruit <- df_interaction_fruit |> mutate(stage = "Fruiting")


plot_df_all <- bind_rows(
  df_cooling_bud, df_comp_bud, df_interaction_bud,
  df_cooling_flower, df_comp_flower, df_interaction_flower,
  df_cooling_fruit, df_comp_fruit, df_interaction_fruit
) |>
  select(region, effect, group, estimate, lower.CL, upper.CL, stars, stage)
plot_df_all


plot_df_all <- plot_df_all |>
  mutate(
    x_group = case_when(
      effect == "Transplantation" ~ group,        # with / without
      effect == "Biotic interactions" ~ group,    # hi / lo
      effect == "Interaction" ~ "all interaction"
    )
  )
plot_df_all

plot_df_all <- plot_df_all |>
  mutate(effect = factor(effect,
                         levels = c("Transplantation",
                                    "Biotic interactions",
                                    "Interaction")))



plot_df_all <- plot_df_all |>
  mutate(
    group = case_when(
      group == "with" ~ "with bi",
      group == "without" ~ "without bi",
      group == "hi" ~ "high site",
      group == "lo" ~ "low site",
      group == "all interactions" ~ "site × bi",
      TRUE ~ group
    )
  ) |>
  mutate(
    group = stringr::str_wrap(group, width = 12)
  ) |>
  mutate(
    group = factor(group,
                   levels = stringr::str_wrap(
                     c("site × bi",
                       "with bi",
                       "without bi",
                       "low site",
                       "high site"),
                     width = 12
                   ))
  )




cooling_all <- ggplot(plot_df_all,
                      aes(x = group, y = estimate, color = stage, shape  = region)) +
  
  geom_point(size = 6,
             position = position_dodge(width = 0.6)) +
  
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0.1, linewidth = 1, 
                position = position_dodge(width = 0.6)) +
  
  geom_text(aes(y = upper.CL, label = stars),
            position = position_dodge(width = 0.6),
            vjust = -0.3,
            size = 7,
            show.legend = FALSE) +
  
  facet_grid(effect ~ region, scales = "free_y") +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  labs(x = NULL,
       y = "Shift in onset (days)",
       color = "Phenological stage")+
  
  scale_color_manual(values = stage_colors)+
  
  theme(
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20)
  ) +
  
  scale_shape_manual(values = c("Norway" = 16, "Switzerland" = 17))+
  guides(shape = "none")

cooling_all

# ggsave(filename = "Output/Onset/Onset_cooling_effect_NOR_CHE_all_interactions.png", 
#        plot = cooling_all, width = 17, height = 16, units = "in")

cooling_all2 <- cooling_all + coord_flip()
cooling_all2

# ggsave(filename = "Output/Onset/Onset_cooling_effect_NOR_CHE_all_interactions2.png", 
#        plot = cooling_all2, width = 17, height = 16, units = "in")



# GDD ---------------------------------------------------------------------

# Combined plot with effect of transplantation, biotic interaction using GDD --------
# source all the effect cooling gdd scripts

source("GDD_Cooling_budding_onset_NOR_CHE.R")

source("GDD_Cooling_flowering_onset_NOR_CHE.R")

source("GDD_Cooling_fruiting_onset_NOR_CHE.R")


stage_colors <- c(
  "Budding"   = "#4F9EC9",   
  "Flowering" = "pink3",   
  "Fruiting"  =  "#F4A636"   
)


plot_df_gdd_bud  <- plot_df_gdd_bud   |> mutate(stage = "Budding")

plot_df_gdd_flower <- plot_df_gdd_flower |> mutate(stage = "Flowering")

plot_df_gdd_fruit <- plot_df_gdd_fruit |> mutate(stage = "Fruiting")


plot_df_all <- bind_rows(
  plot_df_gdd_bud,
  plot_df_gdd_flower,
  plot_df_gdd_fruit
) |>
  select(region, effect, group, estimate, lower.CL, upper.CL, stars, stage)
plot_df_all


# plot_df_all <- plot_df_all |>
#   mutate(
#     effect = case_when(
#       effect == "Transplantation" ~ group,        # with / without
#       effect == "Biotic interactions" ~ group,    # hi / lo
#       effect == "Interaction" ~ "all interactions"
#     )
#   )
# plot_df_all

plot_df_all <- plot_df_all |>
  mutate(effect = factor(effect,
                         levels = c("Transplantation",
                                    "Biotic interactions",
                                    "Interaction")))

plot_df_all <- plot_df_all |>
  mutate(
    group = case_when(
      group == "with" ~ "with bi",
      group == "without" ~ "without bi",
      group == "hi" ~ "high site",
      group == "lo" ~ "low site",
      group == "all interactions" ~ "site × bi",
      TRUE ~ group
    )
  ) |>
  mutate(
    group = stringr::str_wrap(group, width = 12)
  ) |>
  mutate(
    group = factor(group,
                   levels = stringr::str_wrap(
                     c("site × bi",
                       "with bi",
                       "without bi",
                       "low site",
                       "high site"),
                     width = 12
                   ))
  )
plot_df_all

cooling_all_gdd <- ggplot(plot_df_all,
                          aes(x = group, y = estimate, color = stage, shape  = region)) +
  
  geom_point(size = 6,
             position = position_dodge(width = 0.6)) +
  
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0.1, linewidth = 1, 
                position = position_dodge(width = 0.6)) +
  
  geom_text(aes(y = upper.CL, label = stars),
            position = position_dodge(width = 0.6),
            vjust = -0.3,
            size = 7,
            show.legend = FALSE) +
  
  facet_grid(effect ~ region, scales = "free_y") +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  labs(x = NULL,
       y = "Shift in onset (GDD)",
       color = "Phenological stage")+
  
  scale_color_manual(values = stage_colors)+
  
  theme(
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20)
  ) +
  
  scale_shape_manual(values = c("Norway" = 16, "Switzerland" = 17))+
  guides(shape = "none")

cooling_all_gdd



# ggsave(filename = "Output/Onset/GDD_onset_cooling_effect_NOR_CHE_all_interactions.png", 
#        plot = cooling_all_gdd, width = 17, height = 16, units = "in")


# flip axis looks better
cooling_all_gdd2 <- cooling_all_gdd + coord_flip()
cooling_all_gdd2

# ggsave(filename = "Output/Onset/GDD_onset_cooling_effect_NOR_CHE_all_interactions2.png", 
#        plot = cooling_all_gdd2, width = 17, height = 16, units = "in")



