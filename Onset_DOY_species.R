


# fit the models per stage for Norway
m_onset_bud_nor_sp    <- fit_onset_model_species(onset_bud, "Norway")
m_onset_flower_nor_sp <- fit_onset_model_species(onset_flower, "Norway")
m_onset_fruit_nor_sp  <- fit_onset_model_species(onset_fruit, "Norway")
m_onset_seed_nor_sp   <- fit_onset_model_species(onset_seed, "Norway")

# check model output
# bud
summary(m_onset_bud_nor_sp)
anova(m_onset_bud_nor_sp)
model_performance(m_onset_bud_nor_sp)
#check_model(m_onset_bud_nor)

emmeans(m_onset_bud_nor_sp,
         ~ treatment_site_temp * treat_competition * species)

# fit the models per stage for Switzerland
m_onset_bud_che_sp    <- fit_onset_model_species(onset_bud, "Switzerland")
m_onset_flower_che_sp <- fit_onset_model_species(onset_flower, "Switzerland")
m_onset_fruit_che_sp  <- fit_onset_model_species(onset_fruit, "Switzerland")
m_onset_seed_che_sp   <- fit_onset_model_species(onset_seed, "Switzerland")


pred_onset_bud_che_sp <- make_onset_predictions_species_che(m_onset_bud_che_sp)
pred_onset_flower_che_sp <- make_onset_predictions_species_che(m_onset_flower_che_sp)
pred_onset_fruit_che_sp  <- make_onset_predictions_species_che(m_onset_fruit_che_sp)
pred_onset_seed_che_sp   <- make_onset_predictions_species_che(m_onset_seed_che_sp)





# combine predictions into one dataframe ----------------------------------

# nor
plot_df_bud_nor_sp  <- pred_onset_bud_nor_sp   |> 
  mutate(stage = "Budding",
         region = "Norway")

plot_df_flower_nor_sp  <- pred_onset_flower_nor_sp   |> 
  mutate(stage = "Flowering",
         region = "Norway")

plot_df_fruit_nor_sp  <- pred_onset_fruit_nor_sp   |> 
  mutate(stage = "Fruiting",
         region = "Norway")

plot_df_seed_nor_sp  <- pred_onset_seed_nor_sp   |> 
  mutate(stage = "Seeds",
         region = "Norway")

# che
plot_df_bud_che_sp  <- pred_onset_bud_che_sp   |> 
  mutate(stage = "Budding",
         region = "Switzerland")

plot_df_flower_che_sp  <- pred_onset_flower_che_sp   |> 
  mutate(stage = "Flowering",
         region = "Switzerland")

plot_df_fruit_che_sp  <- pred_onset_fruit_che_sp   |> 
  mutate(stage = "Fruiting",
         region = "Switzerland")

plot_df_seed_che_sp  <- pred_onset_seed_che_sp   |> 
  mutate(stage = "Seeds",
         region = "Switzerland")

plot_df_all_sp <- bind_rows(
  plot_df_bud_nor_sp,
  plot_df_flower_nor_sp,
  plot_df_fruit_nor_sp,
  plot_df_seed_nor_sp,
  plot_df_bud_che_sp,
  plot_df_flower_che_sp,
  plot_df_fruit_che_sp,
  plot_df_seed_che_sp
)
plot_df_all_sp


# filter by flowering -----------------------------------------------------
plot_df_all_sp_flower <- plot_df_all_sp |> 
  filter(stage == "Flowering")

plot_df_raw_all_flower <- plot_df_raw_all |> 
  filter(stage == "Flowering")


# NOR ---------------------------------------------------------------------
plot_df_all_sp_flower_nor <- plot_df_all_sp_flower |> 
  filter(region == "Norway")

plot_df_raw_all_flower_nor <- plot_df_raw_all_flower |> 
  filter(region == "Norway")




plot_df_all_sp_flower_nor <- plot_df_all_sp_flower_nor |>
  mutate(
    shape_code = case_when(
      treatment_site_temp == "lo_ambi" & treat_competition == "with"    ~ 16,
      treatment_site_temp == "lo_ambi" & treat_competition == "without" ~ 16,
      
      treatment_site_temp == "hi_ambi" & treat_competition == "with"    ~ 2,
      treatment_site_temp == "hi_ambi" & treat_competition == "without" ~ 2,
      
      treatment_site_temp == "hi_warm" & treat_competition == "with"    ~ 17,
      treatment_site_temp == "hi_warm" & treat_competition == "without" ~ 17
    )
  )
plot_df_all_sp_flower_nor


b_f_fr_sp_fl_nor <- ggplot(plot_df_all_sp_flower_nor, aes(
  x = treatment_site_temp,
  y = fit,
  color = treat_competition,
  shape = treatment_site_temp,
  group = interaction(treatment_site_temp, treat_competition)
)) +
  
  # raw:
  geom_jitter(
    data = plot_df_raw_all_flower_nor,
    aes(
      x = treatment_site_temp,
      y = onset,
      color = treat_competition,
      group =  interaction(treatment_site_temp, treat_competition)
    ),
    position = position_jitterdodge(
      jitter.width = 0.25,
      dodge.width = 0.9
    ),
    alpha = 0.4,
    size = 1.5
  ) +
  
  # model predictions
  geom_point(
    position = pd,
    size = 4,
    stroke = 1.2
  ) +
  
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.2,
    position = pd
  ) +
  
  facet_grid(~species) +
  
  scale_color_manual(values = c(
    "with" = "#528B8B",
    "without" = "#CD950C"
  )) +
  
  scale_shape_manual(values = c(
    "lo_ambi" = 16,
    "hi_ambi" = 17,
    "hi_warm" = 2
  )) +
  
  labs(
    x = "Site temperature treatment",
    y = "Onset (DOY)",
    title = "Effect of transplantation and warming on onset",
    shape = "Treatment site × warming",
    color = "Biotic interactions",
    fill = "Biotic interactions"
  ) +
  guides(shape = "none")
b_f_fr_sp_fl_nor

# ggsave(filename = "Output/Onset/Transplantation_Warming_onset_flower_species_NOR.png", 
#        plot = b_f_fr_sp_fl_nor, width = 18, height = 10, units = "in")


a <- plot_df_raw_all_flower_nor |>
  count(species, treatment_site_temp, treat_competition) |>
  arrange(n)
a

# CHE ---------------------------------------------------------------------
plot_df_all_sp_flower_che <- plot_df_all_sp_flower |> 
  filter(region == "Switzerland")

plot_df_raw_all_flower_che <- plot_df_raw_all_flower |> 
  filter(region == "Switzerland")




plot_df_all_sp_flower_che <- plot_df_all_sp_flower_che |>
  mutate(
    shape_code = case_when(
      treatment_site_temp == "lo_ambi" & treat_competition == "with"    ~ 16,
      treatment_site_temp == "lo_ambi" & treat_competition == "without" ~ 16,
      
      treatment_site_temp == "hi_ambi" & treat_competition == "with"    ~ 2,
      treatment_site_temp == "hi_ambi" & treat_competition == "without" ~ 2,
      
      treatment_site_temp == "hi_warm" & treat_competition == "with"    ~ 17,
      treatment_site_temp == "hi_warm" & treat_competition == "without" ~ 17
    )
  )
plot_df_all_sp_flower_che


b_f_fr_sp_fl_che <- ggplot(plot_df_all_sp_flower_che, aes(
  x = treatment_site_temp,
  y = fit,
  color = treat_competition,
  shape = treatment_site_temp
)) +
  
  geom_jitter(
    data = plot_df_raw_all_flower_che,
    aes(
      x = treatment_site_temp,
      y = onset,
      color = treat_competition
    ),
    position = position_jitterdodge(
      jitter.width = 0.25,
      dodge.width = 0.9
    ),
    alpha = 0.4,
    size = 1.5
  ) +
  
  # model predictions
  geom_point(
    position = pd,
    size = 4,
    stroke = 1.2
  ) +
  
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.2,
    position = pd
  ) +
  
  facet_grid(~species) +
  
  scale_color_manual(values = c(
    "with" = "#528B8B",
    "without" = "#CD950C"
  )) +
  
  scale_shape_manual(values = c(
    "lo_ambi" = 16,
    "hi_ambi" = 17,
    "hi_warm" = 2
  )) +
  
  labs(
    x = "Site temperature treatment",
    y = "Onset (DOY)",
    title = "Effect of transplantation and warming on onset",
    shape = "Treatment site × warming",
    color = "Biotic interactions",
    fill = "Biotic interactions"
  ) +
  guides(shape = "none")
b_f_fr_sp_fl_che

ggsave(filename = "Output/Onset/Transplantation_Warming_onset_flower_species_CHE.png", 
       plot = b_f_fr_sp_fl_che, width = 18, height = 10, units = "in")

b <- plot_df_raw_all_flower_che |>
  count(species, treatment_site_temp, treat_competition) |>
  arrange(n)
b




# Check plalan ------------------------------------------------------------

# NOR ---------------------------------------------------------------------
plot_df_all_sp_flower_nor_pm <- plot_df_all_sp_flower |> 
  filter(species == "plalan")

plot_df_raw_all_flower_nor_pm <- plot_df_raw_all_flower |> 
  filter(species == "plalan")


ggplot(plot_df_all_sp_flower_nor_pm, aes(
  x = treatment_site_temp,
  y = fit,
  color = treat_competition,
  shape = treatment_site_temp,
  group = interaction(treatment_site_temp, treat_competition)
)) +
  
  geom_jitter(
    data = plot_df_raw_all_flower_nor_pm,
    aes(
      x = treatment_site_temp,
      y = onset,
      color = treat_competition
    ),
    position = position_jitterdodge(
      jitter.width = 0.25,
      dodge.width = 0.9
    ),
    alpha = 0.4,
    size = 1.5
  ) +
  
  # model predictions
  geom_point(
    position = pd,
    size = 4,
    stroke = 1.2
  ) +
  
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.2,
    position = pd
  ) +
  
  facet_grid(~species) +
  
  scale_color_manual(values = c(
    "with" = "#528B8B",
    "without" = "#CD950C"
  )) +
  
  scale_shape_manual(values = c(
    "lo_ambi" = 16,
    "hi_ambi" = 17,
    "hi_warm" = 2
  )) +
  
  labs(
    x = "Site temperature treatment",
    y = "Onset (DOY)",
    title = "Effect of transplantation and warming on onset",
    shape = "Treatment site × warming",
    color = "Biotic interactions",
    fill = "Biotic interactions"
  ) +
  guides(shape = "none")

