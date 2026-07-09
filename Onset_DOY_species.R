

# Onset DOY per species ---------------------------------------------------

source("Transplantation_warming_onset_predictions_DOY_NOR_CHE.R")

#source("Functions_onset.R")



# fit the models per stage for Norway
m_onset_bud_nor_sp    <- fit_onset_model(onset_bud, "Norway")
m_onset_flower_nor_sp <- fit_onset_model(onset_flower, "Norway")
m_onset_fruit_nor_sp  <- fit_onset_model(onset_fruit, "Norway")
m_onset_seed_nor_sp   <- fit_onset_model(onset_seed, "Norway")

# check model output
# bud
summary(m_onset_bud_nor_sp)
anova(m_onset_bud_nor_sp)
model_performance(m_onset_bud_nor_sp)
#check_model(m_onset_bud_nor)

emmeans(m_onset_bud_nor_sp,
         ~ treatment_site_temp * treat_competition)

# fit the models per stage for Switzerland
m_onset_bud_che_sp    <- fit_onset_model(onset_bud, "Switzerland")
m_onset_flower_che_sp <- fit_onset_model(onset_flower, "Switzerland")
m_onset_fruit_che_sp  <- fit_onset_model(onset_fruit, "Switzerland")
m_onset_seed_che_sp   <- fit_onset_model(onset_seed, "Switzerland")

##
# m_onset_bud_nor_sp    <- fit_onset_model_species(onset_bud, "Norway")
# m_onset_flower_nor_sp <- fit_onset_model_species(onset_flower, "Norway")
# m_onset_fruit_nor_sp  <- fit_onset_model_species(onset_fruit, "Norway")
# m_onset_seed_nor_sp   <- fit_onset_model_species(onset_seed, "Norway")
# 
# m_onset_bud_che_sp    <- fit_onset_model_species(onset_bud, "Switzerland")
# m_onset_flower_che_sp <- fit_onset_model_species(onset_flower, "Switzerland")
# m_onset_fruit_che_sp  <- fit_onset_model_species(onset_fruit, "Switzerland")
# m_onset_seed_che_sp   <- fit_onset_model_species(onset_seed, "Switzerland")


# make predictions --------------------------------------------------------
pred_onset_bud_nor_sp <- make_onset_predictions_species(m_onset_bud_nor_sp)
pred_onset_flower_nor_sp <- make_onset_predictions_species(m_onset_flower_nor_sp)
pred_onset_fruit_nor_sp  <- make_onset_predictions_species(m_onset_fruit_nor_sp)
pred_onset_seed_nor_sp   <- make_onset_predictions_species(m_onset_seed_nor_sp)



pred_onset_bud_che_sp <- make_onset_predictions_species(m_onset_bud_che_sp)
pred_onset_flower_che_sp <- make_onset_predictions_species(m_onset_flower_che_sp)
pred_onset_fruit_che_sp  <- make_onset_predictions_species(m_onset_fruit_che_sp)
pred_onset_seed_che_sp   <- make_onset_predictions_species(m_onset_seed_che_sp)



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

# ggsave(filename = "Output/Onset/Transplantation_Warming_onset_flower_species_CHE.png", 
#        plot = b_f_fr_sp_fl_che, width = 18, height = 10, units = "in")

b <- plot_df_raw_all_flower_che |>
  count(species, treatment_site_temp, treat_competition) |>
  arrange(n)
b






# check -------------------------------------------------------------------

ranef(m_onset_flower_nor_sp)$species

VarCorr(m_onset_flower_nor_sp)

table(onset_flower$species)

#


#-----------------------------------------------------------------
onset_flower <- onset_flower |>
  mutate(treatment_site_temp= factor(treatment_site_temp,
                                     levels = c("lo_ambi",
                                                "hi_warm",
                                                "hi_ambi")))
species_list <- unique(onset_flower$species)
species_list

d_flower_cn <- filter_species (onset_flower, "cennig")
d_flower_cc <- filter_species (onset_flower, "cyncri")
d_flower_hm <- filter_species (onset_flower, "hypmac")
d_flower_lv <- filter_species (onset_flower, "leuvul")
d_flower_lm <- filter_species (onset_flower, "luzmul")
d_flower_ps <- filter_species (onset_flower, "pimsax")
d_flower_pl <- filter_species (onset_flower, "plalan")
d_flower_sup <- filter_species (onset_flower, "sucpra")
d_flower_tp <- filter_species (onset_flower, "tripra")
d_flower_sd <- filter_species (onset_flower, "sildio")

d_flower_bp <- filter_species (onset_flower, "brapin")
d_flower_be <- filter_species (onset_flower, "broere")
d_flower_cj <- filter_species (onset_flower, "cenjac")
d_flower_dc <- filter_species (onset_flower, "daucar")
d_flower_hp <- filter_species (onset_flower, "hypper")
d_flower_ml <- filter_species (onset_flower, "medlup")
d_flower_sv <- filter_species (onset_flower, "silvul")
d_flower_pm <- filter_species (onset_flower, "plamed")
d_flower_sp <- filter_species (onset_flower, "salpra")
d_flower_sc <- filter_species (onset_flower, "scacol")

# sensitivity models ------------------------------------------------------------------
species_list <- list(
  cennig = d_flower_cn,
  cyncri = d_flower_cc,
  hypmac = d_flower_hm,
  leuvul = d_flower_lv,
  luzmul = d_flower_lm,
  pimsax = d_flower_ps,
  plalan = d_flower_pl,
  sucpra = d_flower_sup,
  tripra = d_flower_tp,
  sildio = d_flower_sd,
  brapin = d_flower_bp,
  broere = d_flower_be,
  cenjac = d_flower_cj,
  daucar = d_flower_dc,
  hypper = d_flower_hp,
  medlup = d_flower_ml,
  silvul = d_flower_sv,
  plamed = d_flower_pm,
  #salpra = d_flower_sp, # fails the loop to model
  scacol = d_flower_sc
)
species_list



# loop through all species ------------------------------------------------
models_species_onset <- map(species_list, fit_species_model)


# get outputs -------------------------------------------------------------
model_summaries <- map(models_species_onset, summary)

model_summaries$cennig

model_summaries$leuvul



# get predictions ----------------------------------------------------------
onset_pred <- map(models_species_onset, make_onset_predictions_species)

onset_pred$cennig
onset_pred$hypmac

onset_pred2 <- map(models_species_onset, make_onset_predictions_species2)

onset_pred2$cennig
onset_pred2$hypmac




# make one dataframe ------------------------------------------------------
onset_pred_all_species <- bind_rows(onset_pred, .id = "species")
onset_pred_all_species


onset_pred_all_species <- onset_pred_all_species |>
  mutate(
    region = ifelse(
      species %in% c(
        "cennig", "cyncri", "hypmac", "leuvul", "luzmul",
        "pimsax", "plalan", "sucpra", "tripra", "sildio"),
      "Norway",
      "Switzerland"))


onset_pred_all_species_nor <- onset_pred_all_species |> 
  filter(region == "Norway")


ggplot(onset_pred_all_species_nor, aes(
  x = treatment_site_temp,
  y = fit,
  color = treat_competition,
  shape = treatment_site_temp)) +

  
  # model predictions
  geom_point(
    position = pd,
    size = 4,
    stroke = 1.2) +
  
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.2,
    position = pd) +
  
  facet_grid(region ~species) +
  
  scale_color_manual(values = c(
    "with" = "#528B8B",
    "without" = "#CD950C")) +
  
  scale_shape_manual(values = c(
    "lo_ambi" = 16,
    "hi_ambi" = 17,
    "hi_warm" = 2)) +
  
  labs(
    x = "Site temperature treatment",
    y = "Onset (DOY)",
    title = "Effect of transplantation and warming on onset NOR",
    shape = "Treatment site × warming",
    color = "Biotic interactions",
    fill = "Biotic interactions") +
  guides(shape = "none")



onset_pred_all_species_che <- onset_pred_all_species |> 
  filter(region == "Switzerland")


ggplot(onset_pred_all_species_che, aes(
  x = treatment_site_temp,
  y = fit,
  color = treat_competition,
  shape = treatment_site_temp)) +
  
  
  # model predictions
  geom_point(
    position = pd,
    size = 4,
    stroke = 1.2) +
  
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.2,
    position = pd) +
  
  facet_grid(region ~species) +
  
  scale_color_manual(values = c(
    "with" = "#528B8B",
    "without" = "#CD950C")) +
  
  scale_shape_manual(values = c(
    "lo_ambi" = 16,
    "hi_ambi" = 17,
    "hi_warm" = 2)) +
  
  labs(
    x = "Site temperature treatment",
    y = "Onset (DOY)",
    title = "Effect of transplantation and warming on onset CHE",
    shape = "Treatment site × warming",
    color = "Biotic interactions",
    fill = "Biotic interactions") +
  guides(shape = "none")




