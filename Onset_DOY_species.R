

# Onset DOY per species ---------------------------------------------------

source("Transplantation_warming_onset_predictions_DOY_NOR_CHE.R")

#source("Functions_onset.R")

theme_set(theme_bw(base_size = 20))

# calculate first onset per species and plot for all stages ----------------
onset_bud    <- get_onset(phenology2, "No_Buds", "jday")
onset_flower <- get_onset(phenology2, "No_FloOpen", "jday")
onset_fruit  <- get_onset(phenology2, "No_FloWithrd", "jday")
onset_seed   <- get_onset(phenology2, "No_Seeds", "jday")



# random factor -----------------------------------------------------------

# predictions per species with species as random factor -------------------
# fit the models per stage for Norway
m_onset_bud_nor_sp3    <- fit_onset_model_species(onset_bud, "Norway")
m_onset_flower_nor_sp3 <- fit_onset_model_species(onset_flower, "Norway")
m_onset_fruit_nor_sp3  <- fit_onset_model_species(onset_fruit, "Norway")
m_onset_seed_nor_sp3   <- fit_onset_model_species(onset_seed, "Norway")

# check model output
# bud
summary(m_onset_bud_nor_sp3)
anova(m_onset_bud_nor_sp3)
model_performance(m_onset_bud_nor_sp3)
#check_model(m_onset_bud_nor)

emmeans(m_onset_bud_nor_sp3,
        ~ treatment_site_temp * treat_competition)

# fit the models per stage for Switzerland
m_onset_bud_che_sp3    <- fit_onset_model_species(onset_bud, "Switzerland")
m_onset_flower_che_sp3 <- fit_onset_model_species(onset_flower, "Switzerland")
m_onset_fruit_che_sp3  <- fit_onset_model_species(onset_fruit, "Switzerland")
m_onset_seed_che_sp3   <- fit_onset_model_species(onset_seed, "Switzerland")








# make predictions --------------------------------------------------------
pred_onset_bud_nor_sp3 <- make_onset_predictions_species3(m_onset_bud_nor_sp3)
pred_onset_flower_nor_sp3 <- make_onset_predictions_species3(m_onset_flower_nor_sp3)
pred_onset_fruit_nor_sp3  <- make_onset_predictions_species3(m_onset_fruit_nor_sp3)
pred_onset_seed_nor_sp3   <- make_onset_predictions_species3(m_onset_seed_nor_sp3)



pred_onset_bud_che_sp3 <- make_onset_predictions_species3(m_onset_bud_che_sp3)
pred_onset_flower_che_sp3 <- make_onset_predictions_species3(m_onset_flower_che_sp3)
pred_onset_fruit_che_sp3  <- make_onset_predictions_species3(m_onset_fruit_che_sp3)
pred_onset_seed_che_sp3   <- make_onset_predictions_species3(m_onset_seed_che_sp3)

bud_nor <- as.data.frame(pred_onset_bud_nor_sp3) |>
  mutate(stage = "Bud")

flower_nor <- as.data.frame(pred_onset_flower_nor_sp3) |>
  mutate(stage = "Flower")

fruit_nor <- as.data.frame(pred_onset_fruit_nor_sp3) |>
  mutate(stage = "Fruit")

seed_nor <- as.data.frame(pred_onset_seed_nor_sp3) |>
  mutate(stage = "Seed")

pred_nor_all <- bind_rows(
  bud_nor,
  flower_nor,
  fruit_nor,
  seed_nor
)
pred_nor_all


bud_che <- as.data.frame(pred_onset_bud_che_sp3) |>
  mutate(stage = "Bud")

flower_che <- as.data.frame(pred_onset_flower_che_sp3) |>
  mutate(stage = "Flower")

fruit_che <- as.data.frame(pred_onset_fruit_che_sp3) |>
  mutate(stage = "Fruit")

seed_che <- as.data.frame(pred_onset_seed_che_sp3) |>
  mutate(stage = "Seed")

pred_che_all <- bind_rows(
  bud_che,
  flower_che,
  fruit_che,
  seed_che
)
pred_che_all

names(pred_che_all)

ggplot(
  pred_che_all,
  aes(
    x = x,
    y = predicted,
    color = group,
    shape = x
  )
) +
  
  geom_point(
    position = pd,
    size = 4,
    stroke = 1.2
  ) +
  
  geom_errorbar(
    aes(
      ymin = conf.low,
      ymax = conf.high
    ),
    width = 0.2,
    position = pd
  ) +
  
  facet_grid(
    stage ~ facet) +
  
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
    color = "Biotic interactions"
  ) +
  
  guides(shape = "none") +
  
  scale_x_discrete(
    guide = guide_axis(n.dodge = 2)
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 8)
  )


theme_set(theme_bw(base_size = 20))


onset_species_che_random <- ggplot(
  flower_che,
  aes(
    x = x,
    y = predicted,
    color = group,
    shape = x
  )
) +
  
  geom_point(
    position = pd,
    size = 4,
    stroke = 1.2
  ) +
  
  geom_errorbar(
    aes(
      ymin = conf.low,
      ymax = conf.high
    ),
    width = 0.2,
    position = pd
  ) +
  
  facet_grid(
    stage ~ facet) +
  
  scale_color_manual(values = c(
    "with" = "#528B8B",
    "without" = "#CD950C"
  )) +
  
  scale_shape_manual(values = c(
    "lo_ambi" = 16,
    "hi_ambi" = 17,
    "hi_warm" = 2
  )) +
  
  labs(title = "Effect of transplantation and warming on flowering onset CHE",
    x = "Site temperature treatment",
    y = "Onset (DOY)",
    color = "Biotic interactions"
  ) +
  
  guides(shape = "none") +
  
  scale_x_discrete(
    guide = guide_axis(n.dodge = 2)
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 12)
  )
onset_species_che_random

# ggsave(filename = "Output/Onset/Onset_DOY_species_CHE_random.png", 
#        plot = onset_species_che_random, width = 20, height = 10, units = "in")


onset_species_nor_random <- ggplot(
  flower_nor,
  aes(
    x = x,
    y = predicted,
    color = group,
    shape = x
  )
) +
  
  geom_point(
    position = pd,
    size = 4,
    stroke = 1.2
  ) +
  
  geom_errorbar(
    aes(
      ymin = conf.low,
      ymax = conf.high
    ),
    width = 0.2,
    position = pd
  ) +
  
  facet_grid(
    stage ~ facet) +
  
  scale_color_manual(values = c(
    "with" = "#528B8B",
    "without" = "#CD950C"
  )) +
  
  scale_shape_manual(values = c(
    "lo_ambi" = 16,
    "hi_ambi" = 17,
    "hi_warm" = 2
  )) +
  
  labs(title = "Effect of transplantation and warming on flowering onset NOR",
    x = "Site temperature treatment",
    y = "Onset (DOY)",
    color = "Biotic interactions"
  ) +
  
  guides(shape = "none") +
  
  scale_x_discrete(
    guide = guide_axis(n.dodge = 2)
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 12)
  )
onset_species_nor_random

# ggsave(filename = "Output/Onset/Onset_DOY_species_NOR_random.png", 
#        plot = onset_species_nor_random, width = 20, height = 10, units = "in")


# one model per species ---------------------------------------------------

# predictions with one model per species ----------------------------------
# flowering
onset_flower <- onset_flower |>
  mutate(treatment_site_temp= factor(treatment_site_temp,
                                     levels = c("lo_ambi",
                                                "hi_warm",
                                                "hi_ambi")))
species_list <- unique(onset_flower$species)
species_list


# filter one dataset per species ------------------------------------------
d_flower_cn <- filter_species(onset_flower, "cennig")
d_flower_cc <- filter_species(onset_flower, "cyncri")
d_flower_hm <- filter_species(onset_flower, "hypmac")
d_flower_lv <- filter_species(onset_flower, "leuvul")
d_flower_lm <- filter_species(onset_flower, "luzmul")
d_flower_ps <- filter_species(onset_flower, "pimsax")
d_flower_pl <- filter_species(onset_flower, "plalan")
d_flower_sup <- filter_species(onset_flower, "sucpra")
d_flower_tp <- filter_species(onset_flower, "tripra")
d_flower_sd <- filter_species(onset_flower, "sildio")

d_flower_bp <- filter_species(onset_flower, "brapin")
d_flower_be <- filter_species(onset_flower, "broere")
d_flower_cj <- filter_species(onset_flower, "cenjac")
d_flower_dc <- filter_species(onset_flower, "daucar")
d_flower_hp <- filter_species(onset_flower, "hypper")
d_flower_ml <- filter_species(onset_flower, "medlup")
d_flower_sv <- filter_species(onset_flower, "silvul")
d_flower_pm <- filter_species(onset_flower, "plamed")
d_flower_sp <- filter_species(onset_flower, "salpra")
d_flower_sc <- filter_species(onset_flower, "scacol")

# make list with all species data sets ------------------------------------------------------------------
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



# fit species onset model ------------------------------------------------
# loop through all species ------------------------------------------------
models_species_onset <- map(species_list, fit_species_model)


# get outputs -------------------------------------------------------------
model_summaries <- map(models_species_onset, summary)

model_summaries$cennig

model_summaries$leuvul



# get onset predictions ----------------------------------------------------------
onset_pred <- map(models_species_onset, make_onset_predictions_species)

onset_pred$cennig
onset_pred$hypmac


# ggpredict ---------------------------------------------------------------
onset_pred2 <- map(models_species_onset, make_onset_predictions_species2)

onset_pred2$cennig
onset_pred2$hypmac




# make one data frame with all species predictions ------------------------------------------------------
onset_pred_all_species <- bind_rows(onset_pred, .id = "species")
onset_pred_all_species


# add region --------------------------------------------------------------
onset_pred_all_species <- onset_pred_all_species |>
  mutate(
    region = ifelse(
      species %in% c(
        "cennig", "cyncri", "hypmac", "leuvul", "luzmul",
        "pimsax", "plalan", "sucpra", "tripra", "sildio"),
      "Norway",
      "Switzerland"))

onset_pred_all_species <- onset_pred_all_species |>
  mutate(treatment_site_temp = factor(treatment_site_temp,
                                      levels = c("lo_ambi", "hi_warm", "hi_ambi")))

# filter by Norway --------------------------------------------------------
onset_pred_all_species_nor <- onset_pred_all_species |> 
  filter(region == "Norway")


# plot Norway ------------------------------------------------------------
onset_species_nor <- ggplot(onset_pred_all_species_nor, aes(
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
    title = "Effect of transplantation and warming on flowering onset NOR",
    shape = "Treatment site × warming",
    color = "Biotic interactions",
    fill = "Biotic interactions") +
  guides(shape = "none")+
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 15))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))
onset_species_nor

# ggsave(filename = "Output/Onset/Onset_DOY_species_NOR.png", 
#        plot = onset_species_nor, width = 20, height = 10, units = "in")


# filter by Switzerland ---------------------------------------------------
onset_pred_all_species_che <- onset_pred_all_species |> 
  filter(region == "Switzerland")


# plot all species Switzerland --------------------------------------------
onset_species_che <- ggplot(onset_pred_all_species_che, aes(
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
    title = "Effect of transplantation and warming on flowering onset CHE",
    shape = "Treatment site × warming",
    color = "Biotic interactions",
    fill = "Biotic interactions") +
  guides(shape = "none")+
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 15))+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))
onset_species_che

# ggsave(filename = "Output/Onset/Onset_DOY_species_CHE.png", 
#        plot = onset_species_che, width = 20, height = 10, units = "in")



# heatmaps ----------------------------------------------------------------

# NOR ---------------------------------------------------------------------
ggplot(onset_pred_all_species_nor,
       aes(treatment_site_temp,
           species,
           fill = fit)) +
  geom_tile() +
  scale_fill_viridis_c()


ggplot(onset_pred_all_species_nor,
       aes(treatment_site_temp,
           species,
           fill = fit)) +
  geom_tile() +
  facet_wrap(~treat_competition)



range(onset_pred_all_species_nor$fit)

onset_pred_all_species_nor |>
  group_by(treat_competition) |>
  summarise(min = min(fit),
            max = max(fit))

# CHE ---------------------------------------------------------------------
ggplot(onset_pred_all_species_che,
       aes(treatment_site_temp,
           species,
           fill = fit)) +
  geom_tile() +
  scale_fill_viridis_c()


ggplot(onset_pred_all_species_che,
       aes(treatment_site_temp,
           species,
           fill = fit)) +
  geom_tile() +
  facet_wrap(~treat_competition)

