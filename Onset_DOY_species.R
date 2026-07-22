

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
  mutate(stage = "Budding")

flower_nor <- as.data.frame(pred_onset_flower_nor_sp3) |>
  mutate(stage = "Flowering")

fruit_nor <- as.data.frame(pred_onset_fruit_nor_sp3) |>
  mutate(stage = "Fruiting")

seed_nor <- as.data.frame(pred_onset_seed_nor_sp3) |>
  mutate(stage = "Seeds")

pred_nor_all <- bind_rows(
  bud_nor,
  flower_nor,
  fruit_nor,
  seed_nor
)
pred_nor_all


bud_che <- as.data.frame(pred_onset_bud_che_sp3) |>
  mutate(stage = "Budding")

flower_che <- as.data.frame(pred_onset_flower_che_sp3) |>
  mutate(stage = "Flowering")

fruit_che <- as.data.frame(pred_onset_fruit_che_sp3) |>
  mutate(stage = "Fruiting")

seed_che <- as.data.frame(pred_onset_seed_che_sp3) |>
  mutate(stage = "Seeds")

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



species_labels_che <- c(
  brapin = "Brachypodium\npinnatum",
  broere = "Bromus\nerectus",
  cenjac = "Centaurea\njacea",
  daucar = "Daucus\ncarota",
  hypper = "Hypericum\nperforatum",
  medlup = "Medicago\nlupulina",
  plamed = "Plantago\nmedia",
  salpra = "Salvia\npratensis",
  scacol = "Scabiosa\ncolumbaria",
  silvul = "Silene\nvulgaris"
)

onset_species_che_random <- onset_species_che_random +
  facet_grid(
    stage ~ facet,
    labeller = labeller(facet = species_labels_che)) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 12),
    strip.text = element_text(size = 13, face = "italic"),
    strip.text.y = element_text(size = 13, face = "plain"))
onset_species_che_random

# ggsave(filename = "Output/Onset/Onset_DOY_species_CHE_random.png", 
#        plot = onset_species_che_random, width = 20, height = 10, units = "in")



species_labels_nor <- c(
  cennig = "Centaurea\nnigra",
  cyncri = "Cynosurus\ncristatus",
  hypmac = "Hypericum\nmaculatum",
  leuvul = "Leucanthemum\nvulgare",
  luzmul = "Luzula\nmultiflora",
  pimsax = "Pimpinella\nsaxifraga",
  plalan = "Plantago\nlanceolata",
  sildio = "Silene\ndioica",
  sucpra = "Succisa\npratensis",
  tripra = "Trifolium\npratense"
)

onset_species_nor_random <- onset_species_nor_random +
  facet_grid(
    stage ~ facet,
    labeller = labeller(facet = species_labels_nor)) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(size = 12),
    strip.text.x = element_text(size = 13, face = "italic"),
    strip.text.y = element_text(size = 13, face = "plain"))
onset_species_nor_random

# ggsave(filename = "Output/Onset/Onset_DOY_species_NOR_random.png", 
#        plot = onset_species_nor_random, width = 20, height = 10, units = "in")

