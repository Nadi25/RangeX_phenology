


source("Temperature_sensitivity_analysis.R")





# sensitivity models including functional groups ------------------------------------------------------------------

# NOR ---------------------------------------------------------------------
# fit the models per stage for Norway
m_sens_bud_gs_nor2    <- fit_sens_model2(sens_bud_gs, "Norway")
m_sens_flower_gs_nor2 <- fit_sens_model2(sens_flower_gs, "Norway")
m_sens_fruit_gs_nor2  <- fit_sens_model2(sens_fruit_gs, "Norway")
m_sens_seed_gs_nor2   <- fit_sens_model2(sens_seed_gs, "Norway")

# check model output
summary(m_sens_bud_gs_nor2)
anova(m_sens_bud_gs_nor2)
model_performance(m_sens_bud_gs_nor2)
#check_model(m_sens_bud_gs_nor2)

summary(m_sens_flower_gs_nor2)
anova(m_sens_flower_gs_nor2)
model_performance(m_sens_flower_gs_nor2)
#check_model(m_sens_flower_gs_nor2)

summary(m_sens_fruit_gs_nor2)
anova(m_sens_fruit_gs_nor2)
model_performance(m_sens_fruit_gs_nor2)
#check_model(m_sens_fruit_gs_nor2)

summary(m_sens_seed_gs_nor2)
anova(m_sens_seed_gs_nor2)
model_performance(m_sens_seed_gs_nor2)
#check_model(m_sens_seed_gs_nor2)




# CHE ---------------------------------------------------------------------
# fit the models per stage for Switzerland
m_sens_bud_gs_che2    <- fit_sens_model2(sens_bud_gs, "Switzerland")
m_sens_flower_gs_che2 <- fit_sens_model2(sens_flower_gs, "Switzerland")
m_sens_fruit_gs_che2  <- fit_sens_model2(sens_fruit_gs, "Switzerland")
m_sens_seed_gs_che2   <- fit_sens_model2(sens_seed_gs, "Switzerland")

# check model output
summary(m_sens_bud_gs_che2)
anova(m_sens_bud_gs_che2)
model_performance(m_sens_bud_gs_che2)
#check_model(m_sens_bud_gs_che2)



summary(m_sens_flower_gs_che2)
anova(m_sens_flower_gs_che2)
model_performance(m_sens_flower_gs_che2)
#check_model(m_sens_flower_gs_che2)

summary(m_sens_fruit_gs_che2)
anova(m_sens_fruit_gs_che2)
model_performance(m_sens_fruit_gs_che2)
#check_model(m_sens_fruit_gs_che2)

summary(m_sens_seed_gs_che2)
anova(m_sens_seed_gs_che2)
model_performance(m_sens_seed_gs_che2)
#check_model(m_sens_seed_gs_che2)



# predict sensitivity for each stage with function ------------------------
# NOR ---------------------------------------------------------------------

# predict
pred_bud_gs_nor_fg    <- make_sens_predictions4(m_sens_bud_gs_nor2)
pred_flower_gs_nor_fg <- make_sens_predictions4(m_sens_flower_gs_nor2)
pred_fruit_gs_nor_fg  <- make_sens_predictions4(m_sens_fruit_gs_nor2)
pred_seed_gs_nor_fg   <- make_sens_predictions4(m_sens_seed_gs_nor2)



# CHE ---------------------------------------------------------------------
# predict
pred_bud_gs_che_fg    <- make_sens_predictions4(m_sens_bud_gs_che2)
pred_flower_gs_che_fg <- make_sens_predictions4(m_sens_flower_gs_che2)
pred_fruit_gs_che_fg  <- make_sens_predictions4(m_sens_fruit_gs_che2)
pred_seed_gs_che_fg   <- make_sens_predictions4(m_sens_seed_gs_che2)


# combine all predictions -------------------------------------------------

plot_predictions_gs_fg <- bind_rows(
  pred_bud_gs_nor_fg    |> mutate(stage = "Budding", region = "Norway"),
  pred_flower_gs_nor_fg |> mutate(stage = "Flowering", region = "Norway"),
  pred_fruit_gs_nor_fg  |> mutate(stage = "Fruiting", region = "Norway"),
  pred_seed_gs_nor_fg   |> mutate(stage = "Seeds", region = "Norway"),
  pred_bud_gs_che_fg    |> mutate(stage = "Budding", region = "Switzerland"),
  pred_flower_gs_che_fg |> mutate(stage = "Flowering", region = "Switzerland"),
  pred_fruit_gs_che_fg  |> mutate(stage = "Fruiting", region = "Switzerland"),
  pred_seed_gs_che_fg   |> mutate(stage = "Seeds", region = "Switzerland")
)
plot_predictions_gs_fg



ggplot() +
  
  geom_violin(
    data = sens_all_gs,
    aes(
      x = treat_competition,
      y = temp_sens,
      fill = treat_competition
    ),
    alpha = 0.25,
    color = NA,
    trim = FALSE
  ) +
  
  geom_point(
    data = plot_predictions_gs_fg,
    aes(
      x = treat_competition,
      y = fit,
      color = treat_competition
    ),
    size = 3,
    stroke = 1.2
  ) +
  
  geom_errorbar(
    data = plot_predictions_gs_fg,
    aes(
      x = treat_competition,
      ymin = lower,
      ymax = upper,
      color = treat_competition
    ),
    width = 0.12
  ) +
  
  facet_grid(
    region ~ stage + functional_group
  ) +
  
  scale_color_manual(values = c(
    "with" = "#528B8B",
    "without" = "#CD950C"
  )) +
  
  scale_fill_manual(values = c(
    "with" = "#528B8B",
    "without" = "#CD950C"
  )) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  theme(legend.position = "none")





# Early vs late flowering -------------------------------------------------

# https://github.com/Between-the-Fjords/seedclim-data/blob/main/databaseUtils/setup-data/moreTraits_table.csv



# Leu.vul	Prestekrage	3	FSo	NA	1	1	1	NA	SSo	20	70	Nem	1	1	1	1	0.5	NA	NA	NA	MBor (NBor)	Grassland

# Tri.pra	R<f8>dkl<f8>ver	3	FSo	NA	1	1	1	NA	SSo	15	50	Nem	1	1	1	1	1	0.5	NA	NA	NBor (LAlp)	Grassland	1	NA	NA	NA	not base-rich	Common	494	NA	NA

# Hyp.mac	Firkantperikum	1	MSo	NA	NA	1	1	NA	SSo	20	100	Nem	1	1	1	1	1	NA	NA	NA	NBor	Grassland	1	NA	NA	NA	not base-rich	Common	544	NA	NA

# Luz.mul	Engfrytle	2	FSo	NA	1	1	NA	NA	MSo	10	40	Nem	1	1	1	1	0.5	0.5	NA	NA	MBor (NBor-LAlp)	Grassland	1	NA	NA	NA	not base-rich	Common	942	NA	Amfi-atlantisk, boreal.

# Suc.pra	Bl<e5>knapp	1	SSo	NA	NA	NA	1	1	H<f8>st	20	60	Nem	1	1	1	1	1	0.5	NA	NA	NBor (LAlp)	Grassland	1	NA	NA	NA	indifferent	Common	735	NA	NA

# Pim.sax	Gjeldkarve	3	FSo	NA	1	1	1	NA	SSo	20	50	Nem	1	1	1	1	0.5	NA	NA	NA	MBor (NBor)	Grassland	1	NA	NA	NA	base-rich	Common	583	NA	NA

# Pla.lan	Smalkjempe	2	FSo	NA	1	1	NA	NA	MSo	10	50	Nem	1	1	1	1	0.5	NA	NA	NA	MBor (NBor)	Grassland	1	NA	NA	NA	not base-rich	Common	724	NA	NA

# 

# not yet
# cennig, cyncri, sildio


# get mean flowering onset per species ------------------------------------

onset_flower_general <- phenology3 |>
  filter(phenology_stage == "No_FloOpen", value > 0) |>
  group_by(region, site, year, treat_competition, species, block_ID, unique_plot_ID, 
           unique_plant_ID, phenology_stage, functional_group) |>
  summarise(onset = min(jday), .groups = "drop") |>
  # remove groups where budding never occurred
  filter(is.finite(onset))
onset_flower_general


mean_onset_flower_species_region <- onset_flower_general |>
  group_by(region, species) |>
  summarise(mean_onset = mean(onset, na.rm = TRUE),
            .groups = "drop")
mean_onset_flower_species_region



# define early and late flowering species ---------------------------------
# based on mean onset per species and region
species_flowering <- mean_onset_flower_species_region |>
  group_by(region) |>
  mutate(cutoff = median(mean_onset),
    flowering_time = if_else(mean_onset <= cutoff, "early", "late"))
species_flowering


table(species_flowering$region, species_flowering$flowering_time)



# join flowering time in on rest of phenology data ------------------------
phenology4 <- phenology3 |> 
  left_join(species_flowering, by = c("region", "species")) |> 
  select(-c(cutoff, mean_onset))




get_mean_onset2 <- function(data, stage_name, onset_type) {
  
  data |>
    filter(
      phenology_stage == stage_name,
      value > 0
    ) |>
    
    # first onset per individual
    group_by(
      region, site, treat_competition,
      species, block_ID, unique_plot_ID,
      unique_plant_ID, functional_group, flowering_time
    ) |>
    summarise(
      onset = min(.data[[onset_type]]),
      .groups = "drop"
    ) #|>
  
  # # mean onset across 3 individuals within plot
  # group_by(
  #   region, site, treat_competition,
  #   species, block_ID, unique_plot_ID
  # ) |>
  # summarise(
  #   onset = mean(first_onset),
  #   .groups = "drop"
  # )
}

# calculate mean onset per species and plot for all stages ----------------
onset_bud2    <- get_mean_onset2(phenology4, "No_Buds", "jday")
onset_flower2 <- get_mean_onset2(phenology4, "No_FloOpen", "jday")
onset_fruit2  <- get_mean_onset2(phenology4, "No_FloWithrd", "jday")
onset_seed2   <- get_mean_onset2(phenology4, "No_Seeds", "jday")


# Average across species, site and treat -----------------------------------------
# this is not needed or at least it doesnt change from get_mean_onset
# because unique_plot_ID is essentially treat_comp
# but we still need an average across the plots
get_species_onset2 <- function(onset_data) {
  onset_data |>
    group_by(region, site, treat_competition, species, block_ID, functional_group, flowering_time) |>
    summarise(
      onset = mean(onset, na.rm = TRUE),
      .groups = "drop"
    )
}

# Average across species, site and treat -----------------------------------------
onset_bud_mean2    <- get_species_onset2(onset_bud2)
onset_flower_mean2 <- get_species_onset2(onset_flower2)
onset_fruit_mean2  <- get_species_onset2(onset_fruit2)
onset_seed_mean2   <- get_species_onset2(onset_seed2)
# this doesn't change the values


# Calculate temperature sensitivity ---------------------------------------
# function to get temp sens
get_temp_sens2 <- function(onset_mean_data, temperature_data) {
  onset_mean_data |>
    left_join(temperature_data, by = c("region","site")) |>
    group_by(region, species, treat_competition, block_ID, functional_group, flowering_time) |>  # or group by site and species?
    pivot_wider(names_from = site,
                values_from = c(onset, Tmean)) |>
    mutate(
      temp_sens = (onset_lo - onset_hi) / (Tmean_lo - Tmean_hi)
    )
}



# Calculate temperature sensitivity ---------------------------------------
# growing season mean ----------------------------------
# get temp sens per stage while using the same delta T 
sens_bud_gs2   <- get_temp_sens(onset_bud_mean2, temperature_mean_gs)
sens_flower_gs2    <- get_temp_sens(onset_flower_mean2, temperature_mean_gs)
sens_fruit_gs2    <- get_temp_sens(onset_fruit_mean2, temperature_mean_gs)
sens_seed_gs2   <- get_temp_sens(onset_seed_mean2, temperature_mean_gs)



# sensitivity models ------------------------------------------------------------------
fit_sens_model3 <- function(sens_data, region_name) {
  sens_region <- sens_data |> 
    filter(region == region_name)
  
  model <- lmerTest::lmer(
    temp_sens ~ treat_competition * flowering_time + (1 | species) + (1 | block_ID),
    data = sens_region
  )
  
  return(model)
}


# NOR ---------------------------------------------------------------------
# fit the models per stage for Norway
m_sens_bud_gs_nor2    <- fit_sens_model3(sens_bud_gs2, "Norway")
m_sens_flower_gs_nor2 <- fit_sens_model3(sens_flower_gs2, "Norway")
m_sens_fruit_gs_nor2  <- fit_sens_model3(sens_fruit_gs2, "Norway")
m_sens_seed_gs_nor2   <- fit_sens_model3(sens_seed_gs2, "Norway")


summary(m_sens_bud_gs_nor2)
summary(m_sens_flower_gs_nor2)
summary(m_sens_fruit_gs_nor2)
summary(m_sens_seed_gs_nor2)



# CHE ---------------------------------------------------------------------
# fit the models per stage for Switzerland
m_sens_bud_gs_che2    <- fit_sens_model3(sens_bud_gs2, "Switzerland")
m_sens_flower_gs_che2 <- fit_sens_model3(sens_flower_gs2, "Switzerland")
m_sens_fruit_gs_che2  <- fit_sens_model3(sens_fruit_gs2, "Switzerland")
m_sens_seed_gs_che2   <- fit_sens_model3(sens_seed_gs2, "Switzerland")

# check model output
summary(m_sens_bud_gs_che2)
summary(m_sens_flower_gs_che2)
summary(m_sens_fruit_gs_che2)
summary(m_sens_seed_gs_che2)




# predict sensitivity for each stage with function ------------------------
make_sens_predictions5 <- function(model) {
  
  newdata <- expand.grid(treat_competition = c("with", "without"),
                         flowering_time = c("early", "late")) |> 
    as_tibble()
  
  pred <- predict(
    model,
    newdata = newdata,
    re.form = NA,
    se.fit = TRUE) |> 
    
    as_tibble() |>
    mutate(upper = fit + 1.96 * se.fit,
           lower = fit - 1.96 * se.fit) |>
    bind_cols(newdata)
  
  return(pred)
}


# NOR ---------------------------------------------------------------------
# predict
pred_bud_gs_nor2    <- make_sens_predictions5(m_sens_bud_gs_nor2)
pred_flower_gs_nor2 <- make_sens_predictions5(m_sens_flower_gs_nor2)
pred_fruit_gs_nor2  <- make_sens_predictions5(m_sens_fruit_gs_nor2)
pred_seed_gs_nor2   <- make_sens_predictions5(m_sens_seed_gs_nor2)



# CHE ---------------------------------------------------------------------
# predict
pred_bud_gs_che2    <- make_sens_predictions5(m_sens_bud_gs_che2)
pred_flower_gs_che2 <- make_sens_predictions5(m_sens_flower_gs_che2)
pred_fruit_gs_che2  <- make_sens_predictions5(m_sens_fruit_gs_che2)
pred_seed_gs_che2   <- make_sens_predictions5(m_sens_seed_gs_che2)


# combine all predictions -------------------------------------------------
plot_predictions_gs2 <- bind_rows(
  pred_bud_gs_nor2    |> mutate(stage = "Budding", region = "Norway"),
  pred_flower_gs_nor2 |> mutate(stage = "Flowering", region = "Norway"),
  pred_fruit_gs_nor2  |> mutate(stage = "Fruiting", region = "Norway"),
  pred_seed_gs_nor2   |> mutate(stage = "Seeds", region = "Norway"),
  pred_bud_gs_che2    |> mutate(stage = "Budding", region = "Switzerland"),
  pred_flower_gs_che2 |> mutate(stage = "Flowering", region = "Switzerland"),
  pred_fruit_gs_che2  |> mutate(stage = "Fruiting", region = "Switzerland"),
  pred_seed_gs_che2   |> mutate(stage = "Seeds", region = "Switzerland")
)
plot_predictions_gs2



ggplot() +
  
  # geom_violin(
  #   data = sens_all_gs,
  #   aes(
  #     x = treat_competition,
  #     y = temp_sens,
  #     fill = treat_competition
  #   ),
  #   alpha = 0.25,
  #   color = NA,
  #   trim = FALSE
  # ) +
  
  geom_point(
    data = plot_predictions_gs2,
    aes(
      x = treat_competition,
      y = fit,
      color = treat_competition
    ),
    size = 3,
    stroke = 1.2
  ) +
  
  geom_errorbar(
    data = plot_predictions_gs2,
    aes(
      x = treat_competition,
      ymin = lower,
      ymax = upper,
      color = treat_competition
    ),
    width = 0.12
  ) +
  
  facet_grid(
    region ~ stage + flowering_time
  ) +
  
  scale_color_manual(values = c(
    "with" = "#528B8B",
    "without" = "#CD950C"
  )) +
  
  scale_fill_manual(values = c(
    "with" = "#528B8B",
    "without" = "#CD950C"
  )) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  theme(legend.position = "none")




