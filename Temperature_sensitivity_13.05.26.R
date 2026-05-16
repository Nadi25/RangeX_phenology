


# load library ------------------------------------------------------------
library(conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
library(tidyverse)
library(lubridate)
library(lme4)
library(performance)
library(see)
library(emmeans)



# calcualte delta T -------------------------------------------------------

# mean temperature per site
# timeframe before the first phenological event


# use climate station data for NOR
# and a combination of tomst and climate station for CHE


# source climate scripts --------------------------------------------------
source("Data_preparation_climate_station_NOR.R")
climate_23 
climate_23_daily

source("Data_preparation_climate_station_CHE.R")
climate_che_combined
# which is daily per site



# filter time period for pre-season window --------------------------------

# 1. 30 days before first mean budding onset
# 1 Norway            172. # 21.06
# 2 Switzerland       156. # 05.06

# preclimate 30 - 30 days ------------------------------------------------------------

# NOR
climate_nor_23_pre30 <- climate_23_daily |>
  filter(
    date >= as.Date("2023-05-21"),
    date <= as.Date("2023-06-20")
  ) |> 
  mutate(region = "Norway")



# CHE
climate_che_22_pre30 <- climate_che_combined |>
  filter(
    date >= as.Date("2022-05-05"),
    date <= as.Date("2022-06-04")
  ) |> 
  mutate(region = "Switzerland")


# combine
climate_all30 <- bind_rows(climate_che_22_pre30, climate_nor_23_pre30)


# calculate mean per site in this time period
pre_climate30 <-  climate_all30 |>
  group_by(region, site) |>
  summarise(Tmean = mean(Tavg))
pre_climate30


# preclimate 14 - 14 days ------------------------------------------------------------

# NOR
climate_nor_23_pre14 <- climate_23_daily |>
  filter(
    date >= as.Date("2023-06-07"),
    date <= as.Date("2023-06-20")
  ) |> 
  mutate(region = "Norway")


# CHE
climate_che_22_pre14 <- climate_che_combined |>
  filter(
    date >= as.Date("2022-05-22"),
    date <= as.Date("2022-06-04")
  ) |> 
  mutate(region = "Switzerland")



# combine
climate_all14 <- bind_rows(climate_che_22_pre14, climate_nor_23_pre14)


# calculate mean per site in this time period
pre_climate14 <-  climate_all14 |>
  group_by(region, site) |>
  summarise(Tmean = mean(Tavg))
pre_climate14



# preclimate 60 - 60 days ------------------------------------------------------------

# NOR
climate_nor_23_pre60 <- climate_23_daily |>
  filter(
    date >= as.Date("2023-04-21"),
    date <= as.Date("2023-06-20")
  ) |> 
  mutate(region = "Norway")


# CHE
climate_che_22_pre60 <- climate_che_combined |>
  filter(
    date >= as.Date("2022-04-05"),
    date <= as.Date("2022-06-04")
  ) |> 
  mutate(region = "Switzerland")


# combine
climate_all60 <- bind_rows(climate_che_22_pre60, climate_nor_23_pre60)


# calculate mean per site in this time period
pre_climate60 <-  climate_all60 |>
  group_by(region, site) |>
  summarise(Tmean = mean(Tavg))
pre_climate60





# preclimate growing season mean ------------------------------------------

# NOR
climate_nor_23_pre_gs <- climate_23_daily |>
  filter(
    date >= as.Date("2023-04-01"),
    date <= as.Date("2023-10-31")
  ) |> 
  mutate(region = "Norway")


# CHE
climate_che_22_pre_gs <- climate_che_combined |>
  filter(
    date >= as.Date("2022-04-01"),
    date <= as.Date("2022-10-31")
  ) |> 
  mutate(region = "Switzerland")


# combine
climate_all_gs <- bind_rows(climate_che_22_pre_gs, climate_nor_23_pre_gs)


# calculate mean per site in this time period
pre_climate_gs <-  climate_all_gs |>
  group_by(region, site) |>
  summarise(Tmean = mean(Tavg))
pre_climate_gs



# load clean phenology data -----------------------------------------------
source("Data_preparation_phenology_NOR_CHE_combined.R")

# use this data set
names(phenology)

# set theme for plots for presentation ------------------------------------
theme_set(theme_bw(base_size = 20))


# rename infructescence stage of NOR to FlowWithrd ------------------------
# this will be the fruiting stage
# combine the fruiting stages nor and che to compare the onset
# this is to be taken with caution because the stages are not the same
# but for the onset it could be comparable
phenology <- phenology |>
  mutate(phenology_stage = recode(phenology_stage,
                                  "No_Infructescences" = "No_FloWithrd"))

# and get julian days ---------------------------------------------------
# yday(date)
# che and nor was measured in two years but if we count the days in each year it should be fine

phenology3 <- phenology |> 
  mutate(
    jday = yday(date_measurement),   # Julian day (1–365)
    jday_scaled = scale(jday))        # optional 



# filter only ambi both sites  -------------------------------------
# to compare low ambi with hi ambi = cooling effect
phenology3 <- phenology3 |> 
  filter(treat_warming == "ambi")



phenology3 <- phenology3 |> 
  mutate(year = if_else(region == "Switzerland", 2022, 2023))


# calculate onsets ---------------------------------------------------------


# bud find mean for pre-climate time---------------------------------------------------------------------
onset_bud <- phenology3 |>
  filter(phenology_stage == "No_Buds", value > 0) |>
  group_by(region, site, year, treat_competition, species, block_ID, unique_plot_ID, unique_plant_ID, phenology_stage) |>
  summarise(onset = min(jday), .groups = "drop") |>
  # remove groups where budding never occurred
  filter(is.finite(onset))

first_onset_bud <- onset_bud |>
  group_by(region, site, year) |>
  summarise(first_onset = min(onset, na.rm = TRUE),
            .groups = "drop")
first_onset_bud

mean_onset_bud_species <- onset_bud |>
  group_by(region, site, year, species) |>
  summarise(mean_onset = mean(onset, na.rm = TRUE),
            .groups = "drop")
mean_onset_bud_species

mean_onset_bud <- onset_bud |>
  group_by(region, site, year) |>
  summarise(mean_onset = mean(onset, na.rm = TRUE),
            .groups = "drop")
mean_onset_bud

# 1 Norway      hi     2023       181. # 30.06
# 2 Norway      lo     2023       164. # 13.06
# 3 Switzerland hi     2022       165. # 14.06
# 4 Switzerland lo     2022       149. # 29.05

mean_onset_bud_site <- onset_bud |>
  group_by(region) |>
  summarise(mean_onset = mean(onset, na.rm = TRUE),
            .groups = "drop")
mean_onset_bud_site

# 1 Norway            172. # 21.06
# 2 Switzerland       156. # 05.06




# function to calculate mean onset per plot ------------------------------
# take the average onset for the three individuals per species in one plot
# first budding date per individual
# then average per plot
get_mean_onset <- function(data, stage_name) {
  
  data |>
    filter(
      phenology_stage == stage_name,
      value > 0
    ) |>
    
    # first onset per individual
    group_by(
      region, site, treat_competition,
      species, block_ID, unique_plot_ID,
      unique_plant_ID
    ) |>
    summarise(
      first_onset = min(jday),
      .groups = "drop"
    ) |>
    
    # mean onset across 3 individuals within plot
    group_by(
      region, site, treat_competition,
      species, block_ID, unique_plot_ID
    ) |>
    summarise(
      onset = mean(first_onset),
      .groups = "drop"
    )
}


# calculate mean onset per species and plot for all stages ----------------
onset_bud    <- get_mean_onset(phenology3, "No_Buds")
onset_flower <- get_mean_onset(phenology3, "No_FloOpen")
onset_fruit  <- get_mean_onset(phenology3, "No_FloWithrd")
onset_seed   <- get_mean_onset(phenology3, "No_Seeds")

# combine the fruiting stages nor and che to compare the onset
# this is to be taken with caution because the stages are not the same
# but for the onset it could be comparable




# Average across species, site and treat -----------------------------------------

get_species_onset <- function(onset_data) {
  onset_data |>
    group_by(region, site, treat_competition, species) |>
    summarise(
      onset = mean(onset, na.rm = TRUE),
      .groups = "drop"
    )
}

onset_bud_mean    <- get_species_onset(onset_bud)
onset_flower_mean <- get_species_onset(onset_flower)
onset_fruit_mean  <- get_species_onset(onset_fruit)
onset_seed_mean   <- get_species_onset(onset_seed)



# Calculate temperature sensitivity ---------------------------------------
# function to get temp sens
get_temp_sens <- function(onset_mean_data, pre_climate_data) {
  onset_mean_data |>
    left_join(pre_climate_data, by = c("region","site")) |>
    group_by(region, species, treat_competition ) |>  # or group by site and species?
    pivot_wider(names_from = site,
                values_from = c(onset, Tmean)) |>
    mutate(
      temp_sens = (onset_lo - onset_hi) / (Tmean_lo - Tmean_hi)
    )
}



# temperature sensitivity 30 days window ----------------------------------
# get temp sens per stage for pre_climate1 = 30 days
sens_bud_30    <- get_temp_sens(onset_bud_mean, pre_climate30)
sens_flower_30    <- get_temp_sens(onset_flower_mean, pre_climate30)
sens_fruit_30    <- get_temp_sens(onset_fruit_mean, pre_climate30)
sens_seed_30   <- get_temp_sens(onset_seed_mean, pre_climate30)



# control plots sensitivity -----------------------------------------------
ggplot(sens_bud_30, aes(x = treat_competition, y = temp_sens, color = species)) +
  geom_point(position = position_jitter(width = 0.15, height = 0), size = 2, alpha = 0.85) +
  facet_wrap(~region)


ggplot(sens_flower_30, aes(x = treat_competition, y = temp_sens, color = species)) +
  geom_point(position = position_jitter(width = 0.15, height = 0), size = 2, alpha = 0.85) +
  facet_wrap(~region)


ggplot(sens_fruit_30, aes(x = treat_competition, y = temp_sens, color = species)) +
  geom_point(position = position_jitter(width = 0.15, height = 0), size = 2, alpha = 0.85) +
  facet_wrap(~region)


ggplot(sens_seed_30, aes(x = treat_competition, y = temp_sens, color = species)) +
  geom_point(position = position_jitter(width = 0.15, height = 0), size = 2, alpha = 0.85) +
  facet_wrap(~region)




# combine sens from all stages -------------------------------------------
sens_all_30 <- bind_rows(
  sens_bud_30   |> mutate(stage = "Budding"),
  sens_flower_30|> mutate(stage = "Flowering"),
  sens_fruit_30 |> mutate(stage = "Fruiting"),
  sens_seed_30  |> mutate(stage = "Seeds")
)
sens_all_30



# quick control plot ------------------------------------------------------
ggplot(sens_all_30,
       aes(x = stage,
           y = temp_sens,
           color = treat_competition)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.4),
             alpha = 0.7) 


# plot the raw sens data 
ggplot(sens_all_30,
       aes(x = stage,
           y = temp_sens,
           color = treat_competition)) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  # individual species
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.1,
      dodge.width = 0.4
    ),
    alpha = 0.3
  ) +
  
  # mean ± 95% CI
  stat_summary(
    fun.data = mean_cl_normal,
    geom = "errorbar",
    position = position_dodge(width = 0.4),
    width = 0.15,
    linewidth = 0.8
  ) +
  
  stat_summary(
    fun = mean,
    geom = "point",
    position = position_dodge(width = 0.4),
    size = 4
  ) +
  
  labs(
    x = "Phenological stage",
    y = expression("Temperature sensitivity (days / °C)"),
    color = "Biotic interactions"
  )+
  facet_wrap(~ region)



# sensitivity models ------------------------------------------------------------------
# function for fitting the model per stage
fit_sens_model <- function(sens_data) {
  
  model <- lmerTest::lmer(
    temp_sens ~ region * treat_competition + (1 | species),
    data = sens_data
  )
  
  return(model)
}

# fit the models
m_sens_bud_30    <- fit_sens_model(sens_bud_30)
m_sens_flower_30 <- fit_sens_model(sens_flower_30)
m_sens_fruit_30  <- fit_sens_model(sens_fruit_30)
m_sens_seed_30   <- fit_sens_model(sens_seed_30)



# check model outputs -----------------------------------------------------
# bud
summary(m_sens_bud_30)
anova(m_sens_bud_30)
model_performance(m_sens_bud_30)
#check_model(m_sens_bud_30)

# flower
summary(m_sens_flower_30)
anova(m_sens_flower_30)
model_performance(m_sens_flower_30)

# fruit
summary(m_sens_fruit_30)
anova(m_sens_fruit_30)
model_performance(m_sens_fruit_30)


# seeds
summary(m_sens_seed_30)
anova(m_sens_seed_30)
model_performance(m_sens_seed_30)






# predict sensitivity for each stage with function ------------------------

make_sens_predictions <- function(model) {
  
  # 1. new data
  newdat <- expand.grid(
    region = c("Norway", "Switzerland"),
    treat_competition = c("with", "without"),
    temp_sens = 0
  )
  
  # 2. fixed-effect predictions
  newdat$temp_sens <- predict(
    model,
    newdata = newdat,
    re.form = NA
  )
  
  # 3. model matrix
  mm <- model.matrix(terms(model), newdat)
  
  # 4. fixed-effect variance
  pvar <- diag(mm %*% tcrossprod(vcov(model), mm))
  
  # 5. confidence intervals
  cmult <- 2
  
  newdat <- newdat |>
    mutate(
      plo = temp_sens - cmult * sqrt(pvar),
      phi = temp_sens + cmult * sqrt(pvar)
    )
  
  return(newdat)
}


pred_bud_30    <- make_sens_predictions(m_sens_bud_30)
pred_flower_30 <- make_sens_predictions(m_sens_flower_30)
pred_fruit_30  <- make_sens_predictions(m_sens_fruit_30)
pred_seed_30   <- make_sens_predictions(m_sens_seed_30)


plot_predictions_30 <- bind_rows(
  pred_bud_30    |> mutate(stage = "Budding"),
  pred_flower_30 |> mutate(stage = "Flowering"),
  pred_fruit_30  |> mutate(stage = "Fruiting"),
  pred_seed_30   |> mutate(stage = "Seeds")
)
plot_predictions_30



# Plot temp sensitivity of all stages -------------------------------------
ggplot(plot_predictions_30,
       aes(x = treat_competition,
           y = temp_sens,
           color = treat_competition)) +
  
  geom_point(size = 3) +
  
  geom_errorbar(
    aes(ymin = plo, ymax = phi),
    width = 0.15
  ) +
  
  facet_grid(region ~ stage) +
  
  scale_color_manual(values = c(
    "with" = "#528B8B",
    "without" = "#CD950C"
  )) +
  
  labs(
    x = "Biotic interactions",
    y = expression("Predicted temperature sensitivity (days/°C)")
  )+
  theme(legend.position = "none")





temp_sens_30 <- ggplot() +
  
  # raw species values
  geom_jitter(
    data = sens_all_30,
    aes(
      x = treat_competition,
      y = temp_sens,
      color = treat_competition
    ),
    width = 0.08,
    alpha = 0.25
  ) +
  
  # model predictions
  geom_point(
    data = plot_predictions_30,
    aes(
      x = treat_competition,
      y = temp_sens,
      color = treat_competition
    ),
    size = 3
  ) +
  
  geom_errorbar(
    data = plot_predictions_30,
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
  
  labs(
    x = "Biotic interactions",
    y = expression("Temperature sensitivity (days/"*degree*"C)")
  ) +
  theme(
    legend.position = "none"
  )+
  geom_hline(yintercept=0, linetype = "dashed")
temp_sens_30

# ggsave(filename = "Output/Sensitivity/Temperature_sensitivity_30_bud_flower_fruit_seed_onset_NOR_CHE.png", 
#        plot = temp_sens_30,
#        width = 15, height = 10, units = "in")






# temperature sensitivity 14 days window ----------------------------------
# get temp sens per stage for pre_climate1 = 30 days
sens_bud_14    <- get_temp_sens(onset_bud_mean, pre_climate14)
sens_flower_14    <- get_temp_sens(onset_flower_mean, pre_climate14)
sens_fruit_14    <- get_temp_sens(onset_fruit_mean, pre_climate14)
sens_seed_14   <- get_temp_sens(onset_seed_mean, pre_climate14)



# control plots sensitivity -----------------------------------------------
ggplot(sens_bud_14, aes(x = treat_competition, y = temp_sens, color = species)) +
  geom_point(position = position_jitter(width = 0.15, height = 0), size = 2, alpha = 0.85) +
  facet_wrap(~region)


ggplot(sens_flower_14, aes(x = treat_competition, y = temp_sens, color = species)) +
  geom_point(position = position_jitter(width = 0.15, height = 0), size = 2, alpha = 0.85) +
  facet_wrap(~region)


ggplot(sens_fruit_14, aes(x = treat_competition, y = temp_sens, color = species)) +
  geom_point(position = position_jitter(width = 0.15, height = 0), size = 2, alpha = 0.85) +
  facet_wrap(~region)


ggplot(sens_seed_14, aes(x = treat_competition, y = temp_sens, color = species)) +
  geom_point(position = position_jitter(width = 0.15, height = 0), size = 2, alpha = 0.85) +
  facet_wrap(~region)




# combine sens from all stages -------------------------------------------
sens_all_14 <- bind_rows(
  sens_bud_14   |> mutate(stage = "Budding"),
  sens_flower_14|> mutate(stage = "Flowering"),
  sens_fruit_14 |> mutate(stage = "Fruiting"),
  sens_seed_14  |> mutate(stage = "Seeds")
)
sens_all_14



# quick control plot ------------------------------------------------------
ggplot(sens_all_14,
       aes(x = stage,
           y = temp_sens,
           color = treat_competition)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.4),
             alpha = 0.7) 


# plot the raw sens data 
ggplot(sens_all_14,
       aes(x = stage,
           y = temp_sens,
           color = treat_competition)) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  # individual species
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.1,
      dodge.width = 0.4
    ),
    alpha = 0.3
  ) +
  
  # mean ± 95% CI
  stat_summary(
    fun.data = mean_cl_normal,
    geom = "errorbar",
    position = position_dodge(width = 0.4),
    width = 0.15,
    linewidth = 0.8
  ) +
  
  stat_summary(
    fun = mean,
    geom = "point",
    position = position_dodge(width = 0.4),
    size = 4
  ) +
  
  labs(
    x = "Phenological stage",
    y = expression("Temperature sensitivity (days / °C)"),
    color = "Biotic interactions"
  )+
  facet_wrap(~ region)



# sensitivity models ------------------------------------------------------------------
# function for fitting the model per stage
fit_sens_model <- function(sens_data) {
  
  model <- lmerTest::lmer(
    temp_sens ~ region * treat_competition + (1 | species),
    data = sens_data
  )
  
  return(model)
}

# fit the models
m_sens_bud_14    <- fit_sens_model(sens_bud_14)
m_sens_flower_14 <- fit_sens_model(sens_flower_14)
m_sens_fruit_14  <- fit_sens_model(sens_fruit_14)
m_sens_seed_14   <- fit_sens_model(sens_seed_14)



# check model outputs -----------------------------------------------------
# bud
summary(m_sens_bud_14)
anova(m_sens_bud_14)
model_performance(m_sens_bud_14)
#check_model(m_sens_bud_14)

# flower
summary(m_sens_flower_14)
anova(m_sens_flower_14)
model_performance(m_sens_flower_14)

# fruit
summary(m_sens_fruit_14)
anova(m_sens_fruit_14)
model_performance(m_sens_fruit_14)


# seeds
summary(m_sens_seed_14)
anova(m_sens_seed_14)
model_performance(m_sens_seed_14)






# predict sensitivity for each stage with function ------------------------

make_sens_predictions <- function(model) {
  
  # 1. new data
  newdat <- expand.grid(
    region = c("Norway", "Switzerland"),
    treat_competition = c("with", "without"),
    temp_sens = 0
  )
  
  # 2. fixed-effect predictions
  newdat$temp_sens <- predict(
    model,
    newdata = newdat,
    re.form = NA
  )
  
  # 3. model matrix
  mm <- model.matrix(terms(model), newdat)
  
  # 4. fixed-effect variance
  pvar <- diag(mm %*% tcrossprod(vcov(model), mm))
  
  # 5. confidence intervals
  cmult <- 2
  
  newdat <- newdat |>
    mutate(
      plo = temp_sens - cmult * sqrt(pvar),
      phi = temp_sens + cmult * sqrt(pvar)
    )
  
  return(newdat)
}


pred_bud_14    <- make_sens_predictions(m_sens_bud_14)
pred_flower_14 <- make_sens_predictions(m_sens_flower_14)
pred_fruit_14  <- make_sens_predictions(m_sens_fruit_14)
pred_seed_14   <- make_sens_predictions(m_sens_seed_14)


plot_predictions_14 <- bind_rows(
  pred_bud_14    |> mutate(stage = "Budding"),
  pred_flower_14 |> mutate(stage = "Flowering"),
  pred_fruit_14  |> mutate(stage = "Fruiting"),
  pred_seed_14   |> mutate(stage = "Seeds")
)
plot_predictions_14



# Plot temp sensitivity of all stages -------------------------------------
ggplot(plot_predictions_14,
       aes(x = treat_competition,
           y = temp_sens,
           color = treat_competition)) +
  
  geom_point(size = 3) +
  
  geom_errorbar(
    aes(ymin = plo, ymax = phi),
    width = 0.15
  ) +
  
  facet_grid(region ~ stage) +
  
  scale_color_manual(values = c(
    "with" = "#528B8B",
    "without" = "#CD950C"
  )) +
  
  labs(
    x = "Biotic interactions",
    y = expression("Predicted temperature sensitivity (days/°C)")
  )+
  theme(legend.position = "none")





temp_sens_14 <- ggplot() +
  
  # raw species values
  geom_jitter(
    data = sens_all_14,
    aes(
      x = treat_competition,
      y = temp_sens,
      color = treat_competition
    ),
    width = 0.08,
    alpha = 0.25
  ) +
  
  # model predictions
  geom_point(
    data = plot_predictions_14,
    aes(
      x = treat_competition,
      y = temp_sens,
      color = treat_competition
    ),
    size = 3
  ) +
  
  geom_errorbar(
    data = plot_predictions_14,
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
  
  labs(
    x = "Biotic interactions",
    y = expression("Temperature sensitivity (days/"*degree*"C)")
  ) +
  theme(
    legend.position = "none"
  )+
  geom_hline(yintercept=0, linetype = "dashed")
temp_sens_14

# ggsave(filename = "Output/Sensitivity/Temperature_sensitivity_14_bud_flower_fruit_seed_onset_NOR_CHE.png", 
#       plot = temp_sens_14,
#        width = 15, height = 10, units = "in")






# temperature sensitivity 60 days window ----------------------------------
# get temp sens per stage for pre_climate1 = 60 days
sens_bud_60    <- get_temp_sens(onset_bud_mean, pre_climate60)
sens_flower_60    <- get_temp_sens(onset_flower_mean, pre_climate60)
sens_fruit_60    <- get_temp_sens(onset_fruit_mean, pre_climate60)
sens_seed_60   <- get_temp_sens(onset_seed_mean, pre_climate60)



# control plots sensitivity -----------------------------------------------
ggplot(sens_bud_60, aes(x = treat_competition, y = temp_sens, color = species)) +
  geom_point(position = position_jitter(width = 0.15, height = 0), size = 2, alpha = 0.85) +
  facet_wrap(~region)


ggplot(sens_flower_60, aes(x = treat_competition, y = temp_sens, color = species)) +
  geom_point(position = position_jitter(width = 0.15, height = 0), size = 2, alpha = 0.85) +
  facet_wrap(~region)


ggplot(sens_fruit_60, aes(x = treat_competition, y = temp_sens, color = species)) +
  geom_point(position = position_jitter(width = 0.15, height = 0), size = 2, alpha = 0.85) +
  facet_wrap(~region)


ggplot(sens_seed_60, aes(x = treat_competition, y = temp_sens, color = species)) +
  geom_point(position = position_jitter(width = 0.15, height = 0), size = 2, alpha = 0.85) +
  facet_wrap(~region)




# combine sens from all stages -------------------------------------------
sens_all_60 <- bind_rows(
  sens_bud_60   |> mutate(stage = "Budding"),
  sens_flower_60|> mutate(stage = "Flowering"),
  sens_fruit_60 |> mutate(stage = "Fruiting"),
  sens_seed_60  |> mutate(stage = "Seeds")
)
sens_all_60



# quick control plot ------------------------------------------------------
ggplot(sens_all_60,
       aes(x = stage,
           y = temp_sens,
           color = treat_competition)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.4),
             alpha = 0.7) 


# plot the raw sens data 
ggplot(sens_all_60,
       aes(x = stage,
           y = temp_sens,
           color = treat_competition)) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  # individual species
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.1,
      dodge.width = 0.4
    ),
    alpha = 0.3
  ) +
  
  # mean ± 95% CI
  stat_summary(
    fun.data = mean_cl_normal,
    geom = "errorbar",
    position = position_dodge(width = 0.4),
    width = 0.15,
    linewidth = 0.8
  ) +
  
  stat_summary(
    fun = mean,
    geom = "point",
    position = position_dodge(width = 0.4),
    size = 4
  ) +
  
  labs(
    x = "Phenological stage",
    y = expression("Temperature sensitivity (days / °C)"),
    color = "Biotic interactions"
  )+
  facet_wrap(~ region)



# sensitivity models ------------------------------------------------------------------
# function for fitting the model per stage
fit_sens_model <- function(sens_data) {
  
  model <- lmerTest::lmer(
    temp_sens ~ region * treat_competition + (1 | species),
    data = sens_data
  )
  
  return(model)
}

# fit the models
m_sens_bud_60    <- fit_sens_model(sens_bud_60)
m_sens_flower_60 <- fit_sens_model(sens_flower_60)
m_sens_fruit_60  <- fit_sens_model(sens_fruit_60)
m_sens_seed_60   <- fit_sens_model(sens_seed_60)



# check model outputs -----------------------------------------------------
# bud
summary(m_sens_bud_60)
anova(m_sens_bud_60)
model_performance(m_sens_bud_60)
#check_model(m_sens_bud_60)

# flower
summary(m_sens_flower_60)
anova(m_sens_flower_60)
model_performance(m_sens_flower_60)

# fruit
summary(m_sens_fruit_60)
anova(m_sens_fruit_60)
model_performance(m_sens_fruit_60)


# seeds
summary(m_sens_seed_60)
anova(m_sens_seed_60)
model_performance(m_sens_seed_60)






# predict sensitivity for each stage with function ------------------------

make_sens_predictions <- function(model) {
  
  # 1. new data
  newdat <- expand.grid(
    region = c("Norway", "Switzerland"),
    treat_competition = c("with", "without"),
    temp_sens = 0
  )
  
  # 2. fixed-effect predictions
  newdat$temp_sens <- predict(
    model,
    newdata = newdat,
    re.form = NA
  )
  
  # 3. model matrix
  mm <- model.matrix(terms(model), newdat)
  
  # 4. fixed-effect variance
  pvar <- diag(mm %*% tcrossprod(vcov(model), mm))
  
  # 5. confidence intervals
  cmult <- 2
  
  newdat <- newdat |>
    mutate(
      plo = temp_sens - cmult * sqrt(pvar),
      phi = temp_sens + cmult * sqrt(pvar)
    )
  
  return(newdat)
}


pred_bud_60    <- make_sens_predictions(m_sens_bud_60)
pred_flower_60 <- make_sens_predictions(m_sens_flower_60)
pred_fruit_60  <- make_sens_predictions(m_sens_fruit_60)
pred_seed_60   <- make_sens_predictions(m_sens_seed_60)


plot_predictions_60 <- bind_rows(
  pred_bud_60    |> mutate(stage = "Budding"),
  pred_flower_60 |> mutate(stage = "Flowering"),
  pred_fruit_60  |> mutate(stage = "Fruiting"),
  pred_seed_60   |> mutate(stage = "Seeds")
)
plot_predictions_60



# Plot temp sensitivity of all stages -------------------------------------
ggplot(plot_predictions_60,
       aes(x = treat_competition,
           y = temp_sens,
           color = treat_competition)) +
  
  geom_point(size = 3) +
  
  geom_errorbar(
    aes(ymin = plo, ymax = phi),
    width = 0.15
  ) +
  
  facet_grid(region ~ stage) +
  
  scale_color_manual(values = c(
    "with" = "#528B8B",
    "without" = "#CD950C"
  )) +
  
  labs(
    x = "Biotic interactions",
    y = expression("Predicted temperature sensitivity (days/°C)")
  )+
  theme(legend.position = "none")





temp_sens_60 <- ggplot() +
  
  # raw species values
  geom_jitter(
    data = sens_all_60,
    aes(
      x = treat_competition,
      y = temp_sens,
      color = treat_competition
    ),
    width = 0.08,
    alpha = 0.25
  ) +
  
  # model predictions
  geom_point(
    data = plot_predictions_60,
    aes(
      x = treat_competition,
      y = temp_sens,
      color = treat_competition
    ),
    size = 3
  ) +
  
  geom_errorbar(
    data = plot_predictions_60,
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
  
  labs(
    x = "Biotic interactions",
    y = expression("Temperature sensitivity (days/"*degree*"C)")
  ) +
  theme(
    legend.position = "none"
  )+
  geom_hline(yintercept=0, linetype = "dashed")
temp_sens_60

# ggsave(filename = "Output/Sensitivity/Temperature_sensitivity_60_bud_flower_fruit_seed_onset_NOR_CHE.png", 
#       plot = temp_sens_60,
#        width = 15, height = 10, units = "in")



# temperature sensitivity growing season mean ----------------------------------
# get temp sens per stage for pre_climate1 = 60 days
sens_bud_gs   <- get_temp_sens(onset_bud_mean, pre_climate_gs)
sens_flower_gs    <- get_temp_sens(onset_flower_mean, pre_climate_gs)
sens_fruit_gs    <- get_temp_sens(onset_fruit_mean, pre_climate_gs)
sens_seed_gs   <- get_temp_sens(onset_seed_mean, pre_climate_gs)



# control plots sensitivity -----------------------------------------------
ggplot(sens_bud_gs, aes(x = treat_competition, y = temp_sens, color = species)) +
  geom_point(position = position_jitter(width = 0.15, height = 0), size = 2, alpha = 0.85) +
  facet_wrap(~region)


ggplot(sens_flower_gs, aes(x = treat_competition, y = temp_sens, color = species)) +
  geom_point(position = position_jitter(width = 0.15, height = 0), size = 2, alpha = 0.85) +
  facet_wrap(~region)


ggplot(sens_fruit_gs, aes(x = treat_competition, y = temp_sens, color = species)) +
  geom_point(position = position_jitter(width = 0.15, height = 0), size = 2, alpha = 0.85) +
  facet_wrap(~region)


ggplot(sens_seed_gs, aes(x = treat_competition, y = temp_sens, color = species)) +
  geom_point(position = position_jitter(width = 0.15, height = 0), size = 2, alpha = 0.85) +
  facet_wrap(~region)




# combine sens from all stages -------------------------------------------
sens_all_gs <- bind_rows(
  sens_bud_gs   |> mutate(stage = "Budding"),
  sens_flower_gs|> mutate(stage = "Flowering"),
  sens_fruit_gs |> mutate(stage = "Fruiting"),
  sens_seed_gs  |> mutate(stage = "Seeds")
)
sens_all_gs



# quick control plot ------------------------------------------------------
ggplot(sens_all_gs,
       aes(x = stage,
           y = temp_sens,
           color = treat_competition)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.4),
             alpha = 0.7) 


# plot the raw sens data 
ggplot(sens_all_gs,
       aes(x = stage,
           y = temp_sens,
           color = treat_competition)) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  # individual species
  geom_point(
    position = position_jitterdodge(
      jitter.width = 0.1,
      dodge.width = 0.4
    ),
    alpha = 0.3
  ) +
  
  # mean ± 95% CI
  stat_summary(
    fun.data = mean_cl_normal,
    geom = "errorbar",
    position = position_dodge(width = 0.4),
    width = 0.15,
    linewidth = 0.8
  ) +
  
  stat_summary(
    fun = mean,
    geom = "point",
    position = position_dodge(width = 0.4),
    size = 4
  ) +
  
  labs(
    x = "Phenological stage",
    y = expression("Temperature sensitivity (days / °C)"),
    color = "Biotic interactions"
  )+
  facet_wrap(~ region)



# sensitivity models ------------------------------------------------------------------
# function for fitting the model per stage
fit_sens_model <- function(sens_data) {
  
  model <- lmerTest::lmer(
    temp_sens ~ region * treat_competition + (1 | species),
    data = sens_data
  )
  
  return(model)
}

# fit the models
m_sens_bud_gs    <- fit_sens_model(sens_bud_gs)
m_sens_flower_gs <- fit_sens_model(sens_flower_gs)
m_sens_fruit_gs  <- fit_sens_model(sens_fruit_gs)
m_sens_seed_gs   <- fit_sens_model(sens_seed_gs)



# check model outputs -----------------------------------------------------
# bud
summary(m_sens_bud_gs)
anova(m_sens_bud_gs)
model_performance(m_sens_bud_gs)
#check_model(m_sens_bud_gs)

# flower
summary(m_sens_flower_gs)
anova(m_sens_flower_gs)
model_performance(m_sens_flower_gs)

# fruit
summary(m_sens_fruit_gs)
anova(m_sens_fruit_gs)
model_performance(m_sens_fruit_gs)


# seeds
summary(m_sens_seed_gs)
anova(m_sens_seed_gs)
model_performance(m_sens_seed_gs)






# predict sensitivity for each stage with function ------------------------

make_sens_predictions <- function(model) {
  
  # 1. new data
  newdat <- expand.grid(
    region = c("Norway", "Switzerland"),
    treat_competition = c("with", "without"),
    temp_sens = 0
  )
  
  # 2. fixed-effect predictions
  newdat$temp_sens <- predict(
    model,
    newdata = newdat,
    re.form = NA
  )
  
  # 3. model matrix
  mm <- model.matrix(terms(model), newdat)
  
  # 4. fixed-effect variance
  pvar <- diag(mm %*% tcrossprod(vcov(model), mm))
  
  # 5. confidence intervals
  cmult <- 2
  
  newdat <- newdat |>
    mutate(
      plo = temp_sens - cmult * sqrt(pvar),
      phi = temp_sens + cmult * sqrt(pvar)
    )
  
  return(newdat)
}


pred_bud_gs    <- make_sens_predictions(m_sens_bud_gs)
pred_flower_gs <- make_sens_predictions(m_sens_flower_gs)
pred_fruit_gs  <- make_sens_predictions(m_sens_fruit_gs)
pred_seed_gs   <- make_sens_predictions(m_sens_seed_gs)


plot_predictions_gs <- bind_rows(
  pred_bud_gs    |> mutate(stage = "Budding"),
  pred_flower_gs |> mutate(stage = "Flowering"),
  pred_fruit_gs  |> mutate(stage = "Fruiting"),
  pred_seed_gs   |> mutate(stage = "Seeds")
)
plot_predictions_gs



# Plot temp sensitivity of all stages -------------------------------------
ggplot(plot_predictions_gs,
       aes(x = treat_competition,
           y = temp_sens,
           color = treat_competition)) +
  
  geom_point(size = 3) +
  
  geom_errorbar(
    aes(ymin = plo, ymax = phi),
    width = 0.15
  ) +
  
  facet_grid(region ~ stage) +
  
  scale_color_manual(values = c(
    "with" = "#528B8B",
    "without" = "#CD950C"
  )) +
  
  labs(
    x = "Biotic interactions",
    y = expression("Predicted temperature sensitivity (days/°C)")
  )+
  theme(legend.position = "none")





temp_sens_gs <- ggplot() +
  
  # raw species values
  geom_jitter(
    data = sens_all_gs,
    aes(
      x = treat_competition,
      y = temp_sens,
      color = treat_competition
    ),
    width = 0.08,
    alpha = 0.25
  ) +
  
  # model predictions
  geom_point(
    data = plot_predictions_gs,
    aes(
      x = treat_competition,
      y = temp_sens,
      color = treat_competition
    ),
    size = 3
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
  
  labs(
    x = "Biotic interactions",
    y = expression("Temperature sensitivity (days/"*degree*"C)")
  ) +
  theme(
    legend.position = "none"
  )+
  geom_hline(yintercept=0, linetype = "dashed")
temp_sens_gs

# ggsave(filename = "Output/Sensitivity/Temperature_sensitivity_gs_bud_flower_fruit_seed_onset_NOR_CHE.png", 
#       plot = temp_sens_gs,
#        width = 15, height = 10, units = "in")


# which pre-climate is best -----------------------------------------------

AIC(m_sens_bud_14, m_sens_bud_30, m_sens_bud_60, m_sens_bud_gs)

AIC(m_sens_flower_14, m_sens_flower_30, m_sens_flower_60, m_sens_flower_gs)

AIC(m_sens_fruit_14, m_sens_fruit_30, m_sens_fruit_60, m_sens_fruit_gs)

AIC(m_sens_seed_14, m_sens_seed_30, m_sens_seed_60, m_sens_seed_gs)


# looks like the 60 day pre-climate window is the best



