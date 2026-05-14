


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

# NOR
climate_nor_23_pre1 <- climate_23_daily |>
  filter(
    date >= as.Date("2023-05-21"),
    date <= as.Date("2023-06-20")
  )

climate_nor_23_pre1 <- climate_nor_23_pre1 |> 
  mutate(region = "Norway")

# CHE
climate_che_22_pre1 <- climate_che_combined |>
  filter(
    date >= as.Date("2022-05-05"),
    date <= as.Date("2022-06-04")
  )

climate_che_22_pre1 <- climate_che_22_pre1 |> 
  mutate(region = "Switzerland")

# combine
climate_all1 <- bind_rows(climate_che_22_pre1, climate_nor_23_pre1)


# calculate mean per site in this time period
pre_climate1 <-  climate_all1 |>
  group_by(region, site) |>
  summarise(Tmean = mean(Tavg))
pre_climate1



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
sens_bud_30    <- get_temp_sens(onset_bud_mean, pre_climate1)
sens_flower_30    <- get_temp_sens(onset_flower_mean, pre_climate1)
sens_fruit_30    <- get_temp_sens(onset_fruit_mean, pre_climate1)
sens_seed_30   <- get_temp_sens(onset_seed_mean, pre_climate1)



# control plots sensitivity -----------------------------------------------
ggplot(sens_bud, aes(x = treat_competition, y = temp_sens, color = species)) +
  geom_point(position = position_jitter(width = 0.15, height = 0), size = 2, alpha = 0.85) +
  facet_wrap(~region)


ggplot(sens_flower, aes(x = treat_competition, y = temp_sens, color = species)) +
  geom_point(position = position_jitter(width = 0.15, height = 0), size = 2, alpha = 0.85) +
  facet_wrap(~region)


ggplot(sens_fruit, aes(x = treat_competition, y = temp_sens, color = species)) +
  geom_point(position = position_jitter(width = 0.15, height = 0), size = 2, alpha = 0.85) +
  facet_wrap(~region)


ggplot(sens_seed, aes(x = treat_competition, y = temp_sens, color = species)) +
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





ggplot() +
  
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




