


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
climate_nor_23_pre1 <- climate_23_daily |>
  filter(
    date >= as.Date("2023-05-14"),
    date <= as.Date("2023-06-12")
  )

climate_nor_23_pre1 <- climate_23_daily |>
  filter(
    date >= as.Date("2023-05-21"),
    date <= as.Date("2023-06-20")
  )

climate_nor_23_pre1 <- climate_nor_23_pre1 |> 
  mutate(region = "Norway")

climate_che_22_pre1 <- climate_che_combined |>
  filter(
    date >= as.Date("2022-04-29"),
    date <= as.Date("2022-05-28")
  )

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
onset_fruit  <- get_mean_onset(phenology3, "No_Fruit")
onset_seed   <- get_mean_onset(phenology3, "No_Seed")

# combine the fruiting stages nor and che to compare the onset
# this is to be taken with caution because the stages are not the same
# but for the onset it could be comparable


# Average onset per site, species and treatment --------------------------------------------------
# could do per plot as well
onset_bud_mean <- onset_bud |>
  group_by(region, site, year, treat_competition, species) |>
  summarise(onset = mean(onset, na.rm = TRUE),
            .groups = "drop")
onset_bud_mean

onset_flower_mean <- onset_flower |>
  group_by(region, site, year, treat_competition, species) |>
  summarise(onset = mean(onset, na.rm = TRUE),
            .groups = "drop")
onset_flower_mean

onset_fruit_mean <- onset_fruit |>
  group_by(region, site, year, treat_competition, species) |>
  summarise(onset = mean(onset, na.rm = TRUE),
            .groups = "drop")
onset_fruit_mean

onset_seed_mean <- onset_seed |>
  group_by(region, site, year, treat_competition, species) |>
  summarise(onset = mean(onset, na.rm = TRUE),
            .groups = "drop")
onset_seed_mean



#


# Calculate temperature sensitivity ---------------------------------------


# bud ---------------------------------------------------------------------
sens_bud <- onset_bud_mean |>
  left_join(pre_climate1, by = c("region","site")) |>
  group_by(site, species, ) |> 
  pivot_wider(names_from = site,
              values_from = c(onset, Tmean)) |>
  mutate(
    temp_sens = (onset_lo - onset_hi) / (Tmean_lo - Tmean_hi)
  )
sens_bud


ggplot(sens_bud, aes(x = treat_competition, y = temp_sens, color = species)) +
  geom_point(position = position_jitter(width = 0.15, height = 0), size = 2, alpha = 0.85) +
  theme_bw()+
  facet_wrap(~region)













