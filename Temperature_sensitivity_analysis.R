
# 02_temp_sens ------------------------------------------------------------

# Temperature sensitivity analysis -------------------------------------------------


# RangeX phenology effect of transplantation on temperature sensitivity ------------

## Data used: RangeX_clean_phenology_2023_NOR.csv
##            
##            RangeX_clean_MetadataFocal_CHE.csv
##            RangeX_metadata_focal_NOR.csv
##            RangeX_clean_climate_station_NOR_2021-2025.csv
## Date:      27.05.26
## Author:    Nadine Arzt
## Purpose:   Effect of transplantation on temperature sensitivity


# load library ------------------------------------------------------------
library(conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
library(tidyverse)
library(lme4)
library(performance)
library(see)
library(emmeans)
library(ggeffects)


# source script with functions --------------------------------------------
source("Temperature_sensitivity_functions.R")


# calculate delta T -------------------------------------------------------
# use growing season mean because this is not specific to species
# we have species with different onset timings

# source climate scripts --------------------------------------------------
source("Data_preparation_climate_station_NOR.R")
climate_23 
climate_23_daily

source("Data_preparation_climate_station_CHE.R")
climate_che_combined
# which is daily per site



# mean temperature growing season mean ------------------------------------------
# NOR
climate_nor_23_gs <- climate_23_daily |>
  filter(
    date >= as.Date("2023-04-01"),
    date <= as.Date("2023-09-30")
  ) |> 
  mutate(region = "Norway")


# CHE
climate_che_22_gs <- climate_che_combined |>
  filter(
    date >= as.Date("2022-04-01"),
    date <= as.Date("2022-09-30")
  ) |> 
  mutate(region = "Switzerland")


# combine
climate_all_gs <- bind_rows(climate_che_22_gs, climate_nor_23_gs)


# calculate mean per site in this time period
temperature_mean_gs <-  climate_all_gs |>
  group_by(region, site) |>
  summarise(Tmean = mean(Tavg))
temperature_mean_gs




# source clean phenology data -----------------------------------------------
source("Data_preparation_phenology_NOR_CHE_combined.R")

# use this data set
names(phenology)



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
# not so necessary because we use growing season mean temperature
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
  group_by(region, site, species) |>
  summarise(mean_onset = mean(onset, na.rm = TRUE),
            .groups = "drop")
mean_onset_bud_species

mean_onset_bud <- onset_bud |>
  group_by(region, site) |>
  summarise(mean_onset = mean(onset, na.rm = TRUE),
            .groups = "drop")
mean_onset_bud

# 1 Norway      hi     2023       181. # 30.06
# 2 Norway      lo     2023       164. # 13.06
# 3 Switzerland hi     2022       165. # 14.06
# 4 Switzerland lo     2022       149. # 29.05

mean_onset_bud_region <- onset_bud |>
  group_by(region) |>
  summarise(mean_onset = mean(onset, na.rm = TRUE),
            .groups = "drop")
mean_onset_bud_region



# calculate mean onset per species and plot for all stages ----------------
onset_bud    <- get_mean_onset(phenology3, "No_Buds", "jday")
onset_flower <- get_mean_onset(phenology3, "No_FloOpen", "jday")
onset_fruit  <- get_mean_onset(phenology3, "No_FloWithrd", "jday")
onset_seed   <- get_mean_onset(phenology3, "No_Seeds", "jday")

# combine the fruiting stages nor and che to compare the onset
# this is to be taken with caution because the stages are not the same
# but for the onset it could be comparable





# Average across species, site and treat -----------------------------------------
onset_bud_mean    <- get_species_onset(onset_bud)
onset_flower_mean <- get_species_onset(onset_flower)
onset_fruit_mean  <- get_species_onset(onset_fruit)
onset_seed_mean   <- get_species_onset(onset_seed)
# this doesn't change the values


# Calculate temperature sensitivity ---------------------------------------
# growing season mean ----------------------------------
# get temp sens per stage while using the same delta T 
sens_bud_gs   <- get_temp_sens(onset_bud_mean, temperature_mean_gs)
sens_flower_gs    <- get_temp_sens(onset_flower_mean, temperature_mean_gs)
sens_fruit_gs    <- get_temp_sens(onset_fruit_mean, temperature_mean_gs)
sens_seed_gs   <- get_temp_sens(onset_seed_mean, temperature_mean_gs)




# control plot sensitivity ------------------------------------------------
p_bud <- control_plot_temp_sens(sens_bud_gs)
p_bud
p_flower <- control_plot_temp_sens(sens_flower_gs)
p_flower
p_fruit <- control_plot_temp_sens(sens_fruit_gs)
p_fruit
p_seed <- control_plot_temp_sens(sens_seed_gs)
p_seed



# combine sens from all stages -------------------------------------------
sens_all_gs <- bind_rows(
  sens_bud_gs   |> mutate(stage = "Budding"),
  sens_flower_gs|> mutate(stage = "Flowering"),
  sens_fruit_gs |> mutate(stage = "Fruiting"),
  sens_seed_gs  |> mutate(stage = "Seeds")
)
sens_all_gs


# quick control plot ------------------------------------------------------
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

# NOR ---------------------------------------------------------------------
# fit the models per stage for Norway
m_sens_bud_gs_nor    <- fit_sens_model(sens_bud_gs, "Norway")
m_sens_flower_gs_nor <- fit_sens_model(sens_flower_gs, "Norway")
m_sens_fruit_gs_nor  <- fit_sens_model(sens_fruit_gs, "Norway")
m_sens_seed_gs_nor   <- fit_sens_model(sens_seed_gs, "Norway")

# check model output
summary(m_sens_bud_gs_nor)
anova(m_sens_bud_gs_nor)
model_performance(m_sens_bud_gs_nor)
#check_model(m_sens_bud_gs_nor)

summary(m_sens_flower_gs_nor)
anova(m_sens_flower_gs_nor)
model_performance(m_sens_flower_gs_nor)
#check_model(m_sens_flower_gs_nor)

summary(m_sens_fruit_gs_nor)
anova(m_sens_fruit_gs_nor)
model_performance(m_sens_fruit_gs_nor)
#check_model(m_sens_fruit_gs_nor)

summary(m_sens_seed_gs_nor)
anova(m_sens_seed_gs_nor)
model_performance(m_sens_seed_gs_nor)
#check_model(m_sens_seed_gs_nor)


m_sens_bud_gs_nor    <- fit_sens_model2(sens_bud_gs, "Norway")


# CHE ---------------------------------------------------------------------
# fit the models per stage for Switzerland
m_sens_bud_gs_che    <- fit_sens_model(sens_bud_gs, "Switzerland")
m_sens_flower_gs_che <- fit_sens_model(sens_flower_gs, "Switzerland")
m_sens_fruit_gs_che  <- fit_sens_model(sens_fruit_gs, "Switzerland")
m_sens_seed_gs_che   <- fit_sens_model(sens_seed_gs, "Switzerland")

# check model output
summary(m_sens_bud_gs_che)
anova(m_sens_bud_gs_che)
model_performance(m_sens_bud_gs_che)
#check_model(m_sens_bud_gs_che)

# the models have a singularity issue
# and species has no variation
# therfore no R2 marginal
VarCorr(m_sens_bud_gs_che)

summary(m_sens_flower_gs_che)
anova(m_sens_flower_gs_che)
model_performance(m_sens_flower_gs_che)
#check_model(m_sens_flower_gs_che)

summary(m_sens_fruit_gs_che)
anova(m_sens_fruit_gs_che)
model_performance(m_sens_fruit_gs_che)
#check_model(m_sens_fruit_gs_che)

summary(m_sens_seed_gs_che)
anova(m_sens_seed_gs_che)
model_performance(m_sens_seed_gs_che)
#check_model(m_sens_seed_gs_che)



# predict sensitivity for each stage with function ------------------------
# NOR ---------------------------------------------------------------------
# manual
pred_bud_gs_nor    <- make_sens_predictions(m_sens_bud_gs_nor)
pred_flower_gs_nor <- make_sens_predictions(m_sens_flower_gs_nor)
pred_fruit_gs_nor  <- make_sens_predictions(m_sens_fruit_gs_nor)
pred_seed_gs_nor   <- make_sens_predictions(m_sens_seed_gs_nor)


# ggpredict
pred_bud_gs_nor2    <- make_sens_predictions2(m_sens_bud_gs_nor)

# predict
pred_bud_gs_nor3    <- make_sens_predictions3(m_sens_bud_gs_nor)



# CHE ---------------------------------------------------------------------
# manual
pred_bud_gs_che    <- make_sens_predictions(m_sens_bud_gs_che)
pred_flower_gs_che <- make_sens_predictions(m_sens_flower_gs_che)
pred_fruit_gs_che  <- make_sens_predictions(m_sens_fruit_gs_che)
pred_seed_gs_che   <- make_sens_predictions(m_sens_seed_gs_che)



# combine all predictions -------------------------------------------------
plot_predictions_gs <- bind_rows(
  pred_bud_gs_nor    |> mutate(stage = "Budding", region = "Norway"),
  pred_flower_gs_nor |> mutate(stage = "Flowering", region = "Norway"),
  pred_fruit_gs_nor  |> mutate(stage = "Fruiting", region = "Norway"),
  pred_seed_gs_nor   |> mutate(stage = "Seeds", region = "Norway"),
  pred_bud_gs_che    |> mutate(stage = "Budding", region = "Switzerland"),
  pred_flower_gs_che |> mutate(stage = "Flowering", region = "Switzerland"),
  pred_fruit_gs_che  |> mutate(stage = "Fruiting", region = "Switzerland"),
  pred_seed_gs_che   |> mutate(stage = "Seeds", region = "Switzerland")
)
plot_predictions_gs





# GDD ---------------------------------------------------------------------

# source script where GDD is combined with phenology ----------------------
source("GDD_Cooling_bud_flower_fruit_predictions.R")
# use 
phenology_gdd_nor_che


# calculate mean onset per species and plot for all stages ----------------
onset_bud_gdd    <- get_mean_onset(phenology_gdd_nor_che, "No_Buds", "GDD_cum")
onset_flower_gdd <- get_mean_onset(phenology_gdd_nor_che, "No_FloOpen", "GDD_cum")
onset_fruit_gdd  <- get_mean_onset(phenology_gdd_nor_che, "No_FloWithrd", "GDD_cum")
onset_seed_gdd   <- get_mean_onset(phenology_gdd_nor_che, "No_Seeds", "GDD_cum")





# Average across species, site and treat -----------------------------------------
onset_bud_mean_gdd    <- get_species_onset(onset_bud_gdd)
onset_flower_mean_gdd <- get_species_onset(onset_flower_gdd)
onset_fruit_mean_gdd  <- get_species_onset(onset_fruit_gdd)
onset_seed_mean_gdd   <- get_species_onset(onset_seed_gdd)
# this doesn't change the values


# Calculate temperature sensitivity ---------------------------------------
# growing season mean ----------------------------------
# get temp sens per stage while using the same delta T 
sens_bud_gs_gdd   <- get_temp_sens(onset_bud_mean_gdd, temperature_mean_gs)
sens_flower_gs_gdd    <- get_temp_sens(onset_flower_mean_gdd, temperature_mean_gs)
sens_fruit_gs_gdd    <- get_temp_sens(onset_fruit_mean_gdd, temperature_mean_gs)
sens_seed_gs_gdd   <- get_temp_sens(onset_seed_mean_gdd, temperature_mean_gs)




# control plot sensitivity ------------------------------------------------
p_bud_gdd <- control_plot_temp_sens(sens_bud_gs_gdd)
p_bud_gdd 
p_flower_gdd  <- control_plot_temp_sens(sens_flower_gs_gdd)
p_flower_gdd 
p_fruit_gdd  <- control_plot_temp_sens(sens_fruit_gs_gdd)
p_fruit_gdd 
p_seed_gdd  <- control_plot_temp_sens(sens_seed_gs_gdd)
p_seed_gdd 



# combine sens from all stages -------------------------------------------
sens_all_gs_gdd <- bind_rows(
  sens_bud_gs_gdd   |> mutate(stage = "Budding"),
  sens_flower_gs_gdd|> mutate(stage = "Flowering"),
  sens_fruit_gs_gdd |> mutate(stage = "Fruiting"),
  sens_seed_gs_gdd  |> mutate(stage = "Seeds")
)
sens_all_gs_gdd


# quick control plot ------------------------------------------------------
# plot the raw sens data 
ggplot(sens_all_gs_gdd,
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

# NOR ---------------------------------------------------------------------
# fit the models per stage for Norway
m_sens_bud_gs_gdd_nor    <- fit_sens_model(sens_bud_gs_gdd, "Norway")
m_sens_flower_gs_gdd_nor <- fit_sens_model(sens_flower_gs_gdd, "Norway")
m_sens_fruit_gs_gdd_nor  <- fit_sens_model(sens_fruit_gs_gdd, "Norway")
m_sens_seed_gs_gdd_nor   <- fit_sens_model(sens_seed_gs_gdd, "Norway")

# check model output
summary(m_sens_bud_gs_gdd_nor)
anova(m_sens_bud_gs_gdd_nor)
model_performance(m_sens_bud_gs_gdd_nor)
#check_model(m_sens_bud_gs_gdd_nor)

summary(m_sens_flower_gs_gdd_nor)
anova(m_sens_flower_gs_gdd_nor)
model_performance(m_sens_flower_gs_gdd_nor)
#check_model(m_sens_flower_gs_gdd_nor)

summary(m_sens_fruit_gs_gdd_nor)
anova(m_sens_fruit_gs_gdd_nor)
model_performance(m_sens_fruit_gs_gdd_nor)
#check_model(m_sens_fruit_gs_gdd_nor)

summary(m_sens_seed_gs_gdd_nor)
anova(m_sens_seed_gs_gdd_nor)
model_performance(m_sens_seed_gs_gdd_nor)
#check_model(m_sens_seed_gs_gdd_nor)




# CHE ---------------------------------------------------------------------
# fit the models per stage for Switzerland
m_sens_bud_gs_gdd_che    <- fit_sens_model(sens_bud_gs_gdd, "Switzerland")
m_sens_flower_gs_gdd_che <- fit_sens_model(sens_flower_gs_gdd, "Switzerland")
m_sens_fruit_gs_gdd_che  <- fit_sens_model(sens_fruit_gs_gdd, "Switzerland")
m_sens_seed_gs_gdd_che   <- fit_sens_model(sens_seed_gs_gdd, "Switzerland")

# check model output
summary(m_sens_bud_gs_gdd_che)
anova(m_sens_bud_gs_gdd_che)
model_performance(m_sens_bud_gs_gdd_che)
#check_model(m_sens_bud_gs_gdd_che)

summary(m_sens_flower_gs_gdd_che)
anova(m_sens_flower_gs_gdd_che)
model_performance(m_sens_flower_gs_gdd_che)
#check_model(m_sens_flower_gs_gdd_che)

summary(m_sens_fruit_gs_gdd_che)
anova(m_sens_fruit_gs_gdd_che)
model_performance(m_sens_fruit_gs_gdd_che)
#check_model(m_sens_fruit_gs_gdd_che)

summary(m_sens_seed_gs_gdd_che)
anova(m_sens_seed_gs_gdd_che)
model_performance(m_sens_seed_gs_gdd_che)
#check_model(m_sens_seed_gs_gdd_che)




# predict sensitivity for each stage with function ------------------------
# NOR ---------------------------------------------------------------------
# manual
pred_bud_gs_gdd_nor    <- make_sens_predictions(m_sens_bud_gs_gdd_nor)
pred_flower_gs_gdd_nor <- make_sens_predictions(m_sens_flower_gs_gdd_nor)
pred_fruit_gs_gdd_nor  <- make_sens_predictions(m_sens_fruit_gs_gdd_nor)
pred_seed_gs_gdd_nor   <- make_sens_predictions(m_sens_seed_gs_gdd_nor)


# ggpredict
pred_bud_gs_gdd_nor2    <- make_sens_predictions2(m_sens_bud_gs_gdd_nor)
pred_flower_gs_gdd_nor2 <- make_sens_predictions2(m_sens_flower_gs_gdd_nor)
pred_fruit_gs_gdd_nor2  <- make_sens_predictions2(m_sens_fruit_gs_gdd_nor)
pred_seed_gs_gdd_nor2   <- make_sens_predictions2(m_sens_seed_gs_gdd_nor)





# CHE ---------------------------------------------------------------------
# manual
pred_bud_gs_gdd_che    <- make_sens_predictions(m_sens_bud_gs_gdd_che)
pred_flower_gs_gdd_che <- make_sens_predictions(m_sens_flower_gs_gdd_che)
pred_fruit_gs_gdd_che  <- make_sens_predictions(m_sens_fruit_gs_gdd_che)
pred_seed_gs_gdd_che   <- make_sens_predictions(m_sens_seed_gs_gdd_che)

# ggpredict
pred_bud_gs_gdd_che2    <- make_sens_predictions2(m_sens_bud_gs_gdd_che)
pred_flower_gs_gdd_che2 <- make_sens_predictions2(m_sens_flower_gs_gdd_che)
pred_fruit_gs_gdd_che2  <- make_sens_predictions2(m_sens_fruit_gs_gdd_che)
pred_seed_gs_gdd_che2   <- make_sens_predictions2(m_sens_seed_gs_gdd_che)


# combine all predictions -------------------------------------------------
plot_predictions_gs_gdd <- bind_rows(
  pred_bud_gs_gdd_nor    |> mutate(stage = "Budding", region = "Norway"),
  pred_flower_gs_gdd_nor |> mutate(stage = "Flowering", region = "Norway"),
  pred_fruit_gs_gdd_nor  |> mutate(stage = "Fruiting", region = "Norway"),
  pred_seed_gs_gdd_nor   |> mutate(stage = "Seeds", region = "Norway"),
  pred_bud_gs_gdd_che    |> mutate(stage = "Budding", region = "Switzerland"),
  pred_flower_gs_gdd_che |> mutate(stage = "Flowering", region = "Switzerland"),
  pred_fruit_gs_gdd_che  |> mutate(stage = "Fruiting", region = "Switzerland"),
  pred_seed_gs_gdd_che   |> mutate(stage = "Seeds", region = "Switzerland")
)
plot_predictions_gs_gdd

# or with the ggpredict function
plot_predictions_gs_gdd2 <- bind_rows(
  pred_bud_gs_gdd_nor2    |> mutate(stage = "Budding", region = "Norway"),
  pred_flower_gs_gdd_nor2 |> mutate(stage = "Flowering", region = "Norway"),
  pred_fruit_gs_gdd_nor2  |> mutate(stage = "Fruiting", region = "Norway"),
  pred_seed_gs_gdd_nor2   |> mutate(stage = "Seeds", region = "Norway"),
  pred_bud_gs_gdd_che2    |> mutate(stage = "Budding", region = "Switzerland"),
  pred_flower_gs_gdd_che2 |> mutate(stage = "Flowering", region = "Switzerland"),
  pred_fruit_gs_gdd_che2  |> mutate(stage = "Fruiting", region = "Switzerland"),
  pred_seed_gs_gdd_che2   |> mutate(stage = "Seeds", region = "Switzerland")
)
plot_predictions_gs_gdd2

# rename column names 
plot_predictions_gs_gdd2 <- plot_predictions_gs_gdd2 |> 
  rename(treat_competition = x,
         temp_sens = predicted,
         plo = conf.low,
         phi = conf.high)



