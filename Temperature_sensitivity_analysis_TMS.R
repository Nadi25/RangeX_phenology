

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


# source script with final tms data ---------------------------------------
# this one sources the NOR and CHE preparation scripts
source("Data_preparation_TMS4_CHE.R")
names(tomst_che_daily)

source("TMS4_Weatherstation_predictions_NOR.R")
names(tms_joined2)


tms_daily_che <- tomst_che_daily |> 
  mutate(temp_mean_used = T3_mean) |> 
  mutate(region = "Switzerland") |> 
  select(region, date, site, treat_warming, treat_competition, added_focals, unique_plot_ID, temp_mean_used)


tms_daily_nor <- tms_joined2 |> 
  mutate(region = "Norway") |> 
  select(region, date, site, treat_warming, treat_competition, added_focals, unique_plot_ID, temp_mean_used)


# combine
tms_daily_nor_che <- bind_rows(tms_daily_nor, tms_daily_che)

# add site_warming treatment ----------------------------------------------
tms_daily_nor_che$treatment_site_temp <- paste(tms_daily_nor_che$site, tms_daily_nor_che$treat_warming, sep = "_")


tms_daily_nor_che <- tms_daily_nor_che |>
  mutate(treatment_site_temp= factor(treatment_site_temp,
                                     levels = c("lo_ambi",
                                                "hi_ambi",
                                                "hi_warm")))
# with without
tms_daily_nor_che <- tms_daily_nor_che |> 
  mutate(treat_competition = recode(treat_competition,
                                    "bare" = "without",
                                    "vege" = "with"))


# mean temperature growing season mean ------------------------------------------
# define growing season length
tms_daily_nor_che_gs <- tms_daily_nor_che |>
  filter((region == "Norway" & date >= as.Date("2023-04-01") & date <= as.Date("2023-09-30")) |
      (region == "Switzerland" & date >= as.Date("2022-04-01") & date <= as.Date("2022-09-30")))

# calculate one mean value per growing seson and treatment and plot
temperature_mean_gs <-  tms_daily_nor_che_gs |>
  group_by(region, treatment_site_temp, treat_competition, unique_plot_ID) |>
  summarise(Tmean = mean(temp_mean_used),
            .groups = "drop")
temperature_mean_gs



# source clean phenology data -----------------------------------------------
source("Data_preparation_phenology_NOR_CHE_combined.R")

# use this data set
names(phenology2)




# calculate mean onset per species and individual for all stages ----------------
onset_bud    <- get_mean_onset(phenology2, "No_Buds", "jday")
onset_flower <- get_mean_onset(phenology2, "No_FloOpen", "jday")
onset_fruit  <- get_mean_onset(phenology2, "No_FloWithrd", "jday")
onset_seed   <- get_mean_onset(phenology2, "No_Seeds", "jday")


# # Average across plot, species and treat -----------------------------------------
# onset_bud_mean_p    <- get_plot_onset(onset_bud)
# onset_flower_mean_p <- get_plot_onset(onset_flower)
# onset_fruit_mean_p  <- get_plot_onset(onset_fruit)
# onset_seed_mean_p   <- get_plot_onset(onset_seed)




# Join phenology with temperature data ------------------------------------
onset_bud_temp <- onset_bud |>
  left_join(temperature_mean_gs,
            by = c("region", "treatment_site_temp", "treat_competition", "unique_plot_ID"))

onset_flower_temp <- onset_flower |>
  left_join(temperature_mean_gs,
            by = c("region", "treatment_site_temp", "treat_competition", "unique_plot_ID"))

onset_fruit_temp <- onset_fruit |>
  left_join(temperature_mean_gs,
            by = c("region", "treatment_site_temp", "treat_competition", "unique_plot_ID"))

onset_seed_temp <- onset_seed |>
  left_join(temperature_mean_gs,
            by = c("region", "treatment_site_temp", "treat_competition", "unique_plot_ID"))


# combine from all stages -------------------------------------------
onset_all_gs <- bind_rows(
  onset_bud_temp   |> mutate(stage = "Budding"),
  onset_flower_temp|> mutate(stage = "Flowering"),
  onset_fruit_temp |> mutate(stage = "Fruiting"),
  onset_seed_temp  |> mutate(stage = "Seeds")
)
onset_all_gs


# quick control plot ------------------------------------------------------
# plot the raw sens data 
ggplot(onset_all_gs,
       aes(x = stage,
           y = onset,
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
m_sens_bud_gs_nor    <- fit_model_sens(onset_all_gs, "Norway", "Budding")
m_sens_flower_gs_nor <- fit_model_sens(onset_all_gs, "Norway", "Flowering")
m_sens_fruit_gs_nor  <- fit_model_sens(onset_all_gs, "Norway", "Fruiting")
m_sens_seed_gs_nor   <- fit_model_sens(onset_all_gs, "Norway", "Seeds")


summary(m_sens_bud_gs_nor)
summary(m_sens_flower_gs_nor)
summary(m_sens_fruit_gs_nor)
summary(m_sens_seed_gs_nor)

summary(m_sens_bud_gs_nor)$coefficients["Tmean", ]




# predict sensitivity for each stage with function ------------------------
# NOR ---------------------------------------------------------------------
#
pred_bud_gs_nor    <- make_sens_predictions(m_sens_bud_gs_nor, onset_bud_temp)
pred_flower_gs_nor <- make_sens_predictions(m_sens_flower_gs_nor)
pred_fruit_gs_nor  <- make_sens_predictions(m_sens_fruit_gs_nor)
pred_seed_gs_nor   <- make_sens_predictions(m_sens_seed_gs_nor)









