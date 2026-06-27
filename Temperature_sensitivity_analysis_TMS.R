

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



source("DOY_GDD_TMS4_NOR_CHE.R")

names(tms_final_nor_che)

# mean temperature growing season mean ------------------------------------------
# define growing season length
tms_final_nor_che_gs <- tms_final_nor_che |>
  filter((region == "Norway" & date >= as.Date("2023-04-01") & date <= as.Date("2023-09-30")) |
           (region == "Switzerland" & date >= as.Date("2022-04-01") & date <= as.Date("2022-09-30")))


# calculate one mean value per growing seson and treatment and plot
temperature_mean_gs <-  tms_final_nor_che_gs |>
  group_by(region, site, treat_warming, treatment_site_temp, treat_competition) |>
  summarise(Tmean = mean(temp_mean),
            .groups = "drop")
temperature_mean_gs

# with without
temperature_mean_gs <- temperature_mean_gs |> 
  mutate(treat_competition = recode(treat_competition,
                                    "bare" = "without",
                                    "vege" = "with"))



# # source script with final tms data ---------------------------------------
# # this one sources the NOR and CHE preparation scripts
# source("Data_preparation_TMS4_CHE.R")
# names(tomst_che_daily)
# 
# source("TMS4_Weatherstation_predictions_NOR.R")
# names(tms_joined2)
# 
# 
# tms_daily_che <- tomst_che_daily |> 
#   mutate(temp_mean_used = T3_mean) |> 
#   mutate(region = "Switzerland") |> 
#   select(region, date, site, treat_warming, treat_competition, added_focals, unique_plot_ID, temp_mean_used)
# 
# 
# tms_daily_nor <- tms_joined2 |> 
#   mutate(region = "Norway") |> 
#   select(region, date, site, treat_warming, treat_competition, added_focals, unique_plot_ID, temp_mean_used)
# 
# 
# # combine
# tms_daily_nor_che <- bind_rows(tms_daily_nor, tms_daily_che)
# 
# # add site_warming treatment ----------------------------------------------
# tms_daily_nor_che$treatment_site_temp <- paste(tms_daily_nor_che$site, tms_daily_nor_che$treat_warming, sep = "_")
# 
# 
# tms_daily_nor_che <- tms_daily_nor_che |>
#   mutate(treatment_site_temp= factor(treatment_site_temp,
#                                      levels = c("lo_ambi",
#                                                 "hi_ambi",
#                                                 "hi_warm")))
# # with without
# tms_daily_nor_che <- tms_daily_nor_che |> 
#   mutate(treat_competition = recode(treat_competition,
#                                     "bare" = "without",
#                                     "vege" = "with"))
# 
# 
# # mean temperature growing season mean ------------------------------------------
# # define growing season length
# tms_daily_nor_che_gs <- tms_daily_nor_che |>
#   filter((region == "Norway" & date >= as.Date("2023-04-01") & date <= as.Date("2023-09-30")) |
#       (region == "Switzerland" & date >= as.Date("2022-04-01") & date <= as.Date("2022-09-30")))
# 
# # calculate one mean value per growing seson and treatment and plot
# temperature_mean_gs <-  tms_daily_nor_che_gs |>
#   group_by(region, treatment_site_temp, treat_competition, unique_plot_ID) |>
#   summarise(Tmean = mean(temp_mean_used),
#             .groups = "drop")
# temperature_mean_gs
# 


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
            by = c("region", "treatment_site_temp", "treat_competition"))

onset_flower_temp <- onset_flower |>
  left_join(temperature_mean_gs,
            by = c("region", "treatment_site_temp", "treat_competition"))

onset_fruit_temp <- onset_fruit |>
  left_join(temperature_mean_gs,
            by = c("region", "treatment_site_temp", "treat_competition"))

onset_seed_temp <- onset_seed |>
  left_join(temperature_mean_gs,
            by = c("region", "treatment_site_temp", "treat_competition"))


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
    y = expression("Onset"),
    color = "Biotic interactions"
  )+
  facet_wrap(~ region)





# sensitivity models ------------------------------------------------------------------

# NOR ---------------------------------------------------------------------
# fit the models per stage for Norway
m_sens_bud_gs_lh_nor    <- fit_model_sens_lo_hi(onset_all_gs, "Norway", "Budding", "ambi")
m_sens_flower_gs_lh_nor <- fit_model_sens_lo_hi(onset_all_gs, "Norway", "Flowering", "ambi")
m_sens_fruit_gs_lh_nor  <- fit_model_sens_lo_hi(onset_all_gs, "Norway", "Fruiting", "ambi")
m_sens_seed_gs_lh_nor   <- fit_model_sens_lo_hi(onset_all_gs, "Norway", "Seeds", "ambi")


summary(m_sens_bud_gs_lh_nor)
summary(m_sens_flower_gs_lh_nor)
summary(m_sens_fruit_gs_lh_nor)
summary(m_sens_seed_gs_lh_nor)

summary(m_sens_bud_gs_lh_nor)$coefficients["Tmean", ]





# predict sensitivity for each stage with function ------------------------
# NOR ---------------------------------------------------------------------
#
md_bud_lh_nor <- onset_all_gs |>
  filter(
    region == "Norway",
    stage == "Budding",
    treat_warming == "ambi"
  )

md_temp_lh_nor <- temperature_mean_gs |>
  filter(
    region == "Norway",
    treat_warming == "ambi"
  )


pred_bud_gs_lh_nor    <- make_sens_predictions_lh(m_sens_bud_gs_lh_nor, md_bud_lh_nor, md_temp_lh_nor)








pred_flower_gs_nor <- make_sens_predictions(m_sens_flower_gs_nor)
pred_fruit_gs_nor  <- make_sens_predictions(m_sens_fruit_gs_nor)
pred_seed_gs_nor   <- make_sens_predictions(m_sens_seed_gs_nor)



#######################################

model_data <- onset_all_gs |>
  filter(
    region == "Norway",
    stage == "Budding",
    treat_warming == "ambi"
  )

m <- lmer(
  onset ~ treat_competition * Tmean +
    (1 | species) + (1 | block_ID),
  data = model_data
)

summary(m)



newdata <- md_temp_lh_nor |>
  select(treat_competition, Tmean)


pred <- predict(m, newdata = newdata, re.form = NA, se.fit = TRUE)


newdata$fit <- pred$fit
newdata$se  <- pred$se.fit

newdata <- newdata |>
  mutate(
    lower = fit - 1.96 * se,
    upper = fit + 1.96 * se
  )
newdata



coef(m)
# slope for "without"
b0 <- fixef(m)["Tmean"]

# interaction term
b_int <- fixef(m)["treat_competitionwithout:Tmean"]

# slope for "with"
b_with <- b0 + b_int


ggplot(newdata, aes(x = Tmean, y = fit, color = treat_competition)) +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = treat_competition),
              alpha = 0.2, color = NA) +
  scale_color_manual(values = c("with" = "#528B8B",
                                "without" = "#CD950C")) +
  scale_fill_manual(values = c("with" = "#528B8B",
                               "without" = "#CD950C")) +
  labs(x = "Mean temperature",
       y = "Predicted onset") 

ggplot(newdata, aes(x = Tmean, y = fit, color = treat_competition)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  scale_color_manual(values = c("with" = "#528B8B",
                                "without" = "#CD950C")) 


slopes <- data.frame(
  treat_competition = c("without", "with"),
  slope = c(
    fixef(m)["Tmean"],
    fixef(m)["Tmean"] + fixef(m)["treat_competitionwithout:Tmean"]
  )
)

ggplot(slopes, aes(x = treat_competition, y = slope, fill = treat_competition)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Temperature sensitivity (days / °C)")


library(emmeans)

emtrends(m, ~ treat_competition, var = "Tmean")

slopes <- as.data.frame(emtrends(m, ~ treat_competition, var = "Tmean"))

ggplot(slopes, aes(x = treat_competition, y = Tmean.trend, color = treat_competition)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Temperature sensitivity (days / °C)")

