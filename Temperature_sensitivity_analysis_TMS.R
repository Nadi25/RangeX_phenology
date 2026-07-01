

# 02_temp_sens TMS4 ------------------------------------------------------------

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

theme_set(theme_bw(base_size = 20))

# source script with functions --------------------------------------------
source("Functions_temperature_sensitivity.R")



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




# source clean phenology data -----------------------------------------------
source("Data_preparation_phenology_NOR_CHE_combined.R")

# use this data set
names(phenology2)




# calculate mean onset per species and individual for all stages ----------------
onset_bud    <- get_mean_onset(phenology2, "No_Buds", "jday")
onset_flower <- get_mean_onset(phenology2, "No_FloOpen", "jday")
onset_fruit  <- get_mean_onset(phenology2, "No_FloWithrd", "jday")
onset_seed   <- get_mean_onset(phenology2, "No_Seeds", "jday")






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



# Filter correct onset dataset low vs hi --------------------------------------------------
# per stage
# only ambi

# NOR
d_bud_lh_nor <- filter_data_ambi(onset_all_gs, "Norway", "Budding", "ambi")
d_flower_lh_nor <- filter_data_ambi(onset_all_gs, "Norway", "Flowering", "ambi")
d_fruit_lh_nor <- filter_data_ambi(onset_all_gs, "Norway", "Fruiting", "ambi")
d_seed_lh_nor <- filter_data_ambi(onset_all_gs, "Norway", "Seeds", "ambi")


# CHE
d_bud_lh_che <- filter_data_ambi(onset_all_gs, "Switzerland", "Budding", "ambi")
d_flower_lh_che <- filter_data_ambi(onset_all_gs, "Switzerland", "Flowering", "ambi")
d_fruit_lh_che <- filter_data_ambi(onset_all_gs, "Switzerland", "Fruiting", "ambi")
d_seed_lh_che <- filter_data_ambi(onset_all_gs, "Switzerland", "Seeds", "ambi")



# sensitivity models ------------------------------------------------------------------

# NOR ---------------------------------------------------------------------
# fit the models per stage for Norway
m_sens_bud_gs_lh_nor    <- fit_model_sens(d_bud_lh_nor)
m_sens_flower_gs_lh_nor <- fit_model_sens(d_flower_lh_nor)
m_sens_fruit_gs_lh_nor  <- fit_model_sens(d_fruit_lh_nor)
m_sens_seed_gs_lh_nor   <- fit_model_sens(d_seed_lh_nor)


summary(m_sens_bud_gs_lh_nor)
summary(m_sens_flower_gs_lh_nor)
summary(m_sens_fruit_gs_lh_nor)
summary(m_sens_seed_gs_lh_nor)

summary(m_sens_bud_gs_lh_nor)$coefficients["Tmean", ]

anova(m_sens_bud_gs_lh_nor)


# CHE ---------------------------------------------------------------------
# fit the models per stage for Switzerland
m_sens_bud_gs_lh_che    <- fit_model_sens(d_bud_lh_che)
m_sens_flower_gs_lh_che <- fit_model_sens(d_flower_lh_che)
m_sens_fruit_gs_lh_che  <- fit_model_sens(d_fruit_lh_che)
m_sens_seed_gs_lh_che   <- fit_model_sens(d_seed_lh_che)


summary(m_sens_bud_gs_lh_che)
summary(m_sens_flower_gs_lh_che)
summary(m_sens_fruit_gs_lh_che)
summary(m_sens_seed_gs_lh_che)

summary(m_sens_bud_gs_lh_che)$coefficients["Tmean", ]



# filter temp data --------------------------------------------------------
# this is to get the newdataframe for the predictions
temp_lh_nor <- temperature_mean_gs |>
  filter(
    region == "Norway",
    treat_warming == "ambi"
  )
temp_lh_nor

temp_lh_che <- temperature_mean_gs |>
  filter(
    region == "Switzerland",
    treat_warming == "ambi"
  )
temp_lh_che


# predict sensitivity for each stage with function ------------------------
# NOR ---------------------------------------------------------------------
#
pred_bud_lh_nor <- make_sens_predictions(temp_lh_nor, m_sens_bud_gs_lh_nor)
pred_flower_lh_nor <- make_sens_predictions(temp_lh_nor, m_sens_flower_gs_lh_nor)
pred_fruit_lh_nor <- make_sens_predictions(temp_lh_nor, m_sens_fruit_gs_lh_nor)
pred_seed_lh_nor <- make_sens_predictions(temp_lh_nor, m_sens_seed_gs_lh_nor)


# CHE ---------------------------------------------------------------------
#
pred_bud_lh_che <- make_sens_predictions(temp_lh_che, m_sens_bud_gs_lh_che)
pred_flower_lh_che <- make_sens_predictions(temp_lh_che, m_sens_flower_gs_lh_che)
pred_fruit_lh_che <- make_sens_predictions(temp_lh_che, m_sens_fruit_gs_lh_che)
pred_seed_lh_che <- make_sens_predictions(temp_lh_che, m_sens_seed_gs_lh_che)



# get the actual temperature sensitivity from coefficients ----------------
# NOR ---------------------------------------------------------------------

ts_bud_lh_nor <- get_temp_sens_coef(m_sens_bud_gs_lh_nor)
ts_flower_lh_nor <- get_temp_sens_coef(m_sens_flower_gs_lh_nor)
ts_fruit_lh_nor <- get_temp_sens_coef(m_sens_fruit_gs_lh_nor)
ts_seed_lh_nor <- get_temp_sens_coef(m_sens_seed_gs_lh_nor)


# CHE ---------------------------------------------------------------------

ts_bud_lh_che <- get_temp_sens_coef(m_sens_bud_gs_lh_che)
ts_flower_lh_che <- get_temp_sens_coef(m_sens_flower_gs_lh_che)
ts_fruit_lh_che <- get_temp_sens_coef(m_sens_fruit_gs_lh_che)
ts_seed_lh_che <- get_temp_sens_coef(m_sens_seed_gs_lh_che)



# combine sens from all stages -------------------------------------------
sens_all <- bind_rows(
  ts_bud_lh_nor    |> mutate(stage = "Budding", region = "Norway"),
  ts_flower_lh_nor |> mutate(stage = "Flowering", region = "Norway"),
  ts_fruit_lh_nor  |> mutate(stage = "Fruiting", region = "Norway"),
  ts_seed_lh_nor   |> mutate(stage = "Seeds", region = "Norway"),
  ts_bud_lh_che    |> mutate(stage = "Budding", region = "Switzerland"),
  ts_flower_lh_che |> mutate(stage = "Flowering", region = "Switzerland"),
  ts_fruit_lh_che  |> mutate(stage = "Fruiting", region = "Switzerland"),
  ts_seed_lh_che   |> mutate(stage = "Seeds", region = "Switzerland")
)
sens_all



# Plot sensitivity --------------------------------------------------------
ts <- ggplot(sens_all, aes(x = treat_competition, y = Tmean.trend, color = treat_competition)) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Temperature sensitivity (days / °C)",
       x = "Biotic interactions",
       title = "Temperature sensitivity low vs high ambient")+
  facet_grid(region ~ stage)+
  theme(legend.position = "none")+
  scale_color_manual(values = c("with" = "#528B8B", "without" = "#CD950C")) 
ts

# ggsave(filename = "Output/Sensitivity/Temperature_sensitivity_TMS_low_high_NOR_CHE.png", 
#       plot = ts,
#        width = 15, height = 10, units = "in")



# Pot raw data ------------------------------------------------------------
# bud
ggplot(d_bud_lh_nor,
       aes(Tmean, onset,
           colour = treat_competition)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm") +
  facet_wrap(~ species)



raw_sens_all <- bind_rows(
  d_bud_lh_nor,
  d_flower_lh_nor,
  d_fruit_lh_nor,
  d_seed_lh_nor,
  d_bud_lh_che,
  d_flower_lh_che,
  d_fruit_lh_che,
  d_seed_lh_che
)


ggplot(raw_sens_all,
       aes(Tmean, onset,
           colour = treat_competition)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_grid(region ~ stage)



ggplot(raw_sens_all,
       aes(Tmean, onset,
           colour = treat_competition,
           shape = region)) +
  geom_point(alpha = 0.7) +
  facet_grid( ~ stage)



ggplot(raw_sens_all,
       aes(Tmean, onset,
           colour = species)) +
  geom_point(alpha = .25) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_grid(region ~ stage)








# Filter correct onset dataset hi ambi vs warm --------------------------------------------------
# per stage
# only high

# NOR
d_bud_aw_nor <- filter_data_hi(onset_all_gs, "Norway", "Budding", "hi")
d_flower_aw_nor <- filter_data_hi(onset_all_gs, "Norway", "Flowering", "hi")
d_fruit_aw_nor <- filter_data_hi(onset_all_gs, "Norway", "Fruiting", "hi")
d_seed_aw_nor <- filter_data_hi(onset_all_gs, "Norway", "Seeds", "hi")


# CHE
d_bud_aw_che <- filter_data_hi(onset_all_gs, "Switzerland", "Budding", "hi")
d_flower_aw_che <- filter_data_hi(onset_all_gs, "Switzerland", "Flowering", "hi")
d_fruit_aw_che <- filter_data_hi(onset_all_gs, "Switzerland", "Fruiting", "hi")
d_seed_aw_che <- filter_data_hi(onset_all_gs, "Switzerland", "Seeds", "hi")



# sensitivity models ambi warm ------------------------------------------------------------------

# NOR ---------------------------------------------------------------------
# fit the models per stage for Norway
m_sens_bud_gs_aw_nor    <- fit_model_sens(d_bud_aw_nor)
m_sens_flower_gs_aw_nor <- fit_model_sens(d_flower_aw_nor)
m_sens_fruit_gs_aw_nor  <- fit_model_sens(d_fruit_aw_nor)
m_sens_seed_gs_aw_nor   <- fit_model_sens(d_seed_aw_nor)


summary(m_sens_bud_gs_aw_nor)
summary(m_sens_flower_gs_aw_nor)
summary(m_sens_fruit_gs_aw_nor)
summary(m_sens_seed_gs_aw_nor)

summary(m_sens_bud_gs_aw_nor)$coefficients["Tmean", ]

anova(m_sens_bud_gs_aw_nor)


# CHE ---------------------------------------------------------------------
# fit the models per stage for Switzerland
m_sens_bud_gs_aw_che    <- fit_model_sens(d_bud_aw_che)
m_sens_flower_gs_aw_che <- fit_model_sens(d_flower_aw_che)
m_sens_fruit_gs_aw_che  <- fit_model_sens(d_fruit_aw_che)
m_sens_seed_gs_aw_che   <- fit_model_sens(d_seed_aw_che)


summary(m_sens_bud_gs_aw_che)
summary(m_sens_flower_gs_aw_che)
summary(m_sens_fruit_gs_aw_che)
summary(m_sens_seed_gs_aw_che)

summary(m_sens_bud_gs_aw_che)$coefficients["Tmean", ]



# filter temp data --------------------------------------------------------
# this is to get the newdataframe for the predictions
temp_aw_nor <- temperature_mean_gs |>
  filter(
    region == "Norway",
    site == "hi"
  )
temp_aw_nor

temp_aw_che <- temperature_mean_gs |>
  filter(
    region == "Switzerland",
    site == "hi"
  )
temp_aw_che


# predict sensitivity for each stage with function ------------------------
# NOR ---------------------------------------------------------------------
#
pred_bud_aw_nor <- make_sens_predictions(temp_aw_nor, m_sens_bud_gs_aw_nor)
pred_flower_aw_nor <- make_sens_predictions(temp_aw_nor, m_sens_flower_gs_aw_nor)
pred_fruit_aw_nor <- make_sens_predictions(temp_aw_nor, m_sens_fruit_gs_aw_nor)
pred_seed_aw_nor <- make_sens_predictions(temp_aw_nor, m_sens_seed_gs_aw_nor)


# CHE ---------------------------------------------------------------------
#
pred_bud_aw_che <- make_sens_predictions(temp_aw_che, m_sens_bud_gs_aw_che)
pred_flower_aw_che <- make_sens_predictions(temp_aw_che, m_sens_flower_gs_aw_che)
pred_fruit_aw_che <- make_sens_predictions(temp_aw_che, m_sens_fruit_gs_aw_che)
pred_seed_aw_che <- make_sens_predictions(temp_aw_che, m_sens_seed_gs_aw_che)



# get the actual temperature sensitivity from coefficients ----------------
# NOR ---------------------------------------------------------------------

ts_bud_aw_nor <- get_temp_sens_coef(m_sens_bud_gs_aw_nor)
ts_flower_aw_nor <- get_temp_sens_coef(m_sens_flower_gs_aw_nor)
ts_fruit_aw_nor <- get_temp_sens_coef(m_sens_fruit_gs_aw_nor)
ts_seed_aw_nor <- get_temp_sens_coef(m_sens_seed_gs_aw_nor)


# CHE ---------------------------------------------------------------------

ts_bud_aw_che <- get_temp_sens_coef(m_sens_bud_gs_aw_che)
ts_flower_aw_che <- get_temp_sens_coef(m_sens_flower_gs_aw_che)
ts_fruit_aw_che <- get_temp_sens_coef(m_sens_fruit_gs_aw_che)
ts_seed_aw_che <- get_temp_sens_coef(m_sens_seed_gs_aw_che)



# combine sens from all stages -------------------------------------------
sens_all_aw <- bind_rows(
  ts_bud_aw_nor    |> mutate(stage = "Budding", region = "Norway"),
  ts_flower_aw_nor |> mutate(stage = "Flowering", region = "Norway"),
  ts_fruit_aw_nor  |> mutate(stage = "Fruiting", region = "Norway"),
  ts_seed_aw_nor   |> mutate(stage = "Seeds", region = "Norway"),
  ts_bud_aw_che    |> mutate(stage = "Budding", region = "Switzerland"),
  ts_flower_aw_che |> mutate(stage = "Flowering", region = "Switzerland"),
  ts_fruit_aw_che  |> mutate(stage = "Fruiting", region = "Switzerland"),
  ts_seed_aw_che   |> mutate(stage = "Seeds", region = "Switzerland")
)
sens_all_aw



# Plot sensitivity --------------------------------------------------------
ts_aw <- ggplot(sens_all_aw, aes(x = treat_competition, y = Tmean.trend, color = treat_competition)) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Temperature sensitivity (days / °C)",
       x = "Biotic interactions",
       title = "Temperature sensitivity high ambient vs warmed")+
  facet_grid(region ~ stage)+
  theme(legend.position = "none")+
  scale_color_manual(values = c("with" = "#528B8B", "without" = "#CD950C")) 
ts_aw

# ggsave(filename = "Output/Sensitivity/Temperature_sensitivity_TMS_hi_ambient_warmed_NOR_CHE.png", 
#       plot = ts_aw,
#        width = 15, height = 10, units = "in")




theme_set(theme_bw())
# Pot raw data ------------------------------------------------------------
# bud
ggplot(d_bud_aw_nor,
       aes(Tmean, onset,
           colour = treat_competition)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm") +
  facet_wrap(~ species)



raw_sens_all_aw <- bind_rows(
  d_bud_aw_nor,
  d_flower_aw_nor,
  d_fruit_aw_nor,
  d_seed_aw_nor,
  d_bud_aw_che,
  d_flower_aw_che,
  d_fruit_aw_che,
  d_seed_aw_che
)


ggplot(raw_sens_all_aw,
       aes(Tmean, onset,
           colour = treat_competition)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_grid(region ~ stage)


ggplot(raw_sens_all,
       aes(Tmean, onset,
           colour = species)) +
  geom_point(alpha = .25) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_grid(region ~ stage)


raw_summary <- raw_sens_all_aw |>
  group_by(region, stage, treatment_site_temp, treat_competition, Tmean) |>
  summarise(
    mean_onset = mean(onset),
    se = sd(onset) / sqrt(n()),
    .groups = "drop"
  )


ggplot(raw_summary,
       aes(x = Tmean,
           y = mean_onset,
           colour = treat_competition,
           shape = treatment_site_temp)) +
  
  geom_point(size = 3) +
  
  geom_errorbar(aes(ymin = mean_onset - 1.96 * se,
                    ymax = mean_onset + 1.96 * se),
                width = 0.05) +
  
  geom_line(aes(group = treat_competition)) +
  
  facet_grid(region ~ stage) +
  
  labs(
    x = "Mean growing-season temperature (°C)",
    y = "Mean onset (DOY)"
  )




# GDD ---------------------------------------------------------------------

source("GDD_Transplantation_onset_predictions_NOR_CHE.R")


# calculate first onset per species and plot for all stages ----------------
onset_bud_gdd    <- get_onset(phenology_gdd_nor_che, "No_Buds", "GDD_cum")
onset_flower_gdd <- get_onset(phenology_gdd_nor_che, "No_FloOpen", "GDD_cum")
onset_fruit_gdd  <- get_onset(phenology_gdd_nor_che, "No_FloWithrd", "GDD_cum")
onset_seed_gdd   <- get_onset(phenology_gdd_nor_che, "No_Seeds", "GDD_cum")




# Join phenology with temperature data ------------------------------------
onset_bud_temp_gdd <- onset_bud_gdd |>
  left_join(temperature_mean_gs,
            by = c("region", "treatment_site_temp", "treat_competition"))

onset_flower_temp_gdd <- onset_flower_gdd |>
  left_join(temperature_mean_gs,
            by = c("region", "treatment_site_temp", "treat_competition"))

onset_fruit_temp_gdd <- onset_fruit_gdd |>
  left_join(temperature_mean_gs,
            by = c("region", "treatment_site_temp", "treat_competition"))

onset_seed_temp_gdd <- onset_seed_gdd |>
  left_join(temperature_mean_gs,
            by = c("region", "treatment_site_temp", "treat_competition"))


# combine from all stages -------------------------------------------
onset_all_gs_gdd <- bind_rows(
  onset_bud_temp_gdd   |> mutate(stage = "Budding"),
  onset_flower_temp_gdd|> mutate(stage = "Flowering"),
  onset_fruit_temp_gdd |> mutate(stage = "Fruiting"),
  onset_seed_temp_gdd |> mutate(stage = "Seeds")
)
onset_all_gs_gdd


# Filter correct onset dataset low vs hi --------------------------------------------------
# per stage
# only ambi

# NOR
d_bud_gdd_lh_nor <- filter_data_ambi(onset_all_gs_gdd, "Norway", "Budding", "ambi")
d_flower_gdd_lh_nor <- filter_data_ambi(onset_all_gs_gdd, "Norway", "Flowering", "ambi")
d_fruit_gdd_lh_nor <- filter_data_ambi(onset_all_gs_gdd, "Norway", "Fruiting", "ambi")
d_seed_gdd_lh_nor <- filter_data_ambi(onset_all_gs_gdd, "Norway", "Seeds", "ambi")


# CHE
d_bud_gdd_lh_che <- filter_data_ambi(onset_all_gs_gdd, "Switzerland", "Budding", "ambi")
d_flower_gdd_lh_che <- filter_data_ambi(onset_all_gs_gdd, "Switzerland", "Flowering", "ambi")
d_fruit_gdd_lh_che <- filter_data_ambi(onset_all_gs_gdd, "Switzerland", "Fruiting", "ambi")
d_seed_gdd_lh_che <- filter_data_ambi(onset_all_gs_gdd, "Switzerland", "Seeds", "ambi")




# sensitivity models ------------------------------------------------------------------

# NOR ---------------------------------------------------------------------
# fit the models per stage for Norway
m_sens_bud_gs_gdd_lh_nor    <- fit_model_sens(d_bud_gdd_lh_nor)
m_sens_flower_gs_gdd_lh_nor <- fit_model_sens(d_flower_gdd_lh_nor)
m_sens_fruit_gs_gdd_lh_nor  <- fit_model_sens(d_fruit_gdd_lh_nor)
m_sens_seed_gs_gdd_lh_nor   <- fit_model_sens(d_seed_gdd_lh_nor)


summary(m_sens_bud_gs_gdd_lh_nor)
summary(m_sens_flower_gs_gdd_lh_nor)
summary(m_sens_fruit_gs_gdd_lh_nor)
summary(m_sens_seed_gs_gdd_lh_nor)

summary(m_sens_bud_gs_gdd_lh_nor)$coefficients["Tmean", ]

anova(m_sens_bud_gs_gdd_lh_nor)


# CHE ---------------------------------------------------------------------
# fit the models per stage for Switzerland
m_sens_bud_gs_gdd_lh_che    <- fit_model_sens(d_bud_gdd_lh_che)
m_sens_flower_gs_gdd_lh_che <- fit_model_sens(d_flower_gdd_lh_che)
m_sens_fruit_gs_gdd_lh_che  <- fit_model_sens(d_fruit_gdd_lh_che)
m_sens_seed_gs_gdd_lh_che   <- fit_model_sens(d_seed_gdd_lh_che)


summary(m_sens_bud_gs_gdd_lh_che)
summary(m_sens_flower_gs_gdd_lh_che)
summary(m_sens_fruit_gs_gdd_lh_che)
summary(m_sens_seed_gs_gdd_lh_che)

summary(m_sens_bud_gs_gdd_lh_che)$coefficients["Tmean", ]



# filter temp data --------------------------------------------------------
# this is to get the newdataframe for the predictions
temp_lh_nor <- temperature_mean_gs |>
  filter(
    region == "Norway",
    treat_warming == "ambi"
  )
temp_lh_nor

temp_lh_che <- temperature_mean_gs |>
  filter(
    region == "Switzerland",
    treat_warming == "ambi"
  )
temp_lh_che


# predict sensitivity for each stage with function ------------------------
# NOR ---------------------------------------------------------------------
#
pred_bud_gdd_lh_nor <- make_sens_predictions(temp_lh_nor, m_sens_bud_gs_gdd_lh_nor)
pred_flower_gdd_lh_nor <- make_sens_predictions(temp_lh_nor, m_sens_flower_gs_gdd_lh_nor)
pred_fruit_gdd_lh_nor <- make_sens_predictions(temp_lh_nor, m_sens_fruit_gs_gdd_lh_nor)
pred_seed_gdd_lh_nor <- make_sens_predictions(temp_lh_nor, m_sens_seed_gs_gdd_lh_nor)


# CHE ---------------------------------------------------------------------
#
pred_bud_gdd_lh_che <- make_sens_predictions(temp_lh_che, m_sens_bud_gs_gdd_lh_che)
pred_flower_gdd_lh_che <- make_sens_predictions(temp_lh_che, m_sens_flower_gs_gdd_lh_che)
pred_fruit_gdd_lh_che <- make_sens_predictions(temp_lh_che, m_sens_fruit_gs_gdd_lh_che)
pred_seed_gdd_lh_che <- make_sens_predictions(temp_lh_che, m_sens_seed_gs_gdd_lh_che)



# get the actual temperature sensitivity from coefficients ----------------
# NOR ---------------------------------------------------------------------

ts_bud_gdd_lh_nor <- get_temp_sens_coef(m_sens_bud_gs_gdd_lh_nor)
ts_flower_gdd_lh_nor <- get_temp_sens_coef(m_sens_flower_gs_gdd_lh_nor)
ts_fruit_gdd_lh_nor <- get_temp_sens_coef(m_sens_fruit_gs_gdd_lh_nor)
ts_seed_gdd_lh_nor <- get_temp_sens_coef(m_sens_seed_gs_gdd_lh_nor)


# CHE ---------------------------------------------------------------------

ts_bud_gdd_lh_che <- get_temp_sens_coef(m_sens_bud_gs_gdd_lh_che)
ts_flower_gdd_lh_che <- get_temp_sens_coef(m_sens_flower_gs_gdd_lh_che)
ts_fruit_gdd_lh_che <- get_temp_sens_coef(m_sens_fruit_gs_gdd_lh_che)
ts_seed_gdd_lh_che <- get_temp_sens_coef(m_sens_seed_gs_gdd_lh_che)



# combine sens from all stages -------------------------------------------
sens_all_gdd <- bind_rows(
  ts_bud_gdd_lh_nor    |> mutate(stage = "Budding", region = "Norway"),
  ts_flower_gdd_lh_nor |> mutate(stage = "Flowering", region = "Norway"),
  ts_fruit_gdd_lh_nor  |> mutate(stage = "Fruiting", region = "Norway"),
  ts_seed_gdd_lh_nor   |> mutate(stage = "Seeds", region = "Norway"),
  ts_bud_gdd_lh_che    |> mutate(stage = "Budding", region = "Switzerland"),
  ts_flower_gdd_lh_che |> mutate(stage = "Flowering", region = "Switzerland"),
  ts_fruit_gdd_lh_che  |> mutate(stage = "Fruiting", region = "Switzerland"),
  ts_seed_gdd_lh_che   |> mutate(stage = "Seeds", region = "Switzerland")
)
sens_all_gdd



# Plot sensitivity --------------------------------------------------------
ts_gdd <- ggplot(sens_all_gdd, aes(x = treat_competition, y = Tmean.trend, color = treat_competition)) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Temperature sensitivity (GDD2 / °C)",
       x = "Biotic interactions",
       title = "Temperature sensitivity low vs high ambient")+
  facet_grid(region ~ stage)+
  theme(legend.position = "none")+
  scale_color_manual(values = c("with" = "#528B8B", "without" = "#CD950C")) 
ts_gdd

# ggsave(filename = "Output/Sensitivity/Temperature_sensitivity_TMS_GDD_low_high_NOR_CHE.png", 
#       plot = ts_gdd,
#        width = 15, height = 10, units = "in")




# GDD hi ambi vs warm -----------------------------------------------------


# Filter correct onset dataset hi ambi vs warm --------------------------------------------------
# per stage
# only hi

# NOR
d_bud_gdd_aw_nor <- filter_data_hi(onset_all_gs_gdd, "Norway", "Budding", "hi")
d_flower_gdd_aw_nor <- filter_data_hi(onset_all_gs_gdd, "Norway", "Flowering", "hi")
d_fruit_gdd_aw_nor <- filter_data_hi(onset_all_gs_gdd, "Norway", "Fruiting", "hi")
d_seed_gdd_aw_nor <- filter_data_hi(onset_all_gs_gdd, "Norway", "Seeds", "hi")


# CHE
d_bud_gdd_aw_che <- filter_data_hi(onset_all_gs_gdd, "Switzerland", "Budding", "hi")
d_flower_gdd_aw_che <- filter_data_hi(onset_all_gs_gdd, "Switzerland", "Flowering", "hi")
d_fruit_gdd_aw_che <- filter_data_hi(onset_all_gs_gdd, "Switzerland", "Fruiting", "hi")
d_seed_gdd_aw_che <- filter_data_hi(onset_all_gs_gdd, "Switzerland", "Seeds", "hi")




# sensitivity models ------------------------------------------------------------------

# NOR ---------------------------------------------------------------------
# fit the models per stage for Norway
m_sens_bud_gs_gdd_aw_nor    <- fit_model_sens(d_bud_gdd_aw_nor)
m_sens_flower_gs_gdd_aw_nor <- fit_model_sens(d_flower_gdd_aw_nor)
m_sens_fruit_gs_gdd_aw_nor  <- fit_model_sens(d_fruit_gdd_aw_nor)
m_sens_seed_gs_gdd_aw_nor   <- fit_model_sens(d_seed_gdd_aw_nor)


summary(m_sens_bud_gs_gdd_aw_nor)
summary(m_sens_flower_gs_gdd_aw_nor)
summary(m_sens_fruit_gs_gdd_aw_nor)
summary(m_sens_seed_gs_gdd_aw_nor)

summary(m_sens_bud_gs_gdd_aw_nor)$coefficients["Tmean", ]

anova(m_sens_bud_gs_gdd_aw_nor)


# CHE ---------------------------------------------------------------------
# fit the models per stage for Switzerland
m_sens_bud_gs_gdd_aw_che    <- fit_model_sens(d_bud_gdd_aw_che)
m_sens_flower_gs_gdd_aw_che <- fit_model_sens(d_flower_gdd_aw_che)
m_sens_fruit_gs_gdd_aw_che  <- fit_model_sens(d_fruit_gdd_aw_che)
m_sens_seed_gs_gdd_aw_che   <- fit_model_sens(d_seed_gdd_aw_che)


summary(m_sens_bud_gs_gdd_aw_che)
summary(m_sens_flower_gs_gdd_aw_che)
summary(m_sens_fruit_gs_gdd_aw_che)
summary(m_sens_seed_gs_gdd_aw_che)

summary(m_sens_bud_gs_gdd_aw_che)$coefficients["Tmean", ]


# filter temp data --------------------------------------------------------
# this is to get the newdataframe for the predictions
# use from above
temp_aw_nor 

temp_aw_che



# predict sensitivity for each stage with function ------------------------
# NOR ---------------------------------------------------------------------
#
pred_bud_gdd_aw_nor <- make_sens_predictions(temp_aw_nor, m_sens_bud_gs_gdd_aw_nor)
pred_flower_gdd_aw_nor <- make_sens_predictions(temp_aw_nor, m_sens_flower_gs_gdd_aw_nor)
pred_fruit_gdd_aw_nor <- make_sens_predictions(temp_aw_nor, m_sens_fruit_gs_gdd_aw_nor)
pred_seed_gdd_aw_nor <- make_sens_predictions(temp_aw_nor, m_sens_seed_gs_gdd_aw_nor)


# CHE ---------------------------------------------------------------------
#
pred_bud_gdd_aw_che <- make_sens_predictions(temp_aw_che, m_sens_bud_gs_gdd_aw_che)
pred_flower_gdd_aw_che <- make_sens_predictions(temp_aw_che, m_sens_flower_gs_gdd_aw_che)
pred_fruit_gdd_aw_che <- make_sens_predictions(temp_aw_che, m_sens_fruit_gs_gdd_aw_che)
pred_seed_gdd_aw_che <- make_sens_predictions(temp_aw_che, m_sens_seed_gs_gdd_aw_che)



# get the actual temperature sensitivity from coefficients ----------------
# NOR ---------------------------------------------------------------------

ts_bud_gdd_aw_nor <- get_temp_sens_coef(m_sens_bud_gs_gdd_aw_nor)
ts_flower_gdd_aw_nor <- get_temp_sens_coef(m_sens_flower_gs_gdd_aw_nor)
ts_fruit_gdd_aw_nor <- get_temp_sens_coef(m_sens_fruit_gs_gdd_aw_nor)
ts_seed_gdd_aw_nor <- get_temp_sens_coef(m_sens_seed_gs_gdd_aw_nor)


# CHE ---------------------------------------------------------------------

ts_bud_gdd_aw_che <- get_temp_sens_coef(m_sens_bud_gs_gdd_aw_che)
ts_flower_gdd_aw_che <- get_temp_sens_coef(m_sens_flower_gs_gdd_aw_che)
ts_fruit_gdd_aw_che <- get_temp_sens_coef(m_sens_fruit_gs_gdd_aw_che)
ts_seed_gdd_aw_che <- get_temp_sens_coef(m_sens_seed_gs_gdd_aw_che)



# combine sens from all stages -------------------------------------------
sens_all_gdd_aw <- bind_rows(
  ts_bud_gdd_aw_nor    |> mutate(stage = "Budding", region = "Norway"),
  ts_flower_gdd_aw_nor |> mutate(stage = "Flowering", region = "Norway"),
  ts_fruit_gdd_aw_nor  |> mutate(stage = "Fruiting", region = "Norway"),
  ts_seed_gdd_aw_nor   |> mutate(stage = "Seeds", region = "Norway"),
  ts_bud_gdd_aw_che    |> mutate(stage = "Budding", region = "Switzerland"),
  ts_flower_gdd_aw_che |> mutate(stage = "Flowering", region = "Switzerland"),
  ts_fruit_gdd_aw_che  |> mutate(stage = "Fruiting", region = "Switzerland"),
  ts_seed_gdd_aw_che   |> mutate(stage = "Seeds", region = "Switzerland")
)
sens_all_gdd_aw



# Plot sensitivity --------------------------------------------------------
ts_gdd_aw <- ggplot(sens_all_gdd_aw, aes(x = treat_competition, y = Tmean.trend, color = treat_competition)) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Temperature sensitivity (GDD2 / °C)",
       x = "Biotic interactions",
       title = "Temperature sensitivity high ambient vs warmed")+
  facet_grid(region ~ stage)+
  theme(legend.position = "none")+
  scale_color_manual(values = c("with" = "#528B8B", "without" = "#CD950C")) 
ts_gdd_aw

# ggsave(filename = "Output/Sensitivity/Temperature_sensitivity_TMS_GDD_ambi_warm_NOR_CHE.png", 
#       plot = ts_gdd_aw,
#       width = 15, height = 10, units = "in")

