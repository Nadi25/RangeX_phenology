


# 02 Onset GDD ------------------------------------------------------------

# Effect of upslope transplantation on onsets GDD  -------------------------------------------------------------

# load library ------------------------------------------------------------
library(conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
library(tidyverse)
library(lubridate)
library(performance)
library(see)
library(emmeans)

# source script with functions --------------------------------------------
source("Functions_onset.R")


# load clean phenology data -----------------------------------------------
source("Data_preparation_phenology_NOR_CHE_combined.R")

theme_set(theme_bw())

# filter only nor ------------------------------------------------------
phenology_nor <- phenology2 |> 
  filter(region == "Norway")


unique(phenology_nor$date_measurement)
# start: 2023-05-12
# end: 2023-09-13


source("Data_preparation_climate_station_NOR.R")

# use climate_gdd instead of climate_gdd_pt because it is without filtering
climate_gdd

# combine gdd_cum and phenology_nor ------------------------------------
phenology_with_gdd_nor <- phenology_nor |> 
  left_join(climate_gdd |> 
              select(site, date_measurement, GDD_cum, Tavg),
            by = c("site", "date_measurement"))



# filter only che ------------------------------------------------------
phenology_che <- phenology2 |> 
  filter(region == "Switzerland")


unique(phenology_che$date_measurement)
# start: "2022-05-04"
# end: "2022-09-27"


source("Data_preparation_climate_station_CHE.R")

# use 
# climate_gdd_che
# uses tms data for beginning until May and then climate station
climate_gdd_che_comb

# combine gdd_cum and phenology_nor ------------------------------------
phenology_with_gdd_che <- phenology_che |> 
  left_join(climate_gdd_che_comb |> 
              select(site, date_measurement, GDD_cum, Tavg),
            by = c("site", "date_measurement"))




# combine nor and che -----------------------------------------------------
phenology_gdd_nor_che <- rbind(phenology_with_gdd_nor, phenology_with_gdd_che)



# filter only ambi both sites  -------------------------------------
# to compare low ambi with hi ambi = cooling effect
phenology_gdd_nor_che <- phenology_gdd_nor_che |> 
  filter(treat_warming == "ambi")





# calculate first onset per species and plot for all stages ----------------
onset_bud_gdd    <- get_onset(phenology_gdd_nor_che, "No_Buds", "GDD_cum")
onset_flower_gdd <- get_onset(phenology_gdd_nor_che, "No_FloOpen", "GDD_cum")
onset_fruit_gdd  <- get_onset(phenology_gdd_nor_che, "No_FloWithrd", "GDD_cum")
onset_seed_gdd   <- get_onset(phenology_gdd_nor_che, "No_Seeds", "GDD_cum")



# fit onset model ---------------------------------------------------------

# NOR ---------------------------------------------------------------------
# fit the models per stage for Norway
m_onset_bud_gdd_nor    <- fit_onset_model(onset_bud_gdd, "Norway")
m_onset_flower_gdd_nor <- fit_onset_model(onset_flower_gdd, "Norway")
m_onset_fruit_gdd_nor  <- fit_onset_model(onset_fruit_gdd, "Norway")
m_onset_seed_gdd_nor   <- fit_onset_model(onset_seed_gdd, "Norway")

# check model output
# bud
summary(m_onset_bud_gdd_nor)
anova(m_onset_bud_gdd_nor)
model_performance(m_onset_bud_gdd_nor)
#check_model(m_onset_bud_gdd_nor)

emmeans(m_onset_bud_gdd_nor,
        pairwise ~ treatment_site_temp * treat_competition)


# flower
summary(m_onset_flower_gdd_nor)
anova(m_onset_flower_gdd_nor)
model_performance(m_onset_flower_gdd_nor)
#check_model(m_onset_flower_gdd_nor)

emmeans(m_onset_flower_gdd_nor,
        pairwise ~ treatment_site_temp * treat_competition)

# fruit
summary(m_onset_fruit_gdd_nor)
anova(m_onset_fruit_gdd_nor)
model_performance(m_onset_fruit_gdd_nor)
#check_model(m_onset_fruit_gdd_nor)

emmeans(m_onset_fruit_gdd_nor,
        pairwise ~ treatment_site_temp * treat_competition)

# seed
summary(m_onset_seed_gdd_nor)
anova(m_onset_seed_gdd_nor)
model_performance(m_onset_seed_gdd_nor)
#check_model(m_onset_seed_gdd_nor)

emmeans(m_onset_seed_gdd_nor,
        pairwise ~ treatment_site_temp * treat_competition)






# CHE ---------------------------------------------------------------------
# fit the models per stage for Switzerland
m_onset_bud_gdd_che    <- fit_onset_model(onset_bud_gdd, "Switzerland")
m_onset_flower_gdd_che <- fit_onset_model(onset_flower_gdd, "Switzerland")
m_onset_fruit_gdd_che  <- fit_onset_model(onset_fruit_gdd, "Switzerland")
m_onset_seed_gdd_che   <- fit_onset_model(onset_seed_gdd, "Switzerland")

# check model output
# bud
summary(m_onset_bud_gdd_che)
anova(m_onset_bud_gdd_che)
model_performance(m_onset_bud_gdd_che)
#check_model(m_onset_bud_gdd_nor)

emmeans(m_onset_bud_gdd_che,
        pairwise ~ treatment_site_temp * treat_competition)


# flower
summary(m_onset_flower_gdd_che)
anova(m_onset_flower_gdd_che)
model_performance(m_onset_flower_gdd_che)
#check_model(m_onset_flower_gdd_che)

emmeans(m_onset_flower_gdd_che,
        pairwise ~ treatment_site_temp * treat_competition)

# fruit
summary(m_onset_fruit_gdd_che)
anova(m_onset_fruit_gdd_che)
model_performance(m_onset_fruit_gdd_che)
#check_model(m_onset_fruit_gdd_che)

emmeans(m_onset_fruit_gdd_che,
        pairwise ~ treatment_site_temp * treat_competition)

# seed
summary(m_onset_seed_gdd_che)
anova(m_onset_seed_gdd_che)
model_performance(m_onset_seed_gdd_che)
#check_model(m_onset_seed_gdd_che)

emmeans(m_onset_seed_gdd_che,
        pairwise ~ treatment_site_temp * treat_competition)




# make onset predictions --------------------------------------------------

# NOR ---------------------------------------------------------------------
# 
pred_onset_bud_gdd_nor    <- make_onset_predictions_gdd(m_onset_bud_gdd_nor)
pred_onset_flower_gdd_nor <- make_onset_predictions_gdd(m_onset_flower_gdd_nor)
pred_onset_fruit_gdd_nor  <- make_onset_predictions_gdd(m_onset_fruit_gdd_nor)
pred_onset_seed_gdd_nor   <- make_onset_predictions_gdd(m_onset_seed_gdd_nor)


# CHE ---------------------------------------------------------------------
# 
pred_onset_bud_gdd_che    <- make_onset_predictions_gdd(m_onset_bud_gdd_che)
pred_onset_flower_gdd_che <- make_onset_predictions_gdd(m_onset_flower_gdd_che)
pred_onset_fruit_gdd_che  <- make_onset_predictions_gdd(m_onset_fruit_gdd_che)
pred_onset_seed_gdd_che   <- make_onset_predictions_gdd(m_onset_seed_gdd_che)





# combine predictions into one dataframe ----------------------------------

# nor
plot_df_bud_gdd_nor  <- pred_onset_bud_gdd_nor   |> 
  mutate(stage = "Budding",
         region = "Norway")

plot_df_flower_gdd_nor  <- pred_onset_flower_gdd_nor   |> 
  mutate(stage = "Flowering",
         region = "Norway")

plot_df_fruit_gdd_nor  <- pred_onset_fruit_gdd_nor   |> 
  mutate(stage = "Fruiting",
         region = "Norway")

plot_df_seed_gdd_nor  <- pred_onset_seed_gdd_nor   |> 
  mutate(stage = "Seeds",
         region = "Norway")


# che
plot_df_bud_gdd_che  <- pred_onset_bud_gdd_che   |> 
  mutate(stage = "Budding",
         region = "Switzerland")

plot_df_flower_gdd_che  <- pred_onset_flower_gdd_che   |> 
  mutate(stage = "Flowering",
         region = "Switzerland")

plot_df_fruit_gdd_che  <- pred_onset_fruit_gdd_che   |> 
  mutate(stage = "Fruiting",
         region = "Switzerland")

plot_df_seed_gdd_che  <- pred_onset_seed_gdd_che   |> 
  mutate(stage = "Seeds",
         region = "Switzerland")


plot_df_all_gdd <- bind_rows(
  plot_df_bud_gdd_nor,
  plot_df_flower_gdd_nor,
  plot_df_fruit_gdd_nor,
  plot_df_seed_gdd_nor,
  plot_df_bud_gdd_che,
  plot_df_flower_gdd_che,
  plot_df_fruit_gdd_che,
  plot_df_seed_gdd_che
)
plot_df_all_gdd



# combine raw onset data into one data frame ----------------------------------
raw_bud_gdd  <- onset_bud_gdd   |> 
  mutate(stage = "Budding")

raw_flower_gdd  <- onset_flower_gdd  |> 
  mutate(stage = "Flowering")

raw_fruit_gdd  <- onset_fruit_gdd   |> 
  mutate(stage = "Fruiting")

raw_seed_gdd  <- onset_seed_gdd   |> 
  mutate(stage = "Seeds")


plot_df_raw_all_gdd <- bind_rows(
  raw_bud_gdd,
  raw_flower_gdd,
  raw_fruit_gdd,
  raw_seed_gdd
)
plot_df_raw_all_gdd












pd <- position_dodge(width = 0.6) 


b_f_fr_gdd <- ggplot(plot_df_all_gdd, aes(
  x = treatment_site_temp,
  y = fit,
  color = treat_competition,
  shape = treatment_site_temp
)) +
  
  # raw:
  geom_violin(
    data = plot_df_raw_all_gdd,
    aes(
      x = treatment_site_temp,
      y = onset,
      fill = treat_competition
      #group = interaction(treatment_site_temp, treat_competition)
    ),
    position = pd,
    alpha = 0.25,
    color = NA,
    trim = FALSE,
    adjust = 1.5 # makes it smoother so it looks less bubbly
  ) +
  
  scale_fill_manual(values = c(
    "with" = "#528B8B",
    "without" = "#CD950C"
  )) +
  
  # model predictions
  geom_point(
    position = pd,
    size = 4,
    stroke = 1.2
  ) +
  
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.2,
    position = pd
  ) +
  
  facet_grid(region ~ stage) +
  
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
    y = "Onset (GDD)",
    title = "Effect of transplantation and warming on onset",
    shape = "Treatment site × warming",
    color = "Biotic interactions",
    fill = "Biotic interactions"
  ) +
  guides(shape = "none")
b_f_fr_gdd

# ggsave(filename = "Output/Onset/GDD_Transplantation_onset_bud_flower_fruit_seeds_predictions_violin.png", 
#        plot = b_f_fr_gdd, width = 18, height = 10, units = "in")





b_f_fr_gdd2 <- ggplot(plot_df_all_gdd, aes(
  x = treat_competition,
  y = fit,
  color = treat_competition,
  shape = treatment_site_temp
)) +
  
  # raw:
  geom_violin(
    data = plot_df_raw_all_gdd,
    aes(
      x = treat_competition,
      y = onset,
      fill = treat_competition
      #group = interaction(treatment_site_temp, treat_competition)
    ),
    position = pd,
    alpha = 0.25,
    color = NA,
    trim = FALSE,
    adjust = 1.5 # makes it smoother so it looks less bubbly
  ) +
  
  scale_fill_manual(values = c(
    "with" = "#528B8B",
    "without" = "#CD950C"
  )) +
  
  # model predictions
  geom_point(
    position = pd,
    size = 4,
    stroke = 1.2
  ) +
  
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.2,
    position = pd
  ) +
  
  facet_grid(region ~ stage) +
  
  scale_color_manual(values = c(
    "with" = "#528B8B",
    "without" = "#CD950C"
  )) +
  
  scale_shape_manual(values = c(
    "lo_ambi" = 16,
    "hi_ambi" = 17
  )) +
  
  labs(
    x = "Site temperature treatment",
    y = "Onset (GDD)",
    title = "Effect of transplantation and warming on onset",
    shape = "Treatment site × warming",
    color = "Biotic interactions",
    fill = "Biotic interactions"
  ) 
b_f_fr_gdd2










