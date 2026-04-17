


# Temperature sensitivity -------------------------------------------------

# load library ---------------------------------------------------------
library(lme4)
library(ggeffects)
library(broom.mixed)
library(emmeans)
library(lubridate)
library(performance)
library(see)


# calcualte delta T -------------------------------------------------------

# mean temperature per site
# timeframe before the first phenological event

# can only do from mid april because che starts from 13.04.22

# when does onset start at each site
# first_onset_bud
# # A tibble: 4 × 4
# region      site   year first_onset
# <chr>       <chr> <dbl>       <dbl>
# 1 Norway      hi     2023         158 = 07.06
# 2 Norway      lo     2023         132 = 12.05
# 3 Switzerland hi     2022         134 = 14.05
# 4 Switzerland lo     2022         124 = 04.05

# first_onset_flower
# # A tibble: 4 × 4
# region      site   year first_onset
# <chr>       <chr> <dbl>       <dbl>
# 1 Norway      hi     2023         158 = 07.06
# 2 Norway      lo     2023         144 = 24.05
# 3 Switzerland hi     2022         134 = 14.05
# 4 Switzerland lo     2022         124 = 04.05

# first_onset_fruit
# # A tibble: 4 × 4
# region      site   year first_onset
# <chr>       <chr> <dbl>       <dbl>
# 1 Norway      hi     2023         171
# 2 Norway      lo     2023         170
# 3 Switzerland hi     2022         150
# 4 Switzerland lo     2022         135 = 15.05


# source climate scripts --------------------------------------------------
source("Data_preparation_climate_station_CHE.R")

source("Data_preparation_climate_station_NOR.R")

# combine nor and che -----------------------------------------------------
climate_23
climate_che_22

# filter time period
climate_nor_23_pre <- climate_23 |>
  filter(
    date_time >= as.Date("2023-05-01"),
    date_time <= as.Date("2023-05-31")
  )

climate_che_22_pre <- climate_che_22 |>
  filter(
    date_time >= as.Date("2022-05-14"),
    date_time <= as.Date("2022-06-13")
  )


# combine
climate_all <- bind_rows(climate_che_22_pre, climate_nor_23_pre)

climate_all <- climate_all |> 
  select(region, site, year, date_time, AirTemp_Avg, Humidity_Avg, WindDir_Avg, WindSpd_Avg, Radiation_Avg, Rainfall)


# calculate mean per site in this time period
pre_climate <-  climate_all |>
  mutate(date = as.Date(date_time)) |>
  group_by(region, site, year, date) |>
  summarise(T_day = mean(AirTemp_Avg, na.rm = TRUE)) |>
  group_by(region, site, year) |>
  summarise(Tmean = mean(T_day))
pre_climate

pre_climate <- pre_climate |>
  mutate(region = case_when(
    region == "NOR" ~ "Norway",
    region == "CHE" ~ "Switzerland",
    TRUE ~ region
  ))
pre_climate



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
    jday_scaled = scale(jday))        # optional scaling if needed



# filter only ambi both sites  -------------------------------------
# to compare low ambi with hi ambi = cooling effect
phenology3 <- phenology3 |> 
  filter(treat_warming == "ambi")



phenology3 <- phenology3 |> 
  mutate(year = if_else(region == "Switzerland", 2022, 2023))


# calculate onsets ---------------------------------------------------------


# bud ---------------------------------------------------------------------
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


# flower ---------------------------------------------------------------------
onset_flower <- phenology3 |>
  filter(phenology_stage == "No_FloOpen", value > 0) |>
  group_by(region, site, year, treat_competition, species, block_ID, unique_plot_ID, unique_plant_ID, phenology_stage) |>
  summarise(onset = min(jday), .groups = "drop") |>
  # remove groups where flowering never occurred
  filter(is.finite(onset))

first_onset_flower <- onset_flower |>
  group_by(region, site, year) |>
  summarise(first_onset = min(onset, na.rm = TRUE),
            .groups = "drop")
first_onset_flower

# fruit ---------------------------------------------------------------------
# combine the fruiting stages nor and che to comapre the onset ------------
# this is to be taken with caution because the stages are not the same
# but for the onset it could be comparable
phenology3 <- phenology3 |>
  mutate(phenology_stage = recode(
    phenology_stage,
    "No_Infructescences" = "No_FloWithrd"
  ))

onset_fruit <- phenology3 |>
  filter(phenology_stage == "No_FloWithrd", value > 0) |>
  group_by(region, site, year, treat_competition, species, block_ID, unique_plot_ID, unique_plant_ID, phenology_stage) |>
  summarise(onset = min(jday), .groups = "drop") |>
  # remove groups where fruiting never occurred
  filter(is.finite(onset))

first_onset_fruit <- onset_fruit |>
  group_by(region, site, year) |>
  summarise(first_onset = min(onset, na.rm = TRUE),
            .groups = "drop")
first_onset_fruit

# seed ---------------------------------------------------------------------
onset_seed <- phenology3 |>
  filter(phenology_stage == "No_Seeds", value > 0) |>
  group_by(region, site, year, treat_competition, species, block_ID, unique_plot_ID, unique_plant_ID, phenology_stage) |>
  summarise(onset = min(jday), .groups = "drop") |>
  # remove groups where seeds never occurred
  filter(is.finite(onset))




# Average onset per site --------------------------------------------------
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



# Calculate temperature sensitivity ---------------------------------------


# bud ---------------------------------------------------------------------
sens_bud <- onset_bud_mean |>
  left_join(pre_climate, by = c("region","site", "year")) |>
  pivot_wider(names_from = site,
              values_from = c(onset, Tmean)) |>
  mutate(
    temp_sens = (onset_lo - onset_hi) / (Tmean_lo - Tmean_hi)
  )
sens_bud


sens_flower <- onset_flower_mean |>
  left_join(pre_climate, by = c("region","site", "year")) |>
  pivot_wider(names_from = site,
              values_from = c(onset, Tmean)) |>
  mutate(
    temp_sens = (onset_lo - onset_hi) / (Tmean_lo - Tmean_hi)
  )
sens_flower

sens_fruit <- onset_fruit_mean |>
  left_join(pre_climate, by = c("region","site", "year")) |>
  pivot_wider(names_from = site,
              values_from = c(onset, Tmean)) |>
  mutate(
    temp_sens = (onset_lo - onset_hi) / (Tmean_lo - Tmean_hi)
  )
sens_fruit

sens_seed <- onset_seed_mean |>
  left_join(pre_climate, by = c("region","site", "year")) |>
  pivot_wider(names_from = site,
              values_from = c(onset, Tmean)) |>
  mutate(
    temp_sens = (onset_lo - onset_hi) / (Tmean_lo - Tmean_hi)
  )
sens_seed



sens_all <- bind_rows(
  sens_bud   |> mutate(stage = "bud"),
  sens_flower|> mutate(stage = "flower"),
  sens_fruit |> mutate(stage = "fruit"),
  sens_seed  |> mutate(stage = "seed")
)
sens_all





# plot  -------------------------------------------------------------------

ggplot(sens_all,
       aes(x = stage,
           y = temp_sens,
           color = treat_competition)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.4),
             alpha = 0.7) 

sens_summary <- sens_all |>
  group_by(region, stage, treat_competition) |>
  summarise(
    mean = mean(temp_sens, na.rm = TRUE),
    se   = sd(temp_sens, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )
sens_summary


region_colors <- c(
  "Norway" = "turquoise4",
  "Switzerland" = "pink4"
)

ggplot(sens_summary,
       aes(x = treat_competition,
           y = mean,
           color = region)) +
  geom_point(position = position_dodge(0.4), size = 3) +
  geom_errorbar(
    aes(ymin = mean - 1.96 * se,
        ymax = mean + 1.96 * se),
    width = 0.2,
    position = position_dodge(0.4)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed")+
  facet_grid( ~ stage)+
  labs(y = "Temperature sensitivity")+
  scale_color_manual(values = region_colors)











