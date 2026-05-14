



# Data preparation tomst CHE ------------------------------------


# load library ------------------------------------------------------------
library(conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
library(tidyverse)
library(lubridate)
library(slider) 



# import climate data -----------------------------------------------------
tms_che <- read.csv("Data/RangeX_clean_EnvTMS4_2021_2023_CHE.csv")



# import metadata ---------------------------------------------------------
meta_p_che <- read.csv("Data/RangeX_clean_MetadataPlot_CHE.csv")


# combine metadata with tms data ------------------------------------------
env_che <- tms_che |> 
  left_join(meta_p_che, by = "unique_plot_ID")


# Date as date ------------------------------------------------------------
env_che <- env_che |> 
  mutate(date_time = parse_date_time(date_time, orders = c("Y-m-d H:M:S", "Y-m-d H:M", "Y-m-d")))

# check if CHE has NAs in date_time ------------------------------------------------
date_na <- env_che |> 
  filter(is.na(date_time))
date_na

# plot_ID_original as chr -------------------------------------------------
env_che <- env_che |> mutate(plot_ID_original = as.character(plot_ID_original))

# block_ID_original as chr -------------------------------------------------
env_che <- env_che |> mutate(block_ID_original  = as.character(block_ID_original ))


# block_ID as chr -------------------------------------------------
env_che <- env_che |> mutate(block_ID  = as.character(block_ID ))

env_che_22 <- env_che |>
  mutate(
    date = as.Date(date_time),
    year = lubridate::year(date_time)
  ) |>
  filter(year == 2022)


# data availability 2022 ---------------------------------------------------
env_che_22_ <- env_che |>
  mutate(
    date = as.Date(date_time),
    year = lubridate::year(date_time)
  ) |>
  filter(year == 2022) |>
  group_by(site, unique_plot_ID, treat_warming, treat_competition) |>
  summarise(
    start_date = min(date, na.rm = TRUE),
    end_date   = max(date, na.rm = TRUE),
    n_records  = n(),
    .groups = "drop"
  ) |>
  arrange(site, start_date)
env_che_22_

# 2022-01-01 2022-12-31
# for 19 plots

env_che_22_ |>
  count(site, treat_warming, treat_competition)

# 1 hi    ambi          bare                  1
# 2 hi    ambi          vege                  4
# 3 hi    warm          bare                  3
# 4 hi    warm          vege                  4
# 5 lo    ambi          bare                  4
# 6 lo    ambi          vege                  3




# Calculate GDD per site --------------------------------------------------
# daily mean temp

# use T3 which is at 15 cm

# climate_tomst_22_daily <- env_che_22 |>
#   filter(
#     treat_warming == "ambi",
#     treat_competition == "vege"
#   ) |>
#   group_by(region, site, treat_warming, treat_competition, added_focals, block_ID, unique_plot_ID, date) |>
#   summarise(
#     Tmax = max(TMS_T3, na.rm = TRUE),
#     Tmin = min(TMS_T3, na.rm = TRUE),
#     .groups = "drop") |>
#   mutate(Tavg_tms = (Tmax + Tmin)/2)
# climate_tomst_22_daily

# calculate mean per site per day
# for only plots with focals to make it comaprable hi and lo
climate_tomst_22_site <- env_che_22 |>
  filter(
    treat_warming == "ambi",
    treat_competition == "vege",
    added_focals == "wf"
  ) |>
  group_by(site, date) |>
  summarise(
    Tmax = max(TMS_T3, na.rm = TRUE),
    Tmin = min(TMS_T3, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    Tavg_tms = (Tmax + Tmin)/2
  )

climate_tomst_22_site <- climate_tomst_22_site |>
  filter(date >= as.Date("2022-01-07"))



# define Tbase and growing season start  -------------------------
Tbase <- 2 # base temperature for plants to grow
Nconsec <- 5 # number of consecutive days above Tbase to define growing season start
# 5 days above 5 degrees



# 1) indicator: Tavg_tms > Tbase
clim_flag_che_tms <- climate_tomst_22_site |> 
  arrange(site, date) |> 
  mutate(is_warm = (Tavg_tms > Tbase))


# find out start of growing season -------------------------------------
# 2) compute run of Nconsec TRUE per site and find first date
season_start_che_tms <- clim_flag_che_tms |> 
  filter(date >= as.Date("2022-01-07")) |> 
  group_by(site) |> 
  mutate(run_N = slider::slide_dbl(is_warm, ~ ifelse(sum(.x) == length(.x), 1, 0),
                                   .before = Nconsec - 1, .complete = TRUE)) |> 
  # run_N == 1 on the last day of a full Nconsec-run of TRUEs
  filter(run_N == 1) |> 
  summarise(season_start = min(date), .groups = "drop")
season_start_che_tms


# 1 hi    2022-03-28  
# 2 lo    2022-03-26 


# calculate gdd per site ----------------------------------------------
# 3) join start dates back and calculate GDD from that start in april
climate_gdd_che_tms <- climate_tomst_22_site |> 
  left_join(season_start_che_tms, by = "site") |> 
  filter(!is.na(season_start)) |> 
  group_by(site) |> 
  arrange(date) |> 
  mutate(
    # only accumulate on/after season_start
    GDD_day = if_else(date >= season_start, pmax(0, Tavg_tms - Tbase), 0),
    GDD_cum = cumsum(GDD_day)
  ) |> 
  ungroup()
climate_gdd_che_tms

# rename date to date_measurement --------------------------------------
# to match phenology
climate_gdd_che_tms <- climate_gdd_che_tms |> 
  rename("date_measurement" = "date")


# get julian day ----------------------------------------------------------
climate_gdd_che_tms <- climate_gdd_che_tms |>
  mutate(jday = lubridate::yday(date_measurement))


# control plot ------------------------------------------------------------
ggplot(climate_gdd_che_tms,
       aes(x = date_measurement,
           y = GDD_cum,
           color = site)) +
  geom_line(size = 1.2) +
  
  geom_vline(
    xintercept = as.Date(c("2022-03-25", "2022-03-28")),
    linetype = "dashed"
  ) +
  
  labs(
    x = "Date",
    y = "Cumulative GDD",
    title = "Swiss cumulative GDD development"
  ) 







