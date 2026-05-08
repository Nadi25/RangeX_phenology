


# Data preparation climate station CHE ------------------------------------


# load library ------------------------------------------------------------
library(conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
library(tidyverse)
library(lubridate)
library(slider) 



# import climate data -----------------------------------------------------
climate_che <- read.csv("Data/Clean/RangeX_WeatherStation_CleanData_20251227_CHE.csv")

climate_che <- climate_che |> 
  select(-X)


# filter 2022 -------------------------------------------------------------
climate_che_22 <- climate_che |> 
  filter(year == 2022)

# calculate Tmax, Tmin and Temp average per day -----------------
climate_22_daily <- climate_che_22 |> 
  mutate(date = as_date(date_time)) |> 
  group_by(site, date) |> 
  summarise(Tmax = max(AirTemp_Avg, na.rm = TRUE),
            Tmin = min(AirTemp_Avg, na.rm = TRUE),
            .groups = "drop") |> 
  mutate(Tavg = (Tmax + Tmin) / 2)



# define Tbase and growing season start  -------------------------
Tbase <- 2 # base temperature for plants to grow
Nconsec <- 5 # number of consecutive days above Tbase to define growing season start
# 5 days above 5 degrees


# 1) indicator: Tavg > Tbase
clim_flag_che <- climate_22_daily |> 
  arrange(site, date) |> 
  mutate(is_warm = (Tavg > Tbase))


# find out start of growing season -------------------------------------
# 2) compute run of Nconsec TRUE per site and find first date
season_start_che <- clim_flag_che |> 
  group_by(site) |> 
  mutate(run_N = slider::slide_dbl(is_warm, ~ ifelse(sum(.x) == length(.x), 1, 0),
                                   .before = Nconsec - 1, .complete = TRUE)) |> 
  # run_N == 1 on the last day of a full Nconsec-run of TRUEs
  filter(run_N == 1) |> 
  summarise(season_start = min(date), .groups = "drop")
season_start_che

# 1 hi    2022-05-18  
# 2 lo    2022-04-22  

# calculate gdd per site ----------------------------------------------
# 3) join start dates back and calculate GDD from that start in april
climate_gdd_che <- climate_22_daily|> 
  left_join(season_start_che, by = "site") |> 
  filter(!is.na(season_start)) |> 
  group_by(site) |> 
  arrange(date) |> 
  mutate(
    # only accumulate on/after season_start
    GDD_day = if_else(date >= season_start, pmax(0, Tavg - Tbase), 0),
    GDD_cum = cumsum(GDD_day)
  ) |> 
  ungroup()
climate_gdd_che

# rename date to date_measurement --------------------------------------
# to match phenology
climate_gdd_che <- climate_gdd_che |> 
  rename("date_measurement" = "date")




# control plot ------------------------------------------------------------
ggplot(climate_gdd_che,
       aes(x = date_measurement,
           y = GDD_cum,
           color = site)) +
  geom_line(size = 1.2) +
  
  geom_vline(
    xintercept = as.Date(c("2022-05-04", "2022-05-14")),
    linetype = "dashed"
  ) +
  
  labs(
    x = "Date",
    y = "Cumulative GDD",
    title = "Swiss cumulative GDD development"
  ) 







