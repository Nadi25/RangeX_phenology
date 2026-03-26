

# Prepare climate station data for growing degreee days -------------------


# import cliamte station data NOR -----------------------------------------
climate <- read.csv("Data/Clean/RangeX_clean_EnvClimateStation_2021-2025_NOR.csv")



# filter year 23 -------------------------------------------------------
climate_23 <- climate |> 
  filter(year == 2023)


# calculate Tmax, Tmin and Temp average per day -----------------
climate_23_daily <- climate_23 |> 
  mutate(date = as_date(date_time)) |> 
  group_by(site, date) |> 
  summarise(Tmax = max(AirTemp_Avg, na.rm = TRUE),
            Tmin = min(AirTemp_Avg, na.rm = TRUE),
            .groups = "drop") |> 
  mutate(Tavg = (Tmax + Tmin) / 2)



# filter only year 2023 -----------------------------------------------------
# has end of 2022 in 
climate_23_daily <- climate_23_daily |> 
  filter(date >= as.Date("2023-01-01"))


# define Tbase and growing season start  -------------------------
Tbase <- 5 # base temperature for plants to grow
Nconsec <- 5 # number of consecutive days above Tbase to define growing season start
# 5 days above 5 degrees


# 1) indicator: Tavg > Tbase
clim_flag <- climate_23_daily |> 
  arrange(site, date) |> 
  mutate(is_warm = (Tavg > Tbase))


# find out start of growing season -------------------------------------
# 2) compute run of Nconsec TRUE per site and find first date
season_start <- clim_flag |> 
  group_by(site) |> 
  mutate(run_N = slider::slide_dbl(is_warm, ~ ifelse(sum(.x) == length(.x), 1, 0),
                                   .before = Nconsec - 1, .complete = TRUE)) |> 
  # run_N == 1 on the last day of a full Nconsec-run of TRUEs
  filter(run_N == 1) |> 
  summarise(season_start = min(date), .groups = "drop")
season_start

# hi    2023-04-22  so its only two days later
# lo    2023-04-20 

# calculate gdd per site ----------------------------------------------
# 3) join start dates back and calculate GDD from that start in april
climate_gdd <- climate_23_daily |> 
  left_join(season_start, by = "site") |> 
  filter(!is.na(season_start)) |> 
  group_by(site) |> 
  arrange(date) |> 
  mutate(
    # only accumulate on/after season_start
    GDD_day = if_else(date >= season_start, pmax(0, Tavg - Tbase), 0),
    GDD_cum = cumsum(GDD_day)
  ) |> 
  ungroup()
climate_gdd


# rename date to date_measurement --------------------------------------
# to match phenology
climate_gdd <- climate_gdd |> 
  rename("date_measurement" = "date")


# filter by timeframe from phenology ----------------------------------
climate_gdd_pt <- climate_gdd |> 
  filter(date_measurement >= as.Date("2023-05-12"),
         date_measurement <= as.Date("2023-10-23"))

























