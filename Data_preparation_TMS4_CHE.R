

# Data preparation tomst CHE ------------------------------------

# Calculate GDD per treatment ---------------------------------------------


## Data used: RangeX_clean_EnvTMS4_2021_2023_CHE.csv,
##            RangeX_clean_MetadataPlot_CHE.csv
## Date:      19.06.2026
## Author:    Nadine Arzt
## Purpose:   Make one dataset with used daily mean per treatment
##            Calculate GDD, tbase = 2


# load library ------------------------------------------------------------
library(conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
library(tidyverse)
library(lubridate)
library(slider) 



# import tomst data -----------------------------------------------------
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
  group_by(site, unique_plot_ID, treat_warming, treat_competition, added_focals) |>
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



# calculate daily mean per treatment --------------------------------------
tms_che_treat <- env_che_22 |>
  group_by(date, site, treat_warming, treat_competition, added_focals) |>
  summarise(
    temp_mean = mean(TMS_T3, na.rm = TRUE),
    .groups = "drop")
tms_che_treat

# cut away the first days of the year to not have GDD start in this warm period before being cold again
tms_che_treat <- tms_che_treat |>
  filter(date >= as.Date("2022-01-07"),
         date<= as.Date("2022-10-23")) # same day as in Norway


# delete control plots ----------------------------------------------------
# we only need the plots with focals for now
tms_final_che <- tms_che_treat |> 
  filter(added_focals == "wf")

tms_final_che$treatment_site_temp_comp <- paste(tms_final_che$site, tms_final_che$treat_warming,tms_final_che$treat_competition, sep = "_")

# column with site_warming treatment --------------------------------------
tms_final_che$treatment_site_temp <- paste(tms_final_che$site, tms_final_che$treat_warming, sep = "_")


tms_final_che <- tms_final_che |>
  mutate(treatment_site_temp= factor(treatment_site_temp,
                                     levels = c("lo_ambi",
                                                "hi_ambi",
                                                "hi_warm")))

ggplot(tms_final_che,
       aes(date,
           temp_mean, colour = treat_competition)) +
  geom_line()+
  facet_grid( ~treatment_site_temp_comp)




# Calculate GDD per treatment ---------------------------------------------

# define Tbase and growing season start  -------------------------
Tbase <- 2 # base temperature for plants to grow
Nconsec <- 5 # number of consecutive days above Tbase to define growing season start
# 5 days above 5 degrees



# 1) indicator: Tavg_tms > Tbase
clim_flag_che_tms <- tms_che_treat |> 
  arrange(site, treat_warming, treat_competition, date) |> 
  mutate(is_warm = (temp_mean > Tbase))


# find out start of growing season -------------------------------------
# 2) compute run of Nconsec TRUE per site and find first date
season_start_che_tms <- clim_flag_che_tms |> 
  filter(date >= as.Date("2022-01-07")) |> 
  group_by(site, treat_warming, treat_competition) |> 
  mutate(run_N = slider::slide_dbl(is_warm, ~ ifelse(sum(.x) == length(.x), 1, 0),
                                   .before = Nconsec - 1, .complete = TRUE)) |> 
  # run_N == 1 on the last day of a full Nconsec-run of TRUEs
  filter(run_N == 1) |> 
  summarise(season_start = min(date), .groups = "drop")
season_start_che_tms


# 1 hi    ambi          bare              2022-04-16  
# 2 hi    ambi          vege              2022-04-14  
# 3 hi    warm          bare              2022-04-16  
# 4 hi    warm          vege              2022-04-13  
# 5 lo    ambi          bare              2022-03-27  
# 6 lo    ambi          vege              2022-03-27



# calculate gdd per treat ----------------------------------------------
# 3) join start dates back and calculate GDD from that start in april
climate_gdd_che_tms <- tms_final_che |> 
  left_join(season_start_che_tms, by = c("site", "treat_warming", "treat_competition")) |> 
  filter(!is.na(season_start)) |> 
  group_by(site, treat_warming, treat_competition) |> 
  arrange(date) |> 
  mutate(
    # only accumulate on/after season_start
    GDD_day = if_else(date >= season_start, pmax(0, temp_mean - Tbase), 0),
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

ggplot(climate_gdd_che_tms,
       aes(x = jday,
           y = GDD_cum,
           color = treatment_site_temp_comp)) +
  geom_line(size = 1.2) +
  labs(
    x = "Day of year (DOY)",
    y = "Cumulative temperature (GDD2)",
    title = "Switzerland cumulative GDD")



# daily means -------------------------------------------------------------

# calculate Temp average per day -----------------
tomst_che_daily <- env_che_22 |>
  filter(!is.na(date_time)) |>
  mutate(date = as.Date(date_time)) |>
  group_by(
    date,
    site,
    treat_warming,
    treat_competition,
    added_focals,
    unique_plot_ID
  ) |>
  summarise(
    T1_mean = mean(TMS_T1, na.rm = TRUE),
    T2_mean = mean(TMS_T2, na.rm = TRUE),
    T3_mean = mean(TMS_T3, na.rm = TRUE),
    .groups = "drop"
  )


