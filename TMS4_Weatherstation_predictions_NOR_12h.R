

# 12 h - only daytime -----------------------------------------------------

# Prepare TMS4 data for GDD - NOR -------------------------------------------------

# predict tomst data with climate station data ----------------------------

# Calculate GDD per treatment ---------------------------------------------


## Data used: Data/Data_tomst_loggers/tomst_2023/,
##            tomst_plot_codes_2023.csv,
##            RangeX_clean_MetadataPlot_NOR.csv
## Date:      19.06.2026
## Author:    Nadine Arzt
## Purpose:   Predict tms data for the beginning of the year with weather station data
##            Make one dataset with used daily mean per treatment
##            Calculate GDD, tbase = 2



# load library ------------------------------------------------------------
library(conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
library(tidyverse)
library(performance)
library(lubridate)


# source tms and climate station scripts ----------------------------------
source("RangeX_data_paper_cleaning_tomst_2023.R")

source("Data_preparation_climate_station_NOR.R")
# use climate

theme_set(theme_bw())


# climate station ---------------------------------------------------------

# filter year 23 -------------------------------------------------------
climate_23 <- climate |> 
  filter(year == 2023)

climate_23 <- climate_23 |>
  mutate(
    date_time = lubridate::parse_date_time(
      date_time,
      orders = c("ymd HMS", "ymd")))

climate_23_daily <- climate_23 |> 
  mutate(date = as_date(date_time),
         hour = hour(date_time)) |> 
  group_by(site, date) |> 
  summarise(
    
    Tmean_24h = mean(AirTemp_Avg, na.rm = TRUE),
    Tmean_12h = mean(AirTemp_Avg[hour >=8 & hour < 20], 
                     na.rm = TRUE),
    
    Tmax = max(AirTemp_Avg, na.rm = TRUE),
    Tmin = min(AirTemp_Avg, na.rm = TRUE),
    Humidity = mean(Humidity_Avg, na.rm = TRUE),
    WindSpd = mean(WindSpd_Avg, na.rm = TRUE),
    Radiation = mean(Radiation_Avg, na.rm = TRUE),
    Rainfall = sum(Rainfall, na.rm = TRUE),
    
    .groups = "drop")


# filter only year 2023 -----------------------------------------------------
# has end of 2022 in 
climate_23_daily <- climate_23_daily |> 
  filter(date >= as.Date("2023-01-01"))



# tomst -------------------------------------------------------------------

# calculate Temp average per day -----------------
tomst_nor_daily <- rx_tomst_23_clean |>
  filter(!is.na(date_time)) |>
  mutate(date = as.Date(date_time),
         hour = hour(date_time)) |>
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
    T3_mean_24h = mean(TMS_T3, na.rm = TRUE),
    T3_mean_12h = mean(TMS_T3[hour >= 8 & hour < 20], 
                       na.rm = TRUE),
    .groups = "drop"
  )




# combine climate station and real tomst data -----------------------------
# with treatments
compare_temp <- tomst_nor_daily |>
  inner_join(climate_23_daily, by = c("date", "site"))



ggplot(compare_temp, aes(x = Tmean_12h, y = Tmean_12h)) +
  geom_point()+
  facet_wrap(~ site)





# one model per logger ----------------------------------------------------
compare_temp |>
  distinct(
    unique_plot_ID,
    site,
    treat_warming,
    treat_competition
  ) |>
  arrange(unique_plot_ID)

# make one model per logger = unique_plot_ID
# with daily mean temp from the climate station and radiation and wind speed
logger_models_12 <- compare_temp |>
  group_by(unique_plot_ID) |>
  group_split() |>
  set_names(
    compare_temp |>
      distinct(unique_plot_ID) |>
      pull(unique_plot_ID)
  ) |>
  map(~ lm(Tmean_12h ~ Tmean_12h + Radiation + WindSpd, data = .x))


names(logger_models_12)

# check some outputs
summary(logger_models_12[[1]])
summary(logger_models_12[[3]])
summary(logger_models_12[[9]])
summary(logger_models_12[[20]])

logger_data_12 <- compare_temp |>
  filter(unique_plot_ID == names(logger_models_12)[1])

head(logger_data_12)





# predict per logger and loop through all ---------------------------------

logger_ids <- unique(compare_temp$unique_plot_ID)

logger_ids
length(logger_ids)

# function to predict logger data
predict_logger_12 <- function(id){
  
  # data from one logger
  dat_logger <- compare_temp |>
    filter(unique_plot_ID == id)
  
  # model
  m_logger_12 <- lm(
    T3_mean_12h ~ Tmean_12h + Radiation + WindSpd,
    data = dat_logger
  )
  
  # metadata (one row per logger)
  meta <- dat_logger |>
    distinct(unique_plot_ID, site, treat_warming, treat_competition, added_focals)
  
  # climate data for that site
  newdata <- climate_23_daily |>
    filter(site == meta$site)
  
  # prediction
  pred <- predict(
    m_logger_12,
    newdata = newdata,
    se.fit = TRUE
  )
  
  # output
  tibble(
    unique_plot_ID = id,
    site = meta$site,
    treat_warming = meta$treat_warming,
    treat_competition = meta$treat_competition,
    added_focals = meta$added_focals,
    date = newdata$date,
    fit = pred$fit,
    se.fit = pred$se.fit
  )
}

# does the prediction
predict_logger(logger_ids[1])


# loops through all loggers
pred_logger_12 <- map_df(logger_ids, predict_logger_12)




# combine real tomst with predicted tomst -------------------------------------------------
# to check how well it is predicted in the part that is overlapping

# add column with type (pred or real) -------------------------------------
# predicted tms data
pred_logger2_12 <- pred_logger_12 |>
  mutate(logger_type = "tms_predicted")|> 
  rename(temp_mean = fit)

tms_pred_clean_12 <- pred_logger2_12 |>
  select(date, site, treat_warming, treat_competition, added_focals, temp_mean, unique_plot_ID, logger_type)


# real tms data
tms_real_12 <- tomst_nor_daily |> 
  mutate(logger_type = "tms_real")|>
  rename(temp_mean = T3_mean_12h) |> 
  select(date, site, treat_warming, treat_competition, added_focals, temp_mean, unique_plot_ID, logger_type)


# tomst_range_by_treat <- tomst_nor_daily |>
#   group_by(site, treat_warming, treat_competition) |>
#   summarise(
#     start_date = min(date, na.rm = TRUE),
#     end_date   = max(date, na.rm = TRUE),
#     n_days      = n(),
#     .groups = "drop"
#   )
# 
# tomst_range_by_treat
# 
# 
# logger_ranges <- tms_real |>
#   group_by(
#     unique_plot_ID,
#     site,
#     treat_warming,
#     treat_competition
#   ) |>
#   summarise(
#     first_date = min(date),
#     last_date = max(date),
#     n_days = n_distinct(date),
#     .groups = "drop"
#   ) |>
#   arrange(
#     site,
#     treat_warming,
#     treat_competition,
#     first_date
#   )
# logger_ranges
# 
# logger_ranges |>
#   group_by(
#     site,
#     treat_warming,
#     treat_competition
#   ) |>
#   summarise(
#     n_loggers = n(),
#     earliest_start = min(first_date),
#     latest_start = max(first_date),
#     .groups = "drop"
#   )


# join tms predicted and real --------------------------------------------- 
tms_pred_real_12 <- bind_rows(tms_pred_clean_12, tms_real_12)

# plot predicted and real -------------------------------------------------
ggplot(tms_pred_real_12,
       aes(x = date, y = temp_mean, color = logger_type, linetype = site)) +
  geom_line() +
  facet_grid(treat_warming  + added_focals ~ treat_competition)


ggplot(tms_pred_real_12,
       aes(x = date, y = temp_mean, color = logger_type, linetype = site)) +
  geom_line() +
  facet_grid(treat_warming + site ~ treat_competition)




# make one dataset with predicted first and then real tms --------------------

# combine pred_logger with tomst_nor_daily

tms_pred_12 <- pred_logger_12 |>
  mutate(logger_type = "tms_predicted")|> 
  rename(temp_mean_pred = fit) |> 
  select(date, site, treat_warming, treat_competition, added_focals, unique_plot_ID, temp_mean_pred, logger_type)


tms_real_12 <- tomst_nor_daily |> 
  mutate(logger_type = "tms_real")|>
  rename(temp_mean_real = T3_mean_12h)|> 
  select(date, site, treat_warming, treat_competition, added_focals, unique_plot_ID, temp_mean_real, logger_type)


tms_joined_12 <- tms_pred_12 |>
  left_join(
    tms_real_12,
    by = c(
      "date",
      "site",
      "treat_warming",
      "treat_competition",
      "added_focals",
      "unique_plot_ID"
    )
  )

# make new column chosing pred or real and correct logger type column
tms_joined_12 <- tms_joined_12 |>
  mutate(temp_mean_used = coalesce(temp_mean_real, temp_mean_pred),
         
         logger_type = if_else(
           !is.na(temp_mean_real),
           "tms_real",
           "tms_predicted"))


# select only necessary columns
tms_joined2_12 <- tms_joined_12 |>
  select(
    date,
    site,
    treat_warming,
    treat_competition,
    added_focals,
    unique_plot_ID,
    temp_mean_pred,
    temp_mean_real,
    temp_mean_used,
    logger_type)



ggplot(tms_joined2_12,
       aes(x = date, y = temp_mean_used, color = logger_type)) +
  geom_line() +
  facet_grid(treat_warming + site ~ treat_competition)




tms_joined2_12 |>
  count(
    unique_plot_ID,
    date
  ) |>
  filter(n > 1)

tms_joined2_12 |>
  count(
    date,
    site,
    treat_warming,
    treat_competition,
    added_focals
  ) |>
  arrange(desc(n))

# this plot shows per logger 
# that is why it has some overlap where it is predicted and real at the same time
# whithin one treatment, loggers were deployed at several dates
ggplot(tms_joined2_12,
       aes(date,
           temp_mean_used,
           colour = logger_type,
           group = unique_plot_ID)) +
  geom_line()+
  facet_grid(treat_warming + site ~ treat_competition)



# calculate mean per treatment --------------------------------------------
# average of several loggers
tms_final_12 <- tms_joined2_12 |>
  group_by(
    date,
    site,
    treat_warming,
    treat_competition,
    added_focals
  ) |>
  summarise(
    temp_mean = mean(temp_mean_used, na.rm = TRUE),
    .groups = "drop"
  )



# Calculate daily mean per treatment including percentage of real data --------
# get number of loggers that are using real and predicted tms data
# and calculate percentage for both
tms_final_12 <- tms_joined2_12 |>
  group_by(
    date,
    site,
    treat_warming,
    treat_competition,
    added_focals
  ) |>
  summarise(
    temp_mean = mean(temp_mean_used, na.rm = TRUE),
    
    n_loggers = n(),
    
    n_real = sum(logger_type == "tms_real"),
    n_pred = sum(logger_type == "tms_predicted"),
    
    perc_real = 100 * n_real / n_loggers,
    perc_pred = 100 * n_pred / n_loggers,
    
    .groups = "drop"
  )




# plot used temperature mean per treatment --------------------------------
ggplot(tms_final_12,
       aes(date,
           temp_mean)) +
  geom_line()+
  facet_grid(treat_warming + site + added_focals ~ treat_competition)


ggplot(tms_final_12,
       aes(date,
           temp_mean, colour = treat_competition, linetype = site)) +
  geom_line()+
  facet_grid(treat_warming + added_focals ~ treat_competition)



# column with site_warming treatment --------------------------------------
tms_final_12$treatment_site_temp <- paste(tms_final_12$site, tms_final_12$treat_warming, sep = "_")


tms_final_12 <- tms_final_12|>
  mutate(treatment_site_temp= factor(treatment_site_temp,
                                     levels = c("lo_ambi",
                                                "hi_ambi",
                                                "hi_warm")))


# control plot
ggplot(tms_final_12,
       aes(date,
           temp_mean, colour = treat_competition, linetype = treatment_site_temp)) +
  geom_line()+
  facet_grid(added_focals ~treatment_site_temp)



# delete control plots ----------------------------------------------------
# we only need the plots with focals for now
tms_final_nor_12 <- tms_final_12 |> 
  filter(added_focals == "wf")

tms_final_nor_12$treatment_site_temp_comp <- paste(tms_final_nor_12$site, 
                                                   tms_final_nor_12$treat_warming,
                                                   tms_final_nor_12$treat_competition, sep = "_")



ggplot(tms_final_nor_12,
       aes(date,
           temp_mean, colour = treat_competition)) +
  geom_line()+
  facet_grid( ~treatment_site_temp_comp)



# Calculate GDD -----------------------------------------------------------

# define Tbase and growing season start  -------------------------
Tbase <- 2 # base temperature for plants to grow
Nconsec <- 5 # number of consecutive days above Tbase to define growing season start
# 5 days above 5 degrees



# 1) indicator: Tavg_tms > Tbase
clim_flag_nor_tms_12 <- tms_final_nor_12 |> 
  arrange(site, treat_warming, treat_competition, date) |> 
  mutate(is_warm = (temp_mean > Tbase))
clim_flag_nor_tms_12

# find out start of growing season -------------------------------------
# 2) compute run of Nconsec TRUE per site and find first date
season_start_nor_tms_12 <- clim_flag_nor_tms_12 |> 
  filter(date >= as.Date("2023-01-01")) |> 
  group_by(site, treat_warming, treat_competition) |> 
  mutate(run_N = slider::slide_dbl(is_warm, ~ ifelse(sum(.x) == length(.x), 1, 0),
                                   .before = Nconsec - 1, .complete = TRUE)) |> 
  # run_N == 1 on the last day of a full Nconsec-run of TRUEs
  filter(run_N == 1) |> 
  summarise(season_start = min(date), .groups = "drop")
season_start_nor_tms_12


# 1 hi    ambi          bare              2023-04-12  
# 2 hi    ambi          vege              2023-04-12  
# 3 hi    warm          bare              2023-04-11  
# 4 hi    warm          vege              2023-04-11  
# 5 lo    ambi          bare              2023-04-11  
# 6 lo    ambi          vege              2023-04-11  


# calculate gdd per site ----------------------------------------------
# 3) join start dates back and calculate GDD from that start in april
climate_gdd_nor_tms_12 <- tms_final_nor_12 |> 
  left_join(season_start_nor_tms_12, by = c("site",
                                         "treat_warming",
                                         "treat_competition")) |> 
  filter(!is.na(season_start)) |> 
  group_by(site,  treat_warming, treat_competition) |> 
  arrange(date) |> 
  mutate(
    # only accumulate on/after season_start
    GDD_day = if_else(date >= season_start, pmax(0, temp_mean - Tbase), 0),
    GDD_cum = cumsum(GDD_day)
  ) |> 
  ungroup()
climate_gdd_nor_tms_12

# rename date to date_measurement --------------------------------------
# to match phenology
climate_gdd_nor_tms_12 <- climate_gdd_nor_tms_12 |> 
  rename("date_measurement" = "date")


# get julian day ----------------------------------------------------------
climate_gdd_nor_tms_12 <- climate_gdd_nor_tms_12 |>
  mutate(jday = lubridate::yday(date_measurement))


# control plot ------------------------------------------------------------
ggplot(climate_gdd_nor_tms_12,
       aes(x = date_measurement,
           y = GDD_cum,
           color = treatment_site_temp)) +
  geom_line(size = 1.2) +
  labs(
    x = "Date",
    y = "Cumulative GDD",
    title = "Norway cumulative GDD development")+
  facet_grid(~ treat_competition)



ggplot(climate_gdd_nor_tms_12,
       aes(x = jday,
           y = GDD_cum,
           color = treatment_site_temp_comp,
           alpha = perc_real)) + # includes percentage of real data
  geom_line(size = 1.2) +
  labs(
    x = "Day of year (DOY)",
    y = "Cumulative temperature (GDD2)",
    title = "Norway cumulative GDD")









