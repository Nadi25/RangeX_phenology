


# Prepare TMS4 data - NOR -------------------------------------------------

# predict tomst data with climate station data


# load library
library(conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
library(tidyverse)
library(performance)

# source tms and climate station scripts
source("RangeX_data_paper_cleaning_tomst_2023.R")

source("Data_preparation_climate_station_NOR.R")
# use climate

theme_set(theme_bw())


# filter year 23 -------------------------------------------------------
climate_23 <- climate |> 
  filter(year == 2023)


climate_23_daily <- climate_23 |> 
  mutate(date = as_date(date_time)) |> 
  group_by(site, date) |> 
  summarise(
    Tmax = max(AirTemp_Avg, na.rm = TRUE),
    Tmin = min(AirTemp_Avg, na.rm = TRUE),
    Tmean = mean(AirTemp_Avg, na.rm = TRUE),
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
  mutate(date = as.Date(date_time)) |>
  group_by(
    date,
    site,
    treat_warming,
    treat_competition,
    unique_plot_ID
  ) |>
  summarise(
    T1_mean = mean(TMS_T1, na.rm = TRUE),
    T2_mean = mean(TMS_T2, na.rm = TRUE),
    T3_mean = mean(TMS_T3, na.rm = TRUE),
    .groups = "drop"
  )




# combine climate station and real tomst data -----------------------------
# with treatments
compare_temp <- tomst_nor_daily |>
  inner_join(climate_23_daily, by = c("date", "site"))



ggplot(compare_temp, aes(x = Tmean, y = T3_mean)) +
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
logger_models <- compare_temp |>
  group_by(unique_plot_ID) |>
  group_split() |>
  set_names(
    compare_temp |>
      distinct(unique_plot_ID) |>
      pull(unique_plot_ID)
  ) |>
  map(~ lm(T3_mean ~ Tmean + Radiation + WindSpd, data = .x))


names(logger_models)

# check some outputs
summary(logger_models[[1]])
summary(logger_models[[3]])
summary(logger_models[[9]])
summary(logger_models[[20]])

logger_data <- compare_temp |>
  filter(unique_plot_ID == names(logger_models)[1])

head(logger_data)





# predict per logger and loop through all ---------------------------------

logger_ids <- unique(compare_temp2$unique_plot_ID)

logger_ids
length(logger_ids)

# function to predict logger data
predict_logger <- function(id){
  
  # data from one logger
  dat_logger <- compare_temp |>
    filter(unique_plot_ID == id)
  
  # model
  m_logger <- lm(
    T3_mean ~ Tmean + Radiation + WindSpd,
    data = dat_logger
  )
  
  # metadata (ONE row per logger)
  meta <- dat_logger |>
    distinct(unique_plot_ID, site, treat_warming, treat_competition)
  
  # climate data for that site
  newdata <- climate_23_daily |>
    filter(site == meta$site)
  
  # prediction
  pred <- predict(
    m_logger,
    newdata = newdata,
    se.fit = TRUE
  )
  
  # output
  tibble(
    unique_plot_ID = id,
    site = meta$site,
    treat_warming = meta$treat_warming,
    treat_competition = meta$treat_competition,
    date = newdata$date,
    fit = pred$fit,
    se.fit = pred$se.fit
  )
}

# does the prediction
predict_logger(logger_ids[1])


# loops through all loggers
pred_logger <- map_df(logger_ids, predict_logger)




# combine real tomst with predicted tomst -------------------------------------------------


# add column with type (pred or real) -------------------------------------
# predicted tms data
pred_logger <- pred_logger |>
  mutate(logger_type = "tms_predicted")|> 
  rename(temp_mean = fit)

tms_pred_clean <- pred_logger |>
  select(date, site, treat_warming, treat_competition, temp_mean, unique_plot_ID, logger_type)




# real tms data
tms_real <- tomst_nor_daily |> 
  mutate(logger_type = "tms_real")|>
  rename(temp_mean = T3_mean) |> 
  select(date, site, treat_warming, treat_competition, temp_mean, unique_plot_ID, logger_type)


tomst_range_by_treat <- tomst_nor_daily |>
  group_by(site, treat_warming, treat_competition) |>
  summarise(
    start_date = min(date, na.rm = TRUE),
    end_date   = max(date, na.rm = TRUE),
    n_days      = n(),
    .groups = "drop"
  )

tomst_range_by_treat


# join tms predicted and real --------------------------------------------- 
tms_pred_real <- bind_rows(tms_pred_clean, tms_real)

# plot predicted and real -------------------------------------------------
ggplot(tms_pred_real,
       aes(x = date, y = temp_mean, color = logger_type, linetype = site)) +
  geom_line() +
  facet_grid(treat_warming ~ treat_competition)


ggplot(tms_pred_real,
       aes(x = date, y = temp_mean, color = logger_type, linetype = site)) +
  geom_line() +
  facet_grid(treat_warming + site ~ treat_competition)




# make one dataset with predicted first and then real tms --------------------
# assign priority for which to keep

# combine pred_logger with tomst_nor_daily

tms_pred <- pred_logger |>
  mutate(logger_type = "tms_predicted")|> 
  rename(temp_mean_pred = fit) |> 
  select(date, site, treat_warming, treat_competition, unique_plot_ID, temp_mean_pred, logger_type)


tms_real <- tomst_nor_daily |> 
  mutate(logger_type = "tms_real")|>
  rename(temp_mean_real = T3_mean)|> 
  select(date, site, treat_warming, treat_competition, unique_plot_ID, temp_mean_real, logger_type)


tms_joined <- tms_pred |>
  left_join(
    tms_real,
    by = c(
      "date",
      "site",
      "treat_warming",
      "treat_competition",
      "unique_plot_ID"
    )
  )

# make new column chosing pred or real and correct logger type column
tms_joined <- tms_joined |>
  mutate(temp_mean_used = coalesce(temp_mean_real, temp_mean_pred),
    
    logger_type = if_else(
      !is.na(temp_mean_real),
      "tms_real",
      "tms_predicted"))


# select only necessary columns
tms_final <- tms_joined |>
  select(
    date,
    site,
    treat_warming,
    treat_competition,
    unique_plot_ID,
    temp_mean_pred,
    temp_mean_real,
    temp_mean_used,
    logger_type)



ggplot(tms_final,
       aes(x = date, y = temp_mean_used, color = logger_type)) +
  geom_line() +
  facet_grid(treat_warming + site ~ treat_competition)




tms_final |>
  count(
    unique_plot_ID,
    date
  ) |>
  filter(n > 1)

tms_final |>
  count(
    date,
    site,
    treat_warming,
    treat_competition
  ) |>
  arrange(desc(n))


ggplot(tms_joined,
       aes(date,
           temp_mean_used,
           colour = logger_type,
           group = unique_plot_ID)) +
  geom_line()+
  facet_grid(treat_warming + site ~ treat_competition)


tms_treatment_daily <- tms_final |>
  group_by(
    date,
    site,
    treat_warming,
    treat_competition
  ) |>
  summarise(
    temp_mean = mean(temp_mean_used, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(tms_treatment_daily,
       aes(date,
           temp_mean)) +
  geom_line()+
  facet_grid(treat_warming + site ~ treat_competition)


