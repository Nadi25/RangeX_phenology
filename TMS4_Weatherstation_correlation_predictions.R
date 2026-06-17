

library(conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
library(tidyverse)
library(performance)

source("Data_preparation_TMS4_NOR.R")

source("Data_preparation_climate_station_NOR.R")

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
    
    Humidity = mean(Humidity_Avg, na.rm = TRUE),
    WindSpd = mean(WindSpd_Avg, na.rm = TRUE),
    Radiation = mean(Radiation_Avg, na.rm = TRUE),
    
    Rainfall = sum(Rainfall, na.rm = TRUE),
    
    .groups = "drop"
  ) |> 
  mutate(
    Tavg = (Tmax + Tmin) / 2
  )

# filter only year 2023 -----------------------------------------------------
# has end of 2022 in 
climate_23_daily <- climate_23_daily |> 
  filter(date >= as.Date("2023-01-01"))



# tomst -------------------------------------------------------------------


# calculate Temp average per day -----------------
env_nor_daily <- env_nor |>
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



# one value per day and plot

range(env_nor_daily$date)

range(climate_23_daily$date)





# with treatments
compare_temp2 <- env_nor_daily |>
  inner_join(climate_23_daily, by = c("date", "site"))



# m_tomst <- lm(T3_mean ~ Tavg + site + treat_warming + treat_competition, data = compare_temp2)
# 
# summary(m_tomst)
# model_performance(m_tomst)
# 
# plot(T3_mean ~ predict(m_tomst),
#      data = compare_temp2)
# 
# 

# 
# dd <- climate_23_daily |> 
#   filter(site == "hi")
# 
# # predict tomst
# newdata <- expand.grid(
#   date = dd$date,
#   site = c("lo", "hi"),
#   treat_warming = c("ambi", "warm"),
#   treat_competition = c("vege", "bare")) |>
#   filter(!(site == "lo" & treat_warming == "warm"))
# 
# 
# newdata <- newdata |>
#   left_join(
#     climate_23_daily,
#     by = c("date", "site")
#   )
# 
# p <- predict(
#   m_tomst,
#   newdata = newdata,
#   re.form = NA,
#   se.fit = TRUE) |> 
#   
#   as_tibble() |>
#   mutate(upper = fit + 1.96 * se.fit,
#          lower = fit - 1.96 * se.fit) |>
#   bind_cols(newdata)
# p

# ggplot(p, aes(x = Tavg, y = fit, colour = treat_competition)) +
#   geom_point()+
#   facet_wrap(~ site)



ggplot(compare_temp2, aes(x = Tavg, y = T3_mean)) +
  geom_point()+
  facet_wrap(~ site)





# one model per treatment -------------------------------------------------

# hi ambi bare
dat_hi_ambi_bare <- compare_temp2 |>
  filter(
    site == "hi",
    treat_warming == "ambi",
    treat_competition == "bare"
  )

m_hi_ambi_bare <- lm(T3_mean ~ Tavg, data = dat_hi_ambi_bare)

summary(m_hi_ambi_bare)
model_performance(m_hi_ambi_bare)

m_hi_ambi_bare2 <- lm(T3_mean ~ Tavg + Radiation, data = dat_hi_ambi_bare)

summary(m_hi_ambi_bare2)
model_performance(m_hi_ambi_bare2)

m_hi_ambi_bare3 <- lm(T3_mean ~ Tavg + Radiation + WindSpd, data = dat_hi_ambi_bare)

summary(m_hi_ambi_bare3)
model_performance(m_hi_ambi_bare3)

AIC(m_hi_ambi_bare, m_hi_ambi_bare2, m_hi_ambi_bare3)


# hi ambi vege
dat_hi_ambi_vege <- compare_temp2 |>
  filter(
    site == "hi",
    treat_warming == "ambi",
    treat_competition == "vege"
  )

m_hi_ambi_vege <- lm(T3_mean ~ Tavg, data = dat_hi_ambi_vege)

summary(m_hi_ambi_vege)
model_performance(m_hi_ambi_vege)

m_hi_ambi_vege2 <- lm(T3_mean ~ Tavg + Radiation, data = dat_hi_ambi_vege)

summary(m_hi_ambi_vege2)
model_performance(m_hi_ambi_vege2)

m_hi_ambi_vege3 <- lm(T3_mean ~ Tavg + Radiation + WindSpd, data = dat_hi_ambi_vege)

summary(m_hi_ambi_vege3)
model_performance(m_hi_ambi_vege3)

AIC(m_hi_ambi_vege, m_hi_ambi_vege2, m_hi_ambi_vege3)




# hi warm vege
dat_hi_warm_vege <- compare_temp2 |>
  filter(
    site == "hi",
    treat_warming == "warm",
    treat_competition == "vege"
  )

m_hi_warm_vege <- lm(T3_mean ~ Tavg, data = dat_hi_warm_vege)

summary(m_hi_warm_vege)
model_performance(m_hi_warm_vege)


m_hi_warm_vege2 <- lm(T3_mean ~ Tavg + Radiation, data = dat_hi_warm_vege)

summary(m_hi_warm_vege2)
model_performance(m_hi_warm_vege2)

m_hi_warm_vege3 <- lm(T3_mean ~ Tavg + Radiation + WindSpd, data = dat_hi_warm_vege)

summary(m_hi_warm_vege3)
model_performance(m_hi_warm_vege3)


AIC(m_hi_warm_vege, m_hi_warm_vege2, m_hi_warm_vege3)


# hi warm bare
dat_hi_warm_bare <- compare_temp2 |>
  filter(
    site == "hi",
    treat_warming == "warm",
    treat_competition == "bare"
  )

m_hi_warm_bare <- lm(T3_mean ~ Tavg, data = dat_hi_warm_bare)

summary(m_hi_warm_bare)
model_performance(m_hi_warm_bare)


m_hi_warm_bare2 <- lm(T3_mean ~ Tavg + Radiation, data = dat_hi_warm_bare)

summary(m_hi_warm_bare2)
model_performance(m_hi_warm_bare2)

m_hi_warm_bare3 <- lm(T3_mean ~ Tavg + Radiation + WindSpd, data = dat_hi_warm_bare)

summary(m_hi_warm_bare3)
model_performance(m_hi_warm_bare3)

AIC(m_hi_warm_bare, m_hi_warm_bare2, m_hi_warm_bare3)


# lo ambi bare
dat_lo_ambi_bare <- compare_temp2 |>
  filter(
    site == "lo",
    treat_warming == "ambi",
    treat_competition == "bare"
  )

m_lo_ambi_bare <- lm(T3_mean ~ Tavg, data = dat_lo_ambi_bare)

summary(m_lo_ambi_bare)

m_lo_ambi_bare2 <- lm(T3_mean ~ Tavg + Radiation , data = dat_lo_ambi_bare)

m_lo_ambi_bare3 <- lm(T3_mean ~ Tavg + Radiation + WindSpd, data = dat_lo_ambi_bare)


AIC(m_lo_ambi_bare, m_lo_ambi_bare2, m_lo_ambi_bare3)
model_performance(m_lo_ambi_bare)
model_performance(m_lo_ambi_bare3)


# lo ambi vege
dat_lo_ambi_vege <- compare_temp2 |>
  filter(
    site == "lo",
    treat_warming == "ambi",
    treat_competition == "vege"
  )

m_lo_ambi_vege <- lm(T3_mean ~ Tavg, data = dat_lo_ambi_vege)

summary(m_lo_ambi_vege)

m_lo_ambi_vege2 <- lm(T3_mean ~ Tavg + Radiation, data = dat_lo_ambi_vege)

m_lo_ambi_vege3 <- lm(T3_mean ~ Tavg + Radiation + WindSpd, data = dat_lo_ambi_vege)

AIC(m_lo_ambi_vege, m_lo_ambi_vege2, m_lo_ambi_vege3)
model_performance(m_lo_ambi_vege2)
model_performance(m_lo_ambi_vege3)


# check  correlation between T and Radiation
cor(compare_temp2$Tavg,
    compare_temp2$Radiation,
    use = "complete.obs")

cor(compare_temp2$Tavg,
    compare_temp2$WindSpd,
    use = "complete.obs")




# predictions -------------------------------------------------------------

new_hi_ambi_vege <- climate_23_daily |>
  filter(site == "hi")

pred_hi_ambi_vege <- predict(
  m_hi_ambi_vege,
  newdata = new_hi_ambi_vege,
  se.fit = TRUE
) |>
  as_tibble() |>
  mutate(
    upper = fit + 1.96 * se.fit,
    lower = fit - 1.96 * se.fit
  ) |>
  bind_cols(new_hi_ambi_vege) |>
  mutate(
    treat_warming = "ambi",
    treat_competition = "vege"
  )
pred_hi_ambi_vege





predict_tomst <- function(model, data, warming, competition) {
  
  predict(
    model,
    newdata = data,
    se.fit = TRUE
  ) |>
    as_tibble() |>
    mutate(
      upper = fit + 1.96 * se.fit,
      lower = fit - 1.96 * se.fit,
      treat_warming = warming,
      treat_competition = competition
    ) |>
    bind_cols(data)
}

new_hi <- climate_23_daily |> filter(site == "hi")
new_lo <- climate_23_daily |> filter(site == "lo")



pred_hi_ambi_vege <- predict_tomst(
  m_hi_ambi_vege3, new_hi, "ambi", "vege"
)

pred_hi_ambi_bare <- predict_tomst(
  m_hi_ambi_bare3, new_hi, "ambi", "bare"
)

pred_hi_warm_vege <- predict_tomst(
  m_hi_warm_vege3, new_hi, "warm", "vege"
)

pred_hi_warm_bare <- predict_tomst(
  m_hi_warm_bare3, new_hi, "warm", "bare"
)



pred_lo_ambi_vege <- predict_tomst(
  m_lo_ambi_vege3, new_lo, "ambi", "vege"
)

pred_lo_ambi_bare <- predict_tomst(
  m_lo_ambi_bare3, new_lo, "ambi", "bare"
)



pred_all <- bind_rows(
  pred_hi_ambi_vege,
  pred_hi_ambi_bare,
  pred_hi_warm_vege,
  pred_hi_warm_bare,
  pred_lo_ambi_vege,
  pred_lo_ambi_bare
)
pred_all

pred_all2 <- pred_all |>
  filter(date >= as.Date("2023-03-01"),
         date <= as.Date("2023-09-30"))

ggplot(pred_all2,
       aes(date, fit,
           colour = treat_competition,
           linetype = treat_warming)) +
  geom_line() +
  facet_wrap(~site)



# compare real with predicted tms data ------------------------------------
obs <- compare_temp2 |>
  select(date, site, treat_warming, treat_competition, T3_mean)


pred <- pred_all |>
  select(date, site, treat_warming, treat_competition, fit)


validation <- obs |>
  inner_join(pred,
             by = c("date", "site", "treat_warming", "treat_competition"))

ggplot(validation,
       aes(x = T3_mean, y = fit,
           colour = treat_competition)) +
  geom_point(alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(~site)

model_performance(lm(T3_mean ~ fit, data = validation))


#####

# combine real tomst with predicted tomst -------------------------------------------------

# rename column of pred values
pred_all2 <- pred_all2  |> 
  rename(tomst_pred = fit)



# add column with type (pred or real) -------------------------------------
# predicted tms data
pred_all2 <- pred_all2 |>
  mutate(logger_type = "tms_predicted")

tms_pred_clean <- pred_all2 |>
  select(date, site, treat_warming, treat_competition, temp_mean = tomst_pred, logger_type)


# real tms data
tms_real_clean <- env_nor_daily |>
  group_by(date, site, treat_warming, treat_competition) |>
  summarise(temp_mean = mean(T3_mean, na.rm = TRUE),
    .groups = "drop") |>
  mutate(logger_type = "tms_real")




# join tms predicted and real ---------------------------------------------
tms_pred_real <- bind_rows(tms_pred_clean, tms_real_clean)



# plot predicted and real -------------------------------------------------

ggplot(tms_pred_real,
       aes(x = date, y = temp_mean, color = logger_type)) +
  geom_line() +
  facet_grid(treat_warming + site ~ treat_competition)


ggplot(tms_pred_real,
       aes(x = date, y = temp_mean, color = logger_type, linetype = site)) +
  geom_line() +
  facet_grid(treat_warming  ~ treat_competition)



# add climate station real data to tms ------------------------------------
climate_real <- climate_23_daily |>
  mutate(
    temp_mean = Tavg,
    logger_type = "climate_station",
    treat_warming = NA_character_,
    treat_competition = NA_character_
  ) |>
  select(date, site, treat_warming, treat_competition,
         temp_mean, logger_type)






# add climate station to pred and real tms data ---------------------------
tms_climate_all <- bind_rows(
  climate_real,
  tms_pred_real
)

# plot everything
ggplot(tms_climate_all,
       aes(x = date, y = temp_mean, color = logger_type)) +
  geom_line() +
  facet_grid(treat_warming ~ site + treat_competition)



ggplot(tms_climate_all,
       aes(x = date, y = temp_mean, color = logger_type, linetype = site)) +
  geom_line() 
  #facet_grid(treat_warming ~ site + treat_competition)






# one model per logger ----------------------------------------------------
compare_temp2 |>
  distinct(
    unique_plot_ID,
    site,
    treat_warming,
    treat_competition
  ) |>
  arrange(unique_plot_ID)


logger_models <- compare_temp2 |>
  group_by(unique_plot_ID) |>
  group_split() |>
  set_names(
    compare_temp2 |>
      distinct(unique_plot_ID) |>
      pull(unique_plot_ID)
  ) |>
  map(~ lm(T3_mean ~ Tavg + Radiation + WindSpd, data = .x))


names(logger_models)

summary(logger_models[[1]])

logger_data <- compare_temp2 |>
  filter(unique_plot_ID == names(logger_models)[1])

head(logger_data)







# test out workflow with one logger as example ----------------------------

logger_ids <- unique(compare_temp2$unique_plot_ID)

logger_ids


id <- logger_ids[1]

dat_logger <- compare_temp2 |>
  filter(unique_plot_ID == id)

m_logger <- lm(T3_mean ~ Tavg + Radiation + WindSpd, data = dat_logger)

summary(m_logger)


unique(dat_logger$site)


newdata <- climate_23_daily |>
  filter(site == "hi")

pred <- predict(
  m_logger,
  newdata = newdata,
  se.fit = TRUE
)

pred_logger <- tibble(
  unique_plot_ID = id,
  date = newdata$date,
  fit = pred$fit,
  se.fit = pred$se.fit
)
pred_logger





# predict per logger and loop through all ---------------------------------

logger_ids <- unique(compare_temp2$unique_plot_ID)

logger_ids
length(logger_ids)


predict_logger <- function(id){
  
  # data from one logger
  dat_logger <- compare_temp2 |>
    filter(unique_plot_ID == id)
  
  # model
  m_logger <- lm(
    T3_mean ~ Tavg + Radiation + WindSpd,
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

predict_logger(logger_ids[1])



pred_logger <- map_df(logger_ids, predict_logger)












