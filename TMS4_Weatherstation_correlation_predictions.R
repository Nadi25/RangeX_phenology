



source("Data_preparation_TMS4_NOR.R")

source("Data_preparation_climate_station_NOR.R")

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



m_tomst <- lm(T3_mean ~ Tavg + site + treat_warming + treat_competition, data = compare_temp2)

summary(m_tomst)
model_performance(m_tomst)

plot(T3_mean ~ predict(m_tomst),
     data = compare_temp2)




dd <- climate_23_daily |> 
  filter(site == "hi")

# predict tomst
newdata <- expand.grid(
  date = dd$date,
  site = c("lo", "hi"),
  treat_warming = c("ambi", "warm"),
  treat_competition = c("vege", "bare")) |>
  filter(!(site == "lo" & treat_warming == "warm"))


newdata <- newdata |>
  left_join(
    climate_23_daily,
    by = c("date", "site")
  )




p <- predict(
  m_tomst,
  newdata = newdata,
  re.form = NA,
  se.fit = TRUE) |> 
  
  as_tibble() |>
  mutate(upper = fit + 1.96 * se.fit,
         lower = fit - 1.96 * se.fit) |>
  bind_cols(newdata)
p

ggplot(p, aes(x = Tavg, y = fit, colour = treat_competition)) +
  geom_point()+
  facet_wrap(~ site)



ggplot(compare_temp2, aes(x = Tavg, y = T3_mean)) +
  geom_point()+
  facet_wrap(~ site)





# one model per treatment -------------------------------------------------

dat_hi_ambi_bare <- compare_temp2 |>
  filter(
    site == "hi",
    treat_warming == "ambi",
    treat_competition == "bare"
  )

m_hi_ambi_bare <- lm(T3_mean ~ Tavg, data = dat_hi_ambi_bare)

summary(m_hi_ambi_bare)
model_performance(m_hi_ambi_bare)


dat_hi_ambi_vege <- compare_temp2 |>
  filter(
    site == "hi",
    treat_warming == "ambi",
    treat_competition == "vege"
  )

m_hi_ambi_vege <- lm(T3_mean ~ Tavg, data = dat_hi_ambi_vege)

summary(m_hi_ambi_vege)
model_performance(m_hi_ambi_vege)


dat_hi_warm_vege <- compare_temp2 |>
  filter(
    site == "hi",
    treat_warming == "warm",
    treat_competition == "vege"
  )

m_hi_warm_vege <- lm(T3_mean ~ Tavg, data = dat_hi_warm_vege)

summary(m_hi_warm_vege)

m_hi_warm_vege2 <- lm(T3_mean ~ Tavg + Radiation + Radiation, data = dat_hi_warm_vege)

summary(m_hi_warm_vege2)

m_hi_warm_vege3 <- lm(T3_mean ~ Tavg + Radiation + WindSpd, data = dat_hi_warm_vege)

summary(m_hi_warm_vege3)


AIC(m_hi_warm_vege, m_hi_warm_vege2, m_hi_warm_vege3)



dat_hi_warm_bare <- compare_temp2 |>
  filter(
    site == "hi",
    treat_warming == "warm",
    treat_competition == "bare"
  )

m_hi_warm_bare <- lm(T3_mean ~ Tavg, data = dat_hi_warm_bare)

summary(m_hi_warm_bare)



dat_lo_ambi_bare <- compare_temp2 |>
  filter(
    site == "lo",
    treat_warming == "ambi",
    treat_competition == "bare"
  )

m_lo_ambi_bare <- lm(T3_mean ~ Tavg, data = dat_lo_ambi_bare)

summary(dat_lo_ambi_bare)


dat_lo_ambi_vege <- compare_temp2 |>
  filter(
    site == "lo",
    treat_warming == "ambi",
    treat_competition == "vege"
  )

m_lo_ambi_vege <- lm(T3_mean ~ Tavg, data = dat_lo_ambi_vege)

summary(dat_lo_ambi_vege)

# predictions
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

