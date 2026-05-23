


# GDD vs DOY --------------------------------------------------------------

# plot gdd against calender day to see how warth is accumumlating


source("Data_preparation_climate_station_CHE.R")
source("Data_preparation_climate_station_NOR.R")

theme_set(theme_bw(base_size = 20))


climate_gdd_che$type <- "Switzerland climate station"
climate_gdd$type <- "Norway climate station"

climate_gdd_che$region <- "Switzerland"
climate_gdd$region <- "Norway"

climate_gdd_nor <- climate_gdd


climate_all <- bind_rows(climate_gdd_che, climate_gdd_nor)


climate_all <- climate_all |>
  mutate(jday = lubridate::yday(date_measurement))



cum_gdd_doy <- ggplot(climate_all,
       aes(x = jday, y = GDD_cum,
           color = region, linetype = site)) +
  
  geom_line(linewidth = 1.2) +
  
  labs(x = "Day of year",
       y = "Cumulative temperature in GDD",
       color = "Site",
       linetype = "Site") +
  
  scale_color_manual(values = c("turquoise4", "pink4"))
cum_gdd_doy


# ggsave(filename = "Output/Onset/GDD_DOY_comparison_NOR_CHE_baseT2.png", 
#        plot = cum_gdd_doy, width = 12, height = 10, units = "in")





# add tomst CHE -----------------------------------------------------------

source("Data_preparation_TMS4_CHE.R")


climate_gdd_che_tms$type <- "Switzerland TMS4"
climate_gdd_che_tms$region <- "Switzerland"




climate_all2 <- bind_rows(climate_gdd_che, climate_gdd_nor, climate_gdd_che_tms)


climate_all2 <- climate_all2 |>
  mutate(jday = lubridate::yday(date_measurement))



cum_gdd_doy2 <- ggplot(climate_all2,
                      aes(x = jday, y = GDD_cum,
                          color = region, linetype = site)) +
  
  geom_line(linewidth = 1.2) +
  
  labs(x = "Day of year",
       y = "Cumulative temperature in GDD",
       color = "Site",
       linetype = "Site") +
  
  scale_color_manual(values = c("turquoise4", "pink4", "blue"))
cum_gdd_doy2



# compare same dates
compare_temp <- climate_tomst_22_site |>
  left_join(climate_22_daily,
            by = c("site", "date")) |>
  mutate(
    diff = Tavg_tms - Tavg
  )

summary(compare_temp$diff)


cor(compare_temp$Tavg_tms,
    compare_temp$Tavg,
    use = "complete.obs")


ggplot(compare_temp,
       aes(x = Tavg, y = Tavg_tms)) +
  geom_point(alpha = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")




# compare daily temps cliamte station and tomst CHE -----------------------



# join datasets
compare_temp2 <- climate_tomst_22_site |>
  left_join(
    climate_22_daily,
    by = c("site", "date")
  )

compare_temp2 <- compare_temp2 |> 
  filter(date >= as.Date("2022-04-01"), 
         date <= as.Date("2022-10-01"))

# plot
ggplot(compare_temp2,
       aes(x = date)) +
  
  # climate station
  geom_line(aes(y = Tavg,
                color = "Climate station"),
            linewidth = 1) +
  
  # TOMST
  geom_line(aes(y = Tavg_tms,
                color = "TOMST"),
            linewidth = 0.7,
            alpha = 0.7) +
  
  facet_wrap(~ site, scales = "free_x") +
  
  scale_color_manual(values = c(
    "Climate station" = "black",
    "TOMST" = "red"
  )) +
  
  labs(
    x = "Date",
    y = "Daily mean temperature (°C)",
    color = "",
    title = "Comparison of daily temperatures: TOMST vs climate station"
  ) 



# add combined tms climate station  ---------------------------------------

climate_gdd_che_comb$type <- "Switzerland comb"
climate_gdd_che_comb$region <- "Switzerland comb" # if we keep that we can have three panels, 2 for che


climate_all3 <- bind_rows(climate_gdd_che, climate_gdd_nor, climate_gdd_che_tms, climate_gdd_che_comb)


climate_all3 <- climate_all3 |>
  mutate(jday = lubridate::yday(date_measurement))



cum_gdd_doy3 <- ggplot(climate_all3,
                       aes(x = jday, y = GDD_cum,
                           color = region, linetype = site)) +
  
  geom_line(linewidth = 1.2) +
  
  labs(x = "Day of year",
       y = "Cumulative temperature in GDD",
       color = "Site",
       linetype = "Site") +
  
  scale_color_manual(values = c("turquoise4", "pink4", "red", "blue"))
cum_gdd_doy3

# ggsave(filename = "Output/Onset/GDD_DOY_Tb2_climate_tms_NOR_CHE.png", 
#       plot = cum_gdd_doy3,
#        width = 15, height = 10, units = "in")





cum_gdd_doy4 <- ggplot(climate_all3,
                       aes(x = jday, y = GDD_cum,
                           color = region, linetype = site)) +
  
  geom_line(linewidth = 1.2) +
  facet_grid(~ type)+
  
  labs(x = "Day of year",
       y = "Cumulative temperature in GDD",
       color = "Site",
       linetype = "Site") +
  
  scale_color_manual(values = c("turquoise4", "pink4", "red", "blue"))
cum_gdd_doy4

# ggsave(filename = "Output/Onset/GDD_DOY_Tb2_climate_tms_NOR_CHE.png", 
#       plot = cum_gdd_doy3,
#        width = 15, height = 10, units = "in")


ggplot(climate_all3,
       aes(x = jday, y = GDD_cum,
           color = type, linetype = site)) +
  
  geom_line(linewidth = 1.2) +
  facet_grid(~ region)+
  
  labs(x = "Day of year",
       y = "Cumulative temperature in GDD",
       color = "Site",
       linetype = "Site") +
  
  scale_color_manual(values = c("turquoise4", "pink4", "red", "blue"))


















