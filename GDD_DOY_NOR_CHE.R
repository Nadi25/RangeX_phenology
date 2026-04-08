


# GDD vs DOY --------------------------------------------------------------

# plot gdd against calender day to see how warth is accumumlating


source("Data_preparation_climate_station_CHE.R")
source("Data_preparation_climate_station_NOR.R")

theme_set(theme_bw(base_size = 20))


climate_gdd_che$region <- "Switzerland"
climate_gdd_pt$region <- "Norway"

climate_gdd_nor <- climate_gdd_pt


climate_all <- bind_rows(climate_gdd_che, climate_gdd_nor)


climate_all <- climate_all |>
  mutate(jday = lubridate::yday(date_measurement))



cum_gdd_doy <- ggplot(climate_all,
       aes(x = jday, y = GDD_cum,
           color = site, linetype = region)) +
  
  geom_line(linewidth = 1.2) +
  
  labs(x = "Day of year",
       y = "Cumulative temperature in GDD",
       color = "Site",
       linetype = "Site") +
  
  scale_color_manual(values = c("blue", "pink2"))
cum_gdd_doy


# ggsave(filename = "Output/Onset/GDD_DOY_comparison_NOR_CHE.png", 
#        plot = cum_gdd_doy, width = 12, height = 10, units = "in")





























