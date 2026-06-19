


# Figure of DOY vs GDD ----------------------------------------------------
library(colorspace)



source("TMS4_Weatherstation_predictions_NOR.R")

source("Data_preparation_TMS4_CHE.R")






# combine tms data --------------------------------------------------------

# use 
tms_final_nor$region <- "Norway" 
# and
tms_final_che$region <- "Switzerland" 




tms_final_nor_che <- bind_rows(tms_final_nor, tms_final_che)




# combine datasets with GDD -----------------------------------------------
climate_gdd_nor_tms$region <- "Norway" 

climate_gdd_che_tms$region <- "Switzerland" 


gdd_nor_che <- bind_rows(climate_gdd_nor_tms, climate_gdd_che_tms)



# plot both regions all treatments GDD vs DOY -----------------------------
# define colors
define_colors <- c(
  "hi_ambi_vege" = "turquoise4",
  "hi_ambi_bare" = lighten("turquoise4", 0.7),
  "hi_warm_vege" = "darkred",
  "hi_warm_bare" = lighten("darkred", 0.7),
  "lo_ambi_vege" = "grey34",
  "lo_ambi_bare" = lighten("grey34", 0.7))

# plot
ggplot(gdd_nor_che,
       aes(x = jday,
           y = GDD_cum,
           color = treatment_site_temp_comp)) +
  geom_line(size = 1.2) +
  labs(
    x = "Day of year (DOY)",
    y = "Cumulative temperature (GDD2)",
    title = "Norway and Switzerland cumulative GDD")+
  facet_grid(~region)+
  scale_color_manual(values = define_colors)




