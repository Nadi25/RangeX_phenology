


# Figure of DOY vs GDD ----------------------------------------------------
library(colorspace)



source("TMS4_Weatherstation_predictions_NOR.R")

source("Data_preparation_TMS4_CHE.R")

theme_set(theme_bw(base_size = 15))




# combine tms data --------------------------------------------------------

# use 
tms_final_nor$region <- "Norway" 
# and
tms_final_che$region <- "Switzerland" 


tms_final_nor_che <- bind_rows(tms_final_nor, tms_final_che)




# add column region -------------------------------------------------------
climate_gdd_nor_tms$region <- "Norway" 

climate_gdd_che_tms$region <- "Switzerland" 


# filter needed columns nor -----------------------------------------------
gdd_nor_tms <- climate_gdd_nor_tms |> 
  select(-n_pred, -n_real, -n_loggers, - perc_pred)


# add column with percentage of real data ---------------------------------
climate_gdd_che_tms$perc_real <- 100 
gdd_che_tms <- climate_gdd_che_tms


# combine datasets with GDD -----------------------------------------------
gdd_nor_che <- bind_rows(gdd_nor_tms, gdd_che_tms)



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
gdd <- ggplot(gdd_nor_che,
       aes(x = jday,
           y = GDD_cum,
           color = treatment_site_temp_comp,
           alpha = perc_real)) +
  geom_line(size = 1.2) +
  labs(
    x = "Day of year (DOY)",
    y = "Cumulative temperature (GDD2)",
    alpha = "% measured\nTMS4 data",
    color = "Treatment",
    title = "Norway and Switzerland cumulative GDD")+
  facet_grid(~region)+
  scale_color_manual(values = define_colors)
gdd

# ggsave(filename = "Output/Temperature/DOY_GDD_TMS4_NOR_CHE.png", 
#       plot = gdd, width = 15, height = 10, units = "in")



