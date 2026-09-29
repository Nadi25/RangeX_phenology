

# Create final tms temp and GDD dataset for both regions --------------------------


# Figure of DOY vs GDD using TMS data ----------------------------------------------------


# load library ------------------------------------------------------------
library(colorspace)



# source data preparation scripts -----------------------------------------
source("TMS4_Weatherstation_predictions_NOR_12h.R")

source("Data_preparation_TMS4_CHE_12h.R")

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



# column with dta type ----------------------------------------------------
# real or predicted
gdd_nor_che <- gdd_nor_che |> 
  mutate(data_type = ifelse(perc_real == 0, "predicted", "measured"))


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
  scale_color_manual(values = define_colors)+
  scale_alpha_continuous(range = c(0.4, 1))
gdd

# ggsave(filename = "Output/Temperature/DOY_GDD_TMS4_NOR_CHE.png", 
#       plot = gdd, width = 15, height = 10, units = "in")



# make dotted line  ----------------------------------
# when at least one of the loggers in the treat are predicted
gdd_nor_che <- gdd_nor_che |> 
  arrange(region, treatment_site_temp_comp, jday) |> 
  group_by(region, treatment_site_temp_comp) |> 
  mutate(
    predicted = perc_real < 100,
    segment = cumsum(predicted != lag(predicted,
                                      default = first(predicted)))
  )

gdd_nor_che <- gdd_nor_che |>
  ungroup()

gdd2 <- ggplot(
  gdd_nor_che,
  aes(
    jday,
    GDD_cum,
    color = treatment_site_temp_comp,
    linetype = predicted,
    group = interaction(treatment_site_temp_comp, segment))) +
  geom_line(linewidth = 1.2) +
  facet_grid(~region) +
  scale_color_manual(values = define_colors)+
  labs(
    x = "Day of year (DOY)",
    y = "Cumulative temperature (GDD5)",
    color = "Treatment",
    title = "Norway and Switzerland cumulative GDD 12h")+
  
  scale_linetype_manual(
    name = "TMS4 data",
    values = c(
      "FALSE" = "solid",
      "TRUE" = "dashed"),
    labels = c(
      "FALSE" = "measured",
      "TRUE" = "(partly) predicted"))
gdd2

# ggsave(filename = "Output/Temperature/DOY_GDD5_TMS4_NOR_CHE_12h.png", 
#       plot = gdd2, width = 15, height = 10, units = "in")





# rename treatments -------------------------------------------------------

gdd_nor_che <- gdd_nor_che |>
  mutate(treatment_site_temp_comp = recode(
      treatment_site_temp_comp,
      "hi_ambi_bare" = "high ambient without",
      "hi_ambi_vege" = "high ambient with",
      "hi_warm_bare" = "high warmed without",
      "hi_warm_vege" = "high warmed with",
      "lo_ambi_bare" = "low ambient without",
      "lo_ambi_vege" = "low ambient with"))

define_colors <- c(
  "high ambient with" = "turquoise4",
  "high ambient without" = lighten("turquoise4", 0.7),
  "high warmed with" = "darkred",
  "high warmed without" = lighten("darkred", 0.7),
  "low ambient with" = "grey34",
  "low ambient without" = lighten("grey34", 0.7))


gdd_nor_che <- gdd_nor_che |>
  mutate(treatment_site_temp_comp = factor(
      treatment_site_temp_comp,
      levels = c(
        "high ambient without",
        "high ambient with",
        "high warmed without",
        "high warmed with",
        "low ambient without",
        "low ambient with")))


theme_set(theme_bw(base_size = 10))

doy_gdd3 <- ggplot(
  gdd_nor_che,
  aes(
    jday,
    GDD_cum,
    color = treatment_site_temp_comp,
    linetype = predicted,
    group = interaction(treatment_site_temp_comp, segment))) +
  geom_line(linewidth = 0.75) +
  facet_grid(~region) +
  scale_color_manual(values = define_colors)+
  labs(
    x = "Day of year (DOY)",
    y = "Cumulative temperature (GDD5)",
    color = "Treatment",
    title = "Norway and Switzerland cumulative GDD 12h")+
  
  scale_linetype_manual(
    name = "Temperature data",
    values = c(
      "FALSE" = "solid",
      "TRUE" = "dotted"),
    labels = c(
      "FALSE" = "measured",
      "TRUE" = "(partly) predicted"))+
  theme(plot.title = element_blank())+
  
  theme(legend.position = c(0.15, 0.75),
        legend.background = element_rect(
          fill = alpha("white", 0.7),
          colour = "NA"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.key.height = unit(0.4, "cm"),
        legend.key.width = unit(0.6, "cm"),
        legend.spacing.y = unit(0.05, "cm"),
        legend.margin = margin(2, 2, 2, 2))
doy_gdd3

# ggsave(filename = "Output/Temperature/Figures_for_manuscript/DOY_GDD5_TMS4_NOR_CHE_12h.png", 
#       plot = doy_gdd3, width = 18, height = 12, units = "cm", dpi = 600)






