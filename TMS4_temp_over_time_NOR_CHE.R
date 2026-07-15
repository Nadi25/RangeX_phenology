

# Create final tms temp and GDD dataset for both regions --------------------------


# Figure of DOY vs GDD using TMS data ----------------------------------------------------


# load library ------------------------------------------------------------
library(colorspace)



# source data preparation scripts -----------------------------------------
source("TMS4_Weatherstation_predictions_NOR.R")

source("Data_preparation_TMS4_CHE.R")

theme_set(theme_bw(base_size = 15))




# combine tms data --------------------------------------------------------

# use 
tms_final_nor$region <- "Norway" 
# and
tms_final_che$region <- "Switzerland" 


tms_final_nor_che <- bind_rows(tms_final_nor, tms_final_che)


tms_final_nor_che_gs <- tms_final_nor_che |>
  filter((region == "Norway" & date >= as.Date("2023-03-15") & date <= as.Date("2023-09-30")) |
           (region == "Switzerland" & date >= as.Date("2022-03-15") & date <= as.Date("2022-09-30")))

tms_final_nor_che_gs <- tms_final_nor_che_gs |> 
  mutate(
    jday = yday(date))   # Julian day (1–365)


define_colors <- c(
  "hi_ambi_vege" = "turquoise4",
  "hi_ambi_bare" = lighten("turquoise4", 0.7),
  "hi_warm_vege" = "darkred",
  "hi_warm_bare" = lighten("darkred", 0.7),
  "lo_ambi_vege" = "grey34",
  "lo_ambi_bare" = lighten("grey34", 0.7))

ggplot(tms_final_nor_che_gs,
       aes(x = jday,
           y = temp_mean,
           color = treatment_site_temp_comp)) +
  geom_line() +
  labs(x = "DOY",
    y = "Daily mean temperature (°C)",
    color = "Treatment")+
  #facet_wrap(~ region, scales = "free_x", ncol = 1)+
  facet_grid(region ~., scales = "free")+
  scale_color_manual(values = define_colors)





# Add season start as lines -----------------------------------------------

season_start_nor_tms <- season_start_nor_tms |>
  mutate(region = "Norway")

season_start_che_tms <- season_start_che_tms |>
  mutate(region = "Switzerland")

season_start_all <- bind_rows(
  season_start_nor_tms,
  season_start_che_tms
) |>
  mutate(jday = yday(season_start))
season_start_all

# add site_warming_comp treatment ----------------------------------------------
season_start_all$treatment_site_temp_comp <- paste(season_start_all$site, season_start_all$treat_warming, 
                                                   season_start_all$treat_competition, sep = "_")




# plot with growing season start ------------------------------------------
temp_tms <- ggplot(tms_final_nor_che_gs,
       aes(x = jday,
           y = temp_mean,
           color = treatment_site_temp_comp)) +
  annotate("rect",
           xmin = 135, # 15.05
           xmax = 244, # 01.09
           ymin = -Inf,
           ymax = Inf,
           fill = "darkseagreen2",
           alpha = 0.2) +
  geom_line() +
  geom_vline(
    data = season_start_all,
    aes(xintercept = jday,
        color = treatment_site_temp_comp),
    linetype = "dashed",
    show.legend = FALSE) +
  geom_hline(yintercept = 0,
             linetype = "dashed")+
  labs(
    x = "DOY",
    y = "Daily mean temperature (°C)",
    color = "Treatment") +
  facet_grid(region ~ ., scales = "free") +
  scale_color_manual(values = define_colors)
temp_tms

# ggsave(filename = "Output/Temperature/Temperature_daily_mean_over_time_NOR_CHE.png", 
#       plot = temp_tms, width = 16, height = 10, units = "in")





# OTC effect: warmed - ambient ------------------------------------------------
tms_final_nor_che_gs



tms_final_nor_che_gs_hi <- tms_final_nor_che_gs |> 
  filter(site == "hi")

tms_final_nor_che_gs_hi <- tms_final_nor_che_gs_hi |> 
  mutate(treat_competition = recode(treat_competition,
                                    "bare" = "without",
                                    "vege" = "with"))


warming_effect <- tms_final_nor_che_gs_hi |>
  select(
    date, jday, region, site,
    treat_warming, treat_competition,
    temp_mean
  ) |>
  pivot_wider(
    names_from = treat_warming,
    values_from = temp_mean
  ) |>
  mutate(temp_diff = warm - ambi)
warming_effect


ggplot(warming_effect,
       aes(x = jday,
           y = temp_diff,
           color = treat_competition)) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "grey50") +
  geom_line(linewidth = 1) +
  facet_grid(region ~ ., scales = "free_x") +
  labs(
    x = "DOY",
    y = expression(Delta * " Temperature (" * degree * "C)"),
    color = "Biotic interactions")+
  scale_color_manual(values = c("with" = "#528B8B", "without" = "#CD950C"))



warming_effect |>
  group_by(region, site, treat_competition) |>
  summarise(
    mean_diff = mean(temp_diff, na.rm = TRUE),
    sd_diff = sd(temp_diff, na.rm = TRUE),
    .groups = "drop"
  )



theme_set(theme_bw(base_size = 20))

ggplot(warming_effect,
       aes(x = treat_competition,
           y = temp_diff,
           fill = treat_competition)) +
  geom_violin(trim = FALSE, 
              alpha = 0.8) +
  geom_boxplot(width = 0.15,
               outlier.shape = NA) +
  facet_grid( ~ region) +
  labs(title = "Temperature difference at high site (T3) with and without OTC",
    x = "Biotic interactions",
    y = "OTC effect (°C)\n(warm - ambient)")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  scale_fill_manual(values = c("with" = "#528B8B", "without" = "#CD950C"))+
  theme(legend.position = "none")




# peak season mean (15.06-15.09) ------------------------------------------

warming_effect_peak <- warming_effect |>
  filter((region == "Norway" & date >= as.Date("2023-06-15") & date <= as.Date("2023-09-15")) |
           (region == "Switzerland" & date >= as.Date("2022-06-15") & date <= as.Date("2022-09-15")))

warming_effect_peak |>
  group_by(region, site, treat_competition) |>
  summarise(
    mean_diff = mean(temp_diff, na.rm = TRUE),
    sd_diff = sd(temp_diff, na.rm = TRUE),
    .groups = "drop")


ggplot(warming_effect_peak,
       aes(x = treat_competition,
           y = temp_diff,
           fill = treat_competition)) +
  geom_violin(trim = FALSE, 
              alpha = 0.7) +
  geom_boxplot(width = 0.15,
               outlier.shape = NA) +
  facet_grid( ~ region) +
  labs(title = "Temperature difference at high site (T3) with and without OTC peak season",
       x = "Biotic interactions",
       y = "Delta temperature (°C)\n(warm - ambient)")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  scale_fill_manual(values = c("with" = "#528B8B", "without" = "#CD950C"))+
  theme(legend.position = "none")

###################

warming_effect_peak2 <- warming_effect |>
  filter((region == "Norway" & date >= as.Date("2023-05-15") & date <= as.Date("2023-09-01")) |
           (region == "Switzerland" & date >= as.Date("2022-05-15") & date <= as.Date("2022-09-01")))


delta_temps <- warming_effect_peak2 |>
  group_by(region, site, treat_competition) |>
  summarise(
    mean_diff = mean(temp_diff, na.rm = TRUE),
    sd_diff = sd(temp_diff, na.rm = TRUE),
    .groups = "drop")
delta_temps


delta_temp_peak <-ggplot(warming_effect_peak2,
       aes(x = treat_competition,
           y = temp_diff,
           fill = treat_competition)) +
  geom_violin(trim = FALSE, 
              alpha = 0.5) +
  geom_boxplot(width = 0.15,
               outlier.shape = NA) +
  facet_grid( ~ region) +
  labs(title = "Temperature difference at high site (T3) with and without OTC peak season",
       x = "Biotic interactions",
       y = "Delta temperature (°C)\n(warm - ambient)")+
  geom_hline(yintercept = 0, linetype = "dashed")+
  scale_fill_manual(values = c("with" = "#528B8B", "without" = "#CD950C"))+
  theme(legend.position = "none")
delta_temp_peak

# ggsave(filename = "Output/Temperature/Temperature_delta_peak_OTC_NOR_CHE.png", 
#       plot = delta_temp_peak, width = 13, height = 10, units = "in")




# Day vs night NOR --------------------------------------------------------

tomst_midday_8_15 <- tomst_23_clean |> 
  mutate(hour = hour(date_time),  # Extract hour from datetime
         day_night = ifelse(hour >= 8 & hour <= 15, "day", "night"))

# filter only day 
tomst_midday_8_15_day <- tomst_midday_8_15 |> 
  filter(day_night == "day")

# filter only night
tomst_midday_8_15_night <- tomst_midday_8_15 |> 
  filter(day_night == "night")





avg_temp_day_long <- tomst_midday_8_15_day |>
  filter(site == "hi") |>
  mutate(date = as.Date(date_time)) |>
  pivot_longer(cols = starts_with("TMS_T"),
               names_to = "sensor",
               values_to = "temperature") |>
  mutate(sensor = recode(sensor,
                         "TMS_T1" = "avg_temp_soil",
                         "TMS_T2" = "avg_temp_surface",
                         "TMS_T3" = "avg_temp_air")) |> 
  group_by(date, treat_warming, sensor) |>
  summarise(mean_temp = mean(temperature, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = treat_warming, values_from = mean_temp) |>
  mutate(delta_temp = warm - ambi) |> 
  mutate(sensor = factor(sensor,
                         levels = c("avg_temp_air", "avg_temp_surface", "avg_temp_soil")))

ggplot(avg_temp_day_long, aes(x = sensor, y = delta_temp, color = sensor)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Sensor", y = "Δ Temperature (warm - ambi)", 
       title = "Daily Temperature Difference 8-15 (High Site)",
       color = "Sensor")





avg_temp_night_long <- tomst_midday_8_15_night |>
  filter(site == "hi") |>
  mutate(date = as.Date(date_time)) |>
  pivot_longer(cols = starts_with("TMS_T"),
               names_to = "sensor",
               values_to = "temperature") |>
  mutate(sensor = recode(sensor,
                         "TMS_T1" = "avg_temp_soil",
                         "TMS_T2" = "avg_temp_surface",
                         "TMS_T3" = "avg_temp_air")) |> 
  group_by(date, treat_warming, sensor) |>
  summarise(mean_temp = mean(temperature, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = treat_warming, values_from = mean_temp) |>
  mutate(delta_temp = warm - ambi) |> 
  mutate(sensor = factor(sensor,
                         levels = c("avg_temp_air", "avg_temp_surface", "avg_temp_soil")))

ggplot(avg_temp_night_long, aes(x = sensor, y = delta_temp, color = sensor)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Sensor", y = "Δ Temperature (warm - ambi)", 
       title = "Daily Temperature Difference night (High Site)",
       color = "Sensor")

# Day vs night CHE --------------------------------------------------------

tomst_midday_8_15_c <- env_che_22 |> 
  mutate(hour = hour(date_time),  # Extract hour from datetime
         day_night = ifelse(hour >= 8 & hour <= 15, "day", "night"))

# filter only day 
tomst_midday_8_15_day_c <- tomst_midday_8_15_c |> 
  filter(day_night == "day")

# filter only night
tomst_midday_8_15_night_c <- tomst_midday_8_15_c |> 
  filter(day_night == "night")





avg_temp_day_long_c <- tomst_midday_8_15_day_c |>
  filter(site == "hi") |>
  mutate(date = as.Date(date_time)) |>
  pivot_longer(cols = starts_with("TMS_T"),
               names_to = "sensor",
               values_to = "temperature") |>
  mutate(sensor = recode(sensor,
                         "TMS_T1" = "avg_temp_soil",
                         "TMS_T2" = "avg_temp_surface",
                         "TMS_T3" = "avg_temp_air")) |> 
  group_by(date, treat_warming, sensor) |>
  summarise(mean_temp = mean(temperature, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = treat_warming, values_from = mean_temp) |>
  mutate(delta_temp = warm - ambi) |> 
  mutate(sensor = factor(sensor,
                         levels = c("avg_temp_air", "avg_temp_surface", "avg_temp_soil")))

ggplot(avg_temp_day_long_c, aes(x = sensor, y = delta_temp, color = sensor)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Sensor", y = "Δ Temperature (warm - ambi)", 
       title = "Daily Temperature Difference 8-15 (High Site)",
       color = "Sensor")





avg_temp_night_long_c <- tomst_midday_8_15_night_c |>
  filter(site == "hi") |>
  mutate(date = as.Date(date_time)) |>
  pivot_longer(cols = starts_with("TMS_T"),
               names_to = "sensor",
               values_to = "temperature") |>
  mutate(sensor = recode(sensor,
                         "TMS_T1" = "avg_temp_soil",
                         "TMS_T2" = "avg_temp_surface",
                         "TMS_T3" = "avg_temp_air")) |> 
  group_by(date, treat_warming, sensor) |>
  summarise(mean_temp = mean(temperature, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = treat_warming, values_from = mean_temp) |>
  mutate(delta_temp = warm - ambi) |> 
  mutate(sensor = factor(sensor,
                         levels = c("avg_temp_air", "avg_temp_surface", "avg_temp_soil")))

ggplot(avg_temp_night_long_c, aes(x = sensor, y = delta_temp, color = sensor)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Sensor", y = "Δ Temperature (warm - ambi)", 
       title = "Daily Temperature Difference night (High Site)",
       color = "Sensor")






# delta temp per stage ambi vs warm ----------------------------------------------------

stage_windows <- read.csv("Data/Stage_windows_NOR_CHE.csv")
stage_windows

stage_windows_long <- stage_windows |>
  select(
    region, site, treat_competition, treat_warming,
    bud_start, bud_end,
    flower_start, flower_end,
    fruit_start, fruit_end,
    seed_start, seed_end
  ) |>
  pivot_longer(
    cols = -(region:treat_warming),
    names_to = c("stage", ".value"),
    names_pattern = "(.*)_(start|end)"
  )

stage_windows_long

stage_windows_long <- stage_windows_long |> 
  mutate(stage = recode(stage,
                        "bud" = "Budding",
         "flower" = "Flowering",
         "fruit" = "Fruiting",
         "seed" = "Seeds"))

tms_final_nor_che

tms_final_nor_che_hi <- tms_final_nor_che |> 
  filter(site == "hi")


tms_final_nor_che_hi <- tms_final_nor_che_hi |> 
  mutate(treat_competition = recode(treat_competition,
                                    "bare" = "without",
                                    "vege" = "with"))



temp_stage <- tms_final_nor_che_hi |>
  left_join(
    stage_windows_long,
    by = c(
      "region",
      "site",
      "treat_competition",
      "treat_warming"
    )
  ) |>
  filter(
    date >= start,
    date <= end
  )
temp_stage



mean_temp_stage <- temp_stage |>
  group_by(
    region,
    site,
    treat_competition,
    treat_warming,
    stage
  ) |>
  summarise(
    mean_temp = mean(temp_mean, na.rm = TRUE),
    sd_temp = sd(temp_mean, na.rm = TRUE),
    n_days = n(),
    .groups = "drop"
  )

mean_temp_stage



delta_temp_stage <- mean_temp_stage |>
  select(
    region,
    site,
    treat_competition,
    treat_warming,
    stage,
    mean_temp
  ) |>
  pivot_wider(
    names_from = treat_warming,
    values_from = mean_temp
  ) |>
  mutate(
    delta_T = warm - ambi
  )

delta_temp_stage

library(gt)

delta_temp_stage |>
  arrange(region, site, treat_competition, stage) |>
  gt() |>
  fmt_number(
    columns = c(ambi, warm, delta_T),
    decimals = 1
  ) |>
  cols_label(
    region = "Region",
    site = "Site",
    treat_competition = "Competition",
    stage = "Phenological stage",
    ambi = "Ambient (°C)",
    warm = "Warming (°C)",
    delta_T = "\u0394T (°C)"
  ) |>
  tab_header(
    title = "Mean temperatures during phenological stages"
  )


delta_stage_plot <- ggplot(
  delta_temp_stage,
  aes(
    x = treat_competition,
    y = delta_T,
    color = treat_competition
  )
) +
  geom_point(
    size = 3,
  ) +
  facet_grid(region
             ~ stage
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_segment(
    aes(
      xend = treat_competition,
      y = 0,
      yend = delta_T
    ),
    colour = "grey70"
  ) +
  scale_color_manual(
    values = c(
      "with" = "#528B8B",
      "without" = "#CD950C"
    )
  ) +
  labs(
    title = "OTC warming effect during each stage",
    x = "Biotic interactions",
    y = "Delta temperature (°C)\n(warm - ambient)"
  ) +
  theme(
    legend.position = "none"
  )

delta_stage_plot





# low high ----------------------------------------------------------------


tms_final_nor_che_ambi <- tms_final_nor_che |> 
  filter(treat_warming == "ambi")


tms_final_nor_che_ambi <- tms_final_nor_che_ambi |> 
  mutate(treat_competition = recode(treat_competition,
                                    "bare" = "without",
                                    "vege" = "with"))



temp_stage_ambi <- tms_final_nor_che_ambi |>
  left_join(
    stage_windows_long,
    by = c(
      "region",
      "site",
      "treat_competition",
      "treat_warming"
    )
  ) |>
  filter(
    date >= start,
    date <= end
  )
temp_stage_ambi



mean_temp_stage_ambi <- temp_stage_ambi |>
  group_by(
    region,
    site,
    treat_competition,
    treat_warming,
    stage
  ) |>
  summarise(
    mean_temp = mean(temp_mean, na.rm = TRUE),
    sd_temp = sd(temp_mean, na.rm = TRUE),
    n_days = n(),
    .groups = "drop"
  )

mean_temp_stage_ambi



delta_temp_stage_ambi <- mean_temp_stage_ambi |>
  select(
    region,
    site,
    treat_competition,
    treat_warming,
    stage,
    mean_temp
  ) |>
  pivot_wider(
    names_from = site,
    values_from = mean_temp
  ) |>
  mutate(
    delta_T = lo - hi
  )

delta_temp_stage_ambi



delta_stage_plot_ambi <- ggplot(
  delta_temp_stage_ambi,
  aes(
    x = treat_competition,
    y = delta_T,
    color = treat_competition
  )
) +
  geom_point(
    size = 3,
  ) +
  facet_grid(region
             ~ stage
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_segment(
    aes(
      xend = treat_competition,
      y = 0,
      yend = delta_T
    ),
    colour = "grey70"
  ) +
  scale_color_manual(
    values = c(
      "with" = "#528B8B",
      "without" = "#CD950C"
    )
  ) +
  labs(
    title = "OTC warming effect during each stage",
    x = "Biotic interactions",
    y = "Delta temperature (°C)\n(warm - ambient)"
  ) +
  theme(
    legend.position = "none"
  )

delta_stage_plot_ambi




delta_temp_stage_ambi$type <- "high low"
delta_temp_stage$type <- "ambient warming"

delta_stage_combined <- bind_rows(delta_temp_stage_ambi, delta_temp_stage)

pd <- position_dodge(width = 0.4)




delta_temp_hl_aw <- ggplot(
  delta_stage_combined,
  aes(
    x = type,
    y = delta_T,
    color = treat_competition,
    shape = type
  )
) +
  geom_point(
    size = 3, 
    position = pd
  ) +
  facet_grid(region
             ~ stage
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  scale_color_manual(
    values = c(
      "with" = "#528B8B",
      "without" = "#CD950C"
    )
  ) +
  labs(
    title = "Temperature difference during each stage",
    x = "Temperature shift",
    y = "Delta temperature (°C)",
    color = "Biotic interactions") +
  theme(
    legend.position = "bottom")+
  guides(shape = "none")
delta_temp_hl_aw


# ggsave(filename = "Output/Temperature/Delta_temperature_stages_NOR_CHE.png", 
#       plot = delta_temp_hl_aw, width = 14, height = 8, units = "in")

# label with temp diff
library(ggrepel)

delta_temp_hl_aw +
  geom_text_repel(
    aes(label = round(delta_T, 1)),
    size = 3,
    show.legend = FALSE,
    position = pd)



t <- bind_rows(
  mutate(delta_temp_stage, type = "Warming"),
  mutate(delta_temp_stage_ambi, type = "Ambient")
) |>
  gt(groupname_col = "type")
t


gtsave(t, "Output/Temperature/Delta_temperature_stages_NOR_CHE_table.docx")





bud_nor <- temp_stage |>
  filter(
    region == "Norway",
    stage == "Budding"
  )

m_bud_nor <- lm(
  temp_mean ~ treat_warming * treat_competition,
  data = bud_nor
)

anova(m_bud_nor)
summary(m_bud_nor)


emm_bud <- emmeans(m_bud_nor,
  ~ treat_warming | treat_competition)

pairs(emm_bud)


flower_nor <- temp_stage |>
  filter(
    region == "Norway",
    stage == "Flowering"
  )

m_flower_nor <- lm(
  temp_mean ~ treat_warming * treat_competition,
  data = flower_nor
)

anova(m_flower_nor)
summary(m_flower_nor)


emm_flower <- emmeans(m_flower_nor,
               ~ treat_warming | treat_competition)

pairs(emm_flower)



fruit_nor <- temp_stage |>
  filter(
    region == "Norway",
    stage == "Fruiting"
  )

m_fruit_nor <- lm(
  temp_mean ~ treat_warming * treat_competition,
  data = fruit_nor
)

anova(m_fruit_nor)
summary(m_fruit_nor)


emm_fruit <- emmeans(m_fruit_nor,
                      ~ treat_warming | treat_competition)

pairs(emm_fruit)





seed_nor <- temp_stage |>
  filter(
    region == "Norway",
    stage == "Seeds"
  )

m_seed_nor <- lm(
  temp_mean ~ treat_warming * treat_competition,
  data = seed_nor
)

anova(m_seed_nor)
summary(m_seed_nor)


emm_seed <- emmeans(m_seed_nor,
                     ~ treat_warming | treat_competition)

pairs(emm_seed)







