
# Growing degree days (GDD) -----------------------------------------------

# RangeX phenology effect of cooling NOR CHE ------------


# load library ------------------------------------------------------------
library(conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
library(tidyverse)
library(lubridate)
library(performance)
library(see)
library(emmeans)

theme_set(theme_bw(base_size = 22))

source("Data_preparation_phenology_NOR_CHE_combined.R")


# filter only nor ------------------------------------------------------
phenology_nor <- phenology |> 
  filter(region == "Norway")


unique(phenology_nor$date_measurement)
# start: 2023-05-12
# end: 2023-09-13


source("Data_preparation_climate_station_NOR.R")

# use 
climate_gdd_pt

# combine gdd_cum and phenology_nor ------------------------------------
phenology_with_gdd <- phenology_nor |> 
  left_join(climate_gdd_pt |> 
              select(site, date_measurement, GDD_cum),
            by = c("site", "date_measurement"))



# filter only che ------------------------------------------------------
phenology_che <- phenology |> 
  filter(region == "Switzerland")


unique(phenology_che$date_measurement)
# start: "2022-05-04"
# end: "2022-09-27"


source("Data_preparation_climate_station_CHE.R")

# use 
climate_gdd_che

# combine gdd_cum and phenology_nor ------------------------------------
phenology_with_gdd_che <- phenology_che |> 
  left_join(climate_gdd_che |> 
              select(site, date_measurement, GDD_cum),
            by = c("site", "date_measurement"))




# combine nor and che -----------------------------------------------------
phenology_gdd_nor_che <- rbind(phenology_with_gdd, phenology_with_gdd_che)



# filter only ambi both sites  -------------------------------------
# to compare low ambi with hi ambi = cooling effect
phenology_gdd_nor_che <- phenology_gdd_nor_che |> 
  filter(treat_warming == "ambi")

# and get julian days --------------------------------------
phenology_gdd_nor_che <- phenology_gdd_nor_che |> 
  mutate(jday = yday(date_measurement))   # Julian day (1–365)       



# calculate flowering onset ------------------------------------------------
flowering_onset_gdd_ambi <- phenology_gdd_nor_che |> 
  filter(phenology_stage == "No_FloOpen", value > 0) |>
  group_by(region, site, species, unique_plant_ID, block_ID, treat_competition) |>
  summarise(onset = min(GDD_cum, na.rm = TRUE), .groups = "drop") |>
  # remove groups where flowering never occurred
  filter(is.finite(onset))



# model flowering onset nor che ------------------------------------------
m_onset_gdd_ambi_flower <- lmerTest::lmer(onset ~ region * site * treat_competition + 
                                            (1|species) + (1|block_ID),
                                   data = flowering_onset_gdd_ambi)
summary(m_onset_gdd_ambi_flower)


model_performance(m_onset_gdd_ambi_flower)
check_collinearity(m_onset_gdd_ambi_flower)
check_model(m_onset_gdd_ambi_flower)


# Plot with interactions --------------------------------------------------

emm_full_gdd_flower <- emmeans(
  m_onset_gdd_ambi_flower,
  ~ site * treat_competition | region)


contr_cooling_gdd_flower <- contrast(
  emm_full_gdd_flower,
  method = list("hi - lo" = c(1, -1)),
  by = c("region", "treat_competition"),
  simple = "site"
)


df_cooling_gdd_flower <- as.data.frame(summary(contr_cooling_gdd_flower, infer = TRUE)) |>
  mutate(effect = "Transplantation",
         group = treat_competition)


contr_comp_gdd_flower <- contrast(
  emm_full_gdd_flower,
  method = list("with - without" = c(1, -1)),
  by = c("region", "site"),
  simple = "treat_competition"
)


df_comp_gdd_flower <- as.data.frame(summary(contr_comp_gdd_flower, infer = TRUE)) |>
  mutate(effect = "Biotic interactions",
         group = site)



contr_interaction_gdd_flower <- contrast(
  emm_full_gdd_flower,
  interaction = "pairwise"
)

df_interaction_gdd_flower <- as.data.frame(summary(contr_interaction_gdd_flower, infer = TRUE)) |>
  mutate(effect = "Interaction",
         group = "all interactions")

#

add_stars <- function(df){
  df |>
    mutate(
      stars = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        TRUE ~ "n.s."
      )
    )
}

df_cooling_gdd_flower     <- add_stars(df_cooling_gdd_flower)
df_comp_gdd_flower        <- add_stars(df_comp_gdd_flower)
df_interaction_gdd_flower <- add_stars(df_interaction_gdd_flower)

plot_df_gdd_flower <- bind_rows(
  df_cooling_gdd_flower     |> select(region, effect, group, estimate, lower.CL, upper.CL, stars),
  df_comp_gdd_flower        |> select(region, effect, group, estimate, lower.CL, upper.CL, stars),
  df_interaction_gdd_flower |> select(region, effect, group, estimate, lower.CL, upper.CL, stars)
)

plot_df_gdd_flower <- plot_df_gdd_flower |>
  mutate(effect = factor(effect,
                         levels = c("Transplantation",
                                    "Biotic interactions",
                                    "Interaction")))
#
delta_raw_gdd_flower <- flowering_onset_gdd_ambi |>
  pivot_wider(names_from = c(site, treat_competition),
              values_from = onset) |>
  mutate(
    cooling_with    = hi_with - lo_with,
    cooling_without = hi_without - lo_without,
    interaction     = cooling_with - cooling_without
  )


cooling_flower_interaction_gdd_flower_nor_che <- ggplot(plot_df_gdd_flower,
                                     aes(x = group, y = estimate, color = group)) +
  
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0.1) +
  
  facet_grid(effect ~ region, scales = "free_y") +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  geom_text(aes(y = upper.CL, label = stars),
            vjust = 0.1,
            size = 8,
            position = position_dodge(width = 0.5),
            show.legend = FALSE) +
  
  theme(
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20)
  ) +
  
  labs(x = NULL,
       y = "Shifted onset (GDD)",
       color = "Condition")
cooling_flower_interaction_gdd_flower_nor_che



# ggsave(filename = "Output/Onset/GDD_cooling_flowering_onset_NOR_CHE.png", 
#        plot = cooling_flower_interaction_gdd_flower_nor_che, 
#        width = 15, height = 10, units = "in")





