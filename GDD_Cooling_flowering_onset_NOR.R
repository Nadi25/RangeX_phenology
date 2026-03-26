
# Growing degree days (GDD) -----------------------------------------------

# RangeX phenology effect of cooling NOR ------------

## Data used: RangeX_clean_phenology_2023_NOR.csv
##            RangeX_clean_MetadataFocal_CHE.csv
##            RangeX_metadata_focal_NOR.csv
##            RangeX_clean_climate_station_NOR_2021-2025.csv
## Date:      08.09.25
## Author:    Nadine Arzt
## Purpose:   Effect of cooling on flowering onset NOR - with GDD


# load library ------------------------------------------------------------
library(conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
library(tidyverse)
library(lubridate)
library(performance)
library(see)


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


# filter only ambi both sites  -------------------------------------
# to compare low ambi with hi ambi = cooling effect
phenology_with_gdd_ambi <- phenology_with_gdd |> 
  filter(treat_warming == "ambi")

# and get julian days --------------------------------------
phenology_with_gdd_ambi <- phenology_with_gdd_ambi |> 
  mutate(jday = yday(date_measurement))   # Julian day (1–365)       



# calculate flowering onset ------------------------------------------------
flowering_onset_gdd_ambi <- phenology_with_gdd_ambi |> 
  filter(phenology_stage == "No_FloOpen", value > 0) |>
  group_by(site, species, unique_plant_ID, block_ID, treat_competition) |>
  summarise(onset = min(GDD_cum, na.rm = TRUE), .groups = "drop") |>
  # remove groups where flowering never occurred
  filter(is.finite(onset))


# model flowering onset nor lmer ------------------------------------------
m_onset_gdd_ambi <- lmerTest::lmer(onset ~ site * treat_competition + (1|species) + (1|block_ID),
                                   data = flowering_onset_gdd_ambi)
summary(m_onset_gdd_ambi)


model_performance(m_onset_gdd_ambi)
check_collinearity(m_onset_gdd_ambi)
check_model(m_onset_gdd_ambi)



# plot point with residuals for effect of warming on flo onset -------------
# estimated marginal means
emm_nor_gdd_ambi <- emmeans(m_onset_gdd_ambi, ~ site | treat_competition)

contr_nor_gdd_ambi <- contrast(emm_nor_gdd_ambi, method = list("hi - lo" = c(1, -1)))

# using summary keeps the p-values
contrast_df_nor_gdd_ambi <- as.data.frame(summary(contr_nor_gdd_ambi, infer = TRUE))


# with raw data points ----------------------------------------------------
# compute mean onset per treatment × group
onset_means_cool_gdd <- flowering_onset_gdd_ambi |>
  group_by(site, species, block_ID, treat_competition) |>
  summarise(mean_onset = mean(onset, na.rm = TRUE), .groups = "drop")

# pivot to get ambi vs warm in same row
delta_onset_cool_gdd <- onset_means_cool_gdd |>
  pivot_wider(names_from = site, values_from = mean_onset) |>
  mutate(delta = hi - lo) |>
  filter(!is.na(delta))

# check result
head(delta_onset_cool_gdd)



# plot --------------------------------------------------------------------
# plot raw deltas + model estimates
nor_che_delta_raw_cool_gdd <- ggplot() +
  # raw deltas (jittered for visibility)
  geom_jitter(data = delta_onset_cool_gdd,
              aes(x = treat_competition, y = delta, color = treat_competition),
              width = 0.1, alpha = 0.4, size = 3) +
  
  # model-based warming effects
  geom_point(data = contrast_df_nor_gdd_ambi, 
             aes(x = treat_competition, y = estimate, color = treat_competition),
             size = 8, position = position_dodge(width = 0.5)) +
  geom_errorbar(data = contrast_df_nor_gdd_ambi,
                aes(x = treat_competition, ymin = lower.CL, ymax = upper.CL, color = treat_competition),
                linewidth = 1, width = 0.1,
                position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  # significance labels
  geom_text(data = contrast_df_nor_gdd_ambi,
            aes(x = treat_competition, y = estimate, 
                label = ifelse(p.value < 0.001, "***",
                               ifelse(p.value < 0.01, "**",
                                      ifelse(p.value < 0.05, "*", "n.s."))),
                color = treat_competition),
            vjust = -1.5, position = position_nudge(x = 0.2),
            show.legend = FALSE, size = 12) +
  
  labs(x = "Biotic interactions",
       y = "Δ GDD shifted flowering onset (high - low)",
       title = "Effect of cooling on flowering onset NOR") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))+
  scale_color_manual(values = c("#528B8B", "#CD950C"))
nor_che_delta_raw_cool_gdd


# ggsave(filename = "Output/Onset/GDD_cooling_flowering_onset_NOR.png", 
#        plot = nor_che_delta_raw_cool_gdd, width = 10, height = 8, units = "in")




