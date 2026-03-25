

# Effect of transplantation on Fruiting onset with and without competition ---------------------------------------

## Data used: RangeX_clean_Phenology_2022_CHE.csv
##            RangeX_clean_phenology_2023_NOR.csv
##            RangeX_clean_MetadataFocal_CHE.csv
##            RangeX_metadata_focal_NOR.csv
## Date:      02.03.26
## Author:    Nadine Arzt
## Purpose:   Effect of transplantation on Fruiting onset NOR and CHE


# load library ---------------------------------------------------------
library(lme4)
library(ggeffects)
library(broom.mixed)
library(emmeans)
library(lubridate)
library(performance)
library(see)


# load clean phenology data -----------------------------------------------
source("Data_preparation_phenology_NOR_CHE_combined.R")

# use this data set
names(phenology)

# set theme for plots for presentation ------------------------------------
theme_set(theme_bw(base_size = 20))



# combine the fruiting stages nor and che to comapre the onset ------------
# this is to be taken with caution because the stages are not the same
# but for the onset it could be comparable
phenology <- phenology |>
  mutate(phenology_stage = recode(
    phenology_stage,
    "No_Infructescences" = "No_FloWithrd"
  ))


# filter only hi ambi and lo -----------------------------------------------
# should maybe do that in the model

# and get julian days ---------------------------------------------------
# yday(date)
# che and nor was measured in two years but if we count the days in each year it should be fine

phenology_cool <- phenology |> 
  filter(treat_warming == "ambi") |> 
  mutate(
    jday = yday(date_measurement),   # Julian day (1–365)
    jday_scaled = scale(jday))        # optional scaling if needed


# calculate fruiting onset ------------------------------------------------
onset_fruit_cool <- phenology_cool |> 
  filter(phenology_stage == "No_FloWithrd", value > 0) |>
  group_by(region, site, treat_competition, species, block_ID, unique_plot_ID, 
           unique_plant_ID, phenology_stage) |>
  summarise(onset = min(jday, na.rm = TRUE), .groups = "drop") |>
  # remove groups where fruits never occurred
  filter(is.finite(onset))

# model with region for fruiting onset lmer ----------------------------------
m_onset_fruit_cooling <- lmerTest::lmer(onset ~ region * site * treat_competition + 
                                        (1|species) + (1|block_ID), 
                                      data = onset_fruit_cool)

summary(m_onset_fruit_cooling)

model_performance(m_onset_fruit_cooling)
check_collinearity(m_onset_fruit_cooling)
check_model(m_onset_fruit_cooling)


# emmeans -----------------------------------------------------------------
# get emmeans for warming within each region × competition
# this calculates marginal means per site
emm_fruit_cool <- emmeans(m_onset_fruit_cooling, ~ site | region * treat_competition)



# contrasts high-low ------------------------------------------------------
# compute contrasts (high - low) within each competition level

contr_fruit_cool <- contrast(
  emm_fruit_cool,
  method = list("hi - lo" = c(1, -1)),
  by = c("region", "treat_competition")
)



# using summary keeps the p-values
contrast_df_fruit_cool <- as.data.frame(summary(contr_fruit_cool, infer = TRUE))



# plot with raw data points ----------------------------------------------------
# compute mean onset per treatment × group
onset_fruit_means_cool <- onset_fruit_cool |>
  group_by(region, site, treat_competition, species, block_ID) |>
  summarise(mean_onset = mean(onset, na.rm = TRUE), .groups = "drop")

# pivot to get ambi vs warm in same row
delta_onset_fruit_cool <- onset_fruit_means_cool |>
  pivot_wider(names_from = site, values_from = mean_onset) |>
  mutate(delta = hi - lo) |>
  filter(!is.na(delta))

# check result
head(delta_onset_fruit_cool)



# plot raw deltas + model estimates ---------------------------------------

nor_che_delta_raw_cool_fruit <- ggplot() +
  # raw deltas
  geom_jitter(
    data = delta_onset_fruit_cool,
    aes(
      x = region,
      y = delta,
      color = treat_competition
    ),
    width = 0.1, alpha = 0.4, size = 3
  ) +
  
  # model-based contrasts
  geom_point(
    data = contrast_df_fruit_cool,
    aes(
      x = region,
      y = estimate,
      color = treat_competition,
      shape = region
    ),
    size = 8,
    position = position_dodge(width = 0.5)
  ) +
  geom_errorbar(
    data = contrast_df_fruit_cool,
    aes(
      x = region,
      ymin = lower.CL,
      ymax = upper.CL,
      color = treat_competition
    ),
    linewidth = 1,
    width = 0.1,
    position = position_dodge(width = 0.5)
  ) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  # significance labels
  geom_text(
    data = contrast_df_fruit_cool,
    aes(
      x = region,
      y = estimate,
      label = ifelse(p.value < 0.001, "***",
                     ifelse(p.value < 0.01, "**",
                            ifelse(p.value < 0.05, "*", "n.s."))),
      color = treat_competition
    ),
    vjust = -1.5,
    position = position_dodge(width = 0.5),
    show.legend = FALSE,
    size = 10
  ) +
  
  labs(
    x = "Region",
    y = "Δ days shifted fruiting onset (high - low)",
    color = "Biotic interactions",
    shape = "Region",
    title = "Effect of cooling on fruit onset across regions"
  ) +
  
  theme(
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20)
  ) +
  
  scale_color_manual(values = c(
    "with" = "#528B8B",
    "without" = "#CD950C"
  )) +
  
  scale_shape_manual(values = c("Norway" = 16, "Switzerland" = 17))+
  guides(shape = "none")+
  scale_y_continuous(
    limits = c(-30, 90),
    breaks = seq(-30, 90, by = 20))

nor_che_delta_raw_cool_fruit

# 
# ggsave(filename = "Output/Onset/Fruiting_onset_cooling_effect_NOR_CHE.png", 
#        plot = nor_che_delta_raw_cool_fruit, width = 12, height = 8, units = "in")






















































