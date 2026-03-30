

# Effect of transplantation on Budding onset with and without competition ---------------------------------------

## Data used: RangeX_clean_Phenology_2022_CHE.csv
##            RangeX_clean_phenology_2023_NOR.csv
##            RangeX_clean_MetadataFocal_CHE.csv
##            RangeX_metadata_focal_NOR.csv
## Date:      02.03.26
## Author:    Nadine Arzt
## Purpose:   Effect of transplantation on Budding onset NOR and CHE


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



# calculate budding onset ------------------------------------------------
onset_bud_cool <- phenology_cool |> 
  filter(phenology_stage == "No_Buds", value > 0) |>
  group_by(region, site, treat_competition, species, block_ID, unique_plot_ID, 
           unique_plant_ID, phenology_stage) |>
  summarise(onset = min(jday, na.rm = TRUE), .groups = "drop") |>
  # remove groups where budding never occurred
  filter(is.finite(onset))

# model with region for budding onset lmer ----------------------------------
m_onset_bud_cooling <- lmerTest::lmer(onset ~ region * site * treat_competition + 
                                        (1|species) + (1|block_ID), 
                                      data = onset_bud_cool)

summary(m_onset_bud_cooling)

model_performance(m_onset_bud_cooling)
check_collinearity(m_onset_bud_cooling)
check_model(m_onset_bud_cooling)

# check residuals
plot(m_onset_bud_cooling)
qqnorm(residuals(m_onset_bud_cooling)); qqline(residuals(m_onset_bud_cooling))
hist(residuals(m_onset_bud_cooling))

isSingular(m_onset_bud_cooling)


# emmeans -----------------------------------------------------------------
# get emmeans for warming within each region × competition
# this calculates marginal means per site
emm_bud_cool <- emmeans(m_onset_bud_cooling, ~ site | region * treat_competition)



# contrasts high-low ------------------------------------------------------
# compute contrasts (high - low) within each competition level

contr_bud_cool <- contrast(
  emm_bud_cool,
  method = list("hi - lo" = c(1, -1)),
  by = c("region", "treat_competition")
)



# using summary keeps the p-values
contrast_df_bud_cool <- as.data.frame(summary(contr_bud_cool, infer = TRUE))



# plot with raw data points ----------------------------------------------------
# compute mean onset per treatment × group
onset_bud_means_cool <- onset_bud_cool |>
  group_by(region, site, treat_competition, species, block_ID) |>
  summarise(mean_onset = mean(onset, na.rm = TRUE), .groups = "drop")

# pivot to get ambi vs warm in same row
delta_onset_bud_cool <- onset_bud_means_cool |>
  pivot_wider(names_from = site, values_from = mean_onset) |>
  mutate(delta = hi - lo) |>
  filter(!is.na(delta))

# check result
head(delta_onset_bud_cool)



# plot raw deltas + model estimates ---------------------------------------

nor_che_delta_raw_cool_bud <- ggplot() +
  # raw deltas
  geom_jitter(
    data = delta_onset_bud_cool,
    aes(
      x = region,
      y = delta,
      color = treat_competition
    ),
    width = 0.1, alpha = 0.4, size = 3
  ) +
  
  # model-based contrasts
  geom_point(
    data = contrast_df_bud_cool,
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
    data = contrast_df_bud_cool,
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
    data = contrast_df_bud_cool,
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
    y = "Δ days shifted budding onset (high - low)",
    color = "Biotic interactions",
    shape = "Region",
    title = "Effect of cooling on budding onset across regions"
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

nor_che_delta_raw_cool_bud

# 
# ggsave(filename = "Output/Onset/Budding_onset_cooling_effect_NOR_CHE.png", 
#        plot = nor_che_delta_raw_cool_bud, width = 12, height = 8, units = "in")



# Plot with interactions --------------------------------------------------

emm_full_bud <- emmeans(
  m_onset_bud_cooling,
  ~ site * treat_competition | region)


contr_cooling_bud <- contrast(
  emm_full_bud,
  method = list("hi - lo" = c(1, -1)),
  by = c("region", "treat_competition"),
  simple = "site"
)


df_cooling_bud <- as.data.frame(summary(contr_cooling_bud, infer = TRUE)) |>
  mutate(effect = "Transplantation",
         group = treat_competition)


contr_comp_bud <- contrast(
  emm_full_bud,
  method = list("with - without" = c(1, -1)),
  by = c("region", "site"),
  simple = "treat_competition"
)


df_comp_bud <- as.data.frame(summary(contr_comp_bud, infer = TRUE)) |>
  mutate(effect = "Biotic interactions",
         group = site)



contr_interaction_bud <- contrast(
  emm_full_bud,
  interaction = "pairwise"
)

df_interaction_bud <- as.data.frame(summary(contr_interaction_bud, infer = TRUE)) |>
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

df_cooling_bud     <- add_stars(df_cooling_bud)
df_comp_bud        <- add_stars(df_comp_bud)
df_interaction_bud <- add_stars(df_interaction_bud)

plot_df_bud <- bind_rows(
  df_cooling_bud     |> select(region, effect, group, estimate, lower.CL, upper.CL, stars),
  df_comp_bud        |> select(region, effect, group, estimate, lower.CL, upper.CL, stars),
  df_interaction_bud |> select(region, effect, group, estimate, lower.CL, upper.CL, stars)
)

plot_df_bud <- plot_df_bud |>
  mutate(effect = factor(effect,
                         levels = c("Transplantation",
                                    "Biotic interactions",
                                    "Interaction")))

# delta_raw_bud <- m_onset_bud_cooling |>
#   pivot_wider(names_from = c(site, treat_competition),
#               values_from = onset) |>
#   mutate(
#     cooling_with    = hi_with - lo_with,
#     cooling_without = hi_without - lo_without,
#     interaction     = cooling_with - cooling_without
#   )


cooling_budding_interaction <- ggplot(plot_df_bud,
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
       y = "Effect size (days)",
       color = "Condition")
cooling_budding_interaction


# ggsave(filename = "Output/Onset/Budding_onset_cooling_effect_NOR_CHE_interaction.png", 
#        plot = cooling_budding_interaction, width = 12, height = 15, units = "in")



















































