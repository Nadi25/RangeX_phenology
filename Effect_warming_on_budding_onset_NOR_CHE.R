
# Effect of warming on budding onset ---------------------------------------

## Data used: RangeX_clean_Phenology_2022_CHE.csv
##            RangeX_clean_phenology_2023_NOR.csv
##            RangeX_clean_MetadataFocal_CHE.csv
##            RangeX_metadata_focal_NOR.csv
## Date:      02.09.25
## Author:    Nadine Arzt
## Purpose:   Effect of warming on budding onset NOR and CHE
##            with and without biotic interactions


# load library ------------------------------------------------------------
library(lme4)
library(ggeffects)
library(broom.mixed)
library(emmeans)
library(lubridate)
library(performance)
library(see)

source("Data_preparation_phenology_NOR_CHE_combined.R")

# set theme for plots for presentation ------------------------------------
theme_set(theme_bw(base_size = 20))



# use this data set
names(phenology)




# NOR and CHE together ---------------------------------------------------

# filter only hi ambi and lo -----------------------------------------------
phenology_high <- phenology |> 
  filter(site == "hi")


# julian days -------------------------------------------------------------
# che and nor was measured in two years but if we count the days in each year it should be fine
phenology_high$jday <- yday(phenology_high$date_measurement)


# calculate budding onset ------------------------------------------------
budding_onset_high <- phenology_high |> 
  filter(phenology_stage == "No_Buds", value > 0) |>
  group_by(region, species, unique_plant_ID, block_ID, treat_warming, treat_competition) |>
  summarise(onset = min(jday, na.rm = TRUE), .groups = "drop") |>
  # remove groups where flowering never occurred
  filter(is.finite(onset))


budding_onset_high$treat_warming <- factor(
  budding_onset_high$treat_warming,
  levels = c("ambi", "warm"))


# model with region for budding onset lmer ----------------------------------
m_onset_budding_warming <- lmerTest::lmer(onset ~ region * treat_warming * treat_competition + 
                                           (1|species) + (1|block_ID), 
                                         data = budding_onset_high)

summary(m_onset_budding_warming)

model_performance(m_onset_budding_warming)
check_collinearity(m_onset_budding_warming)
check_model(m_onset_budding_warming)

# check residuals
plot(m_onset_budding_warming)
qqnorm(residuals(m_onset_budding_warming)); qqline(residuals(m_onset_budding_warming))
hist(residuals(m_onset_budding_warming))


# calculate emmeans -------------------------------------------------------
# get emmeans for warming within each region × competition
emm_warm_bud <- emmeans(m_onset_budding_warming, ~ treat_warming | region * treat_competition)
# this calculates marginal means per site
# for each region * competition combination, we have one mean for ambient and warmed that we can then compare


# get contrasts warmed - ambient ------------------------------------------------
# compute contrasts (warmed - ambient ) within each competition level
contr_warm_bud <- contrast(
  emm_warm_bud,
  method = list("warm - ambi" = c(-1, 1))
)

# using summary keeps the p-values
contr_df_warm_bud <- as.data.frame(summary(contr_warm_bud, infer = TRUE))

# this would add the significances between with and without biotic --------
contr_df_warm_bud <- contr_df_warm_bud |>
  mutate(
    stars = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE ~ "n.s."
    )
  )
contr_df_warm_bud



# plot with raw data points ----------------------------------------------------
# compute mean onset per treatment × group
onset_means_warm_bud <- budding_onset_high |>
  group_by(region, treat_warming, species, block_ID, treat_competition) |>
  summarise(mean_onset = mean(onset, na.rm = TRUE), .groups = "drop")

onset_means_warm_bud$treat_warming <- factor(
  onset_means_warm_bud$treat_warming,
  levels = c("ambi", "warm")
)


# pivot to get high vs low in same row
delta_onset_warm_bud <- onset_means_warm_bud |>
  pivot_wider(names_from = treat_warming, values_from = mean_onset) |>
  mutate(delta = warm - ambi) |>
  filter(!is.na(delta))

# check result
head(delta_onset_warm_bud)


# plot raw deltas + model estimates
nor_che_delta_raw_warm_bud <- ggplot() +
  # raw deltas (jittered for visibility)
  geom_jitter(data = delta_onset_warm_bud,
              aes(x = region, y = delta, color = treat_competition),
              width = 0.1, alpha = 0.4, size = 3) +
  
  # model-based warming effects
  geom_point(data = contr_df_warm_bud, 
             aes(x = region, y = estimate, color = treat_competition, shape = region),
             size = 6, position = position_dodge(width = 0.5)) +
  geom_errorbar(data = contr_df_warm_bud,
                aes(x = region, ymin = lower.CL, ymax = upper.CL, color = treat_competition),
                linewidth = 1, width = 0.1,
                position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  # significance labels
  geom_text(data = contr_df_warm_bud,
            aes(x = region, y = estimate, 
                label = ifelse(p.value < 0.001, "***",
                               ifelse(p.value < 0.01, "**",
                                      ifelse(p.value < 0.05, "*", "n.s."))),
                color = treat_competition),
            vjust = -2, position = position_dodge(width = 0.7),
            show.legend = FALSE, size = 12) +
  
  labs(x = "Region",
       y = "Δ days shifted budding onset (warmed - ambient)",
       title = "Effect of warming  on budding onset across regions",
       color = "Biotic interactions") +
  
  theme(
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20)
  ) +
  
  scale_color_manual(values = c("#528B8B", "#CD950C"))+
  scale_shape_manual(values = c("Norway" = 16, "Switzerland" = 17))+
  guides(shape = "none")+
  scale_y_continuous(
    limits = c(-70, 50),
    breaks = seq(-70, 50, by = 20))

nor_che_delta_raw_warm_bud

# 
# ggsave(filename = "Output/Onset/Budding_onset_warming_effect_NOR_CHE.png", 
#        plot = nor_che_delta_raw_warm_bud, width = 12, height = 8, units = "in")




# Plot with interactions --------------------------------------------------

emm_full_warm_bud <- emmeans(
  m_onset_budding_warming,
  ~ treat_warming * treat_competition | region)


contr_warm_bud <- contrast(
  emm_full_warm_bud,
  method = list("warm - ambi" = c(-1, 1)),
  by = c("region", "treat_competition"),
  simple = "treat_warming"
)


df_warm_bud <- as.data.frame(summary(contr_warm_bud, infer = TRUE)) |>
  mutate(effect = "Warming",
         group = treat_competition)


contr_comp_warm_bud <- contrast(
  emm_full_warm_bud,
  method = list("with - without" = c(1, -1)),
  by = c("region", "treat_warming"),
  simple = "treat_competition"
)


df_comp_warm_bud <- as.data.frame(summary(contr_comp_warm_bud, infer = TRUE)) |>
  mutate(effect = "Biotic interactions",
         group = treat_warming)



contr_interaction_warm_bud <- contrast(
  emm_full_warm_bud,
  interaction = "pairwise"
)

df_interaction_warm_bud <- as.data.frame(summary(contr_interaction_warm_bud, infer = TRUE)) |>
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

df_warm_bud    <- add_stars(df_warm_bud)
df_comp_warm_bud       <- add_stars(df_comp_warm_bud)
df_interaction_warm_bud <- add_stars(df_interaction_warm_bud)

plot_df_warm_bud <- bind_rows(
  df_warm_bud     |> select(region, effect, group, estimate, lower.CL, upper.CL, stars),
  df_comp_warm_bud        |> select(region, effect, group, estimate, lower.CL, upper.CL, stars),
  df_interaction_warm_bud |> select(region, effect, group, estimate, lower.CL, upper.CL, stars)
)

plot_df_warm_bud <- plot_df_warm_bud |>
  mutate(effect = factor(effect,
                         levels = c("Warming",
                                    "Biotic interactions",
                                    "Interaction")))



warming_budding_interaction <- ggplot(plot_df_warm_bud,
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
       y = "Shifted onset (days)",
       color = "Condition")
warming_budding_interaction

warming_budding_interaction <- warming_budding_interaction +
  coord_flip()
warming_budding_interaction


# ggsave(filename = "Output/Onset/Budding_onset_warming_effect_NOR_CHE_interaction.png", 
#        plot = warming_budding_interaction, width = 12, height = 15, units = "in")

















