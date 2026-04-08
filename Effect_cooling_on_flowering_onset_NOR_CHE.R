
# Effect of transplantation on flowering onset ---------------------------------------

## Data used: RangeX_clean_Phenology_2022_CHE.csv
##            RangeX_clean_phenology_2023_NOR.csv
##            RangeX_clean_MetadataFocal_CHE.csv
##            RangeX_metadata_focal_NOR.csv
## Date:      02.09.25
## Author:    Nadine Arzt
## Purpose:   Effect of transplantation on flowering onset NOR and CHE
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
phenology_cool <- phenology |> 
  filter(treat_warming == "ambi")


# julian days -------------------------------------------------------------
# che and nor was measured in two years but if we count the days in each year it should be fine
phenology_cool$jday <- yday(phenology_cool$date_measurement)


# calculate flowering onset ------------------------------------------------
flowering_onset_n_c_cool <- phenology_cool |> 
  filter(phenology_stage == "No_FloOpen", value > 0) |>
  group_by(region, site, species, unique_plant_ID, block_ID, treat_competition) |>
  summarise(onset = min(jday, na.rm = TRUE), .groups = "drop") |>
  # remove groups where flowering never occurred
  filter(is.finite(onset))




# model with region for flowering onset lmer ----------------------------------
m_onset_n_c_cooling <- lmerTest::lmer(onset ~ region * site * treat_competition + 
                                        (1|species) + (1|block_ID), 
                                      data = flowering_onset_n_c_cool)

summary(m_onset_n_c_cooling)

model_performance(m_onset_n_c_cooling)
check_collinearity(m_onset_n_c_cooling)
check_model(m_onset_n_c_cooling)

# check residuals
plot(m_onset_n_c_cooling)
qqnorm(residuals(m_onset_n_c_cooling)); qqline(residuals(m_onset_n_c_cooling))
hist(residuals(m_onset_n_c_cooling))


# calculate emmeans -------------------------------------------------------
# get emmeans for warming within each region × competition
emm_n_c_cool <- emmeans(m_onset_n_c_cooling, ~ site | region * treat_competition)
# this calculates marginal means per site
# for each region * competition combination, we have one mean for low and high that we can then compare


# get contrasts high - low ------------------------------------------------
# compute contrasts (high - low) within each competition level
contr_n_c_cool<- contrast(emm_n_c_cool, method = list("hi - lo" = c(1, -1)))


# using summary keeps the p-values
contrast_df_n_c_cool <- as.data.frame(summary(contr_n_c_cool, infer = TRUE))

# this would add the significances between with and without biotic --------
contrast_df_n_c_cool <- contrast_df_n_c_cool |>
  mutate(
    stars = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE ~ "n.s."
    )
  )




# plot with raw data points ----------------------------------------------------
# compute mean onset per treatment × group
onset_means_cool <- flowering_onset_n_c_cool |>
  group_by(region, site, species, block_ID, treat_competition) |>
  summarise(mean_onset = mean(onset, na.rm = TRUE), .groups = "drop")

# pivot to get high vs low in same row
delta_onset_cool <- onset_means_cool |>
  pivot_wider(names_from = site, values_from = mean_onset) |>
  mutate(delta = hi - lo) |>
  filter(!is.na(delta))

# check result
head(delta_onset_cool)


# plot raw deltas + model estimates
nor_che_delta_raw_cool <- ggplot() +
  # raw deltas (jittered for visibility)
  geom_jitter(data = delta_onset_cool,
              aes(x = region, y = delta, color = treat_competition),
              width = 0.1, alpha = 0.4, size = 3) +
  
  # model-based warming effects
  geom_point(data = contrast_df_n_c_cool, 
             aes(x = region, y = estimate, color = treat_competition, shape = region),
             size = 6, position = position_dodge(width = 0.5)) +
  geom_errorbar(data = contrast_df_n_c_cool,
                aes(x = region, ymin = lower.CL, ymax = upper.CL, color = treat_competition),
                linewidth = 1, width = 0.1,
                position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  # significance labels
  geom_text(data = contrast_df_n_c_cool,
            aes(x = region, y = estimate, 
                label = ifelse(p.value < 0.001, "***",
                               ifelse(p.value < 0.01, "**",
                                      ifelse(p.value < 0.05, "*", "n.s."))),
                color = treat_competition),
            vjust = -2, position = position_dodge(width = 0.7),
            show.legend = FALSE, size = 12) +
  
  labs(x = "Region",
       y = "Δ days shifted flowering onset (high - low)",
       title = "Effect of cooling through transplantation on flowering onset across regions",
       color = "Biotic interactions") +
  
  theme(
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20)
  ) +
  
  scale_color_manual(values = c("#528B8B", "#CD950C"))+
  scale_shape_manual(values = c("Norway" = 16, "Switzerland" = 17))+
  guides(shape = "none")+
  scale_y_continuous(
    limits = c(-30, 90),
    breaks = seq(-30, 90, by = 20))

nor_che_delta_raw_cool

# 
# ggsave(filename = "Output/Onset/Flowering_onset_cooling_effect_NOR_CHE.png", 
#        plot = nor_che_delta_raw_cool, width = 12, height = 8, units = "in")





# Plot with interactions --------------------------------------------------

emm_full <- emmeans(
  m_onset_n_c_cooling,
  ~ site * treat_competition | region)


contr_cooling <- contrast(
  emm_full,
  method = list("hi - lo" = c(1, -1)),
  by = c("region", "treat_competition"),
  simple = "site"
)


df_cooling <- as.data.frame(summary(contr_cooling, infer = TRUE)) |>
  mutate(effect = "Transplantation",
         group = treat_competition)


contr_comp <- contrast(
  emm_full,
  method = list("with - without" = c(1, -1)),
  by = c("region", "site"),
  simple = "treat_competition"
)


df_comp <- as.data.frame(summary(contr_comp, infer = TRUE)) |>
  mutate(effect = "Biotic interactions",
         group = site)


levels(flowering_onset_n_c_cool$site)
levels(flowering_onset_n_c_cool$treat_competition)


contr_interaction <- contrast(
  emm_full,
  interaction = "pairwise"
)

df_interaction <- as.data.frame(summary(contr_interaction, infer = TRUE)) |>
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

df_cooling     <- add_stars(df_cooling)
df_comp        <- add_stars(df_comp)
df_interaction <- add_stars(df_interaction)

plot_df <- bind_rows(
  df_cooling     |> select(region, effect, group, estimate, lower.CL, upper.CL, stars),
  df_comp        |> select(region, effect, group, estimate, lower.CL, upper.CL, stars),
  df_interaction |> select(region, effect, group, estimate, lower.CL, upper.CL, stars)
)

plot_df <- plot_df |>
  mutate(effect = factor(effect,
                         levels = c("Transplantation",
                                    "Biotic interactions",
                                    "Interaction")))
#
delta_raw <- flowering_onset_n_c_cool |>
  pivot_wider(names_from = c(site, treat_competition),
              values_from = onset) |>
  mutate(
    cooling_with    = hi_with - lo_with,
    cooling_without = hi_without - lo_without,
    interaction     = cooling_with - cooling_without
  )


cooling_flower_interaction <- ggplot(plot_df,
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
cooling_flower_interaction


# ggsave(filename = "Output/Onset/Flowering_onset_cooling_effect_NOR_CHE_interaction.png", 
#        plot = cooling_flower_interaction, width = 12, height = 15, units = "in")


















