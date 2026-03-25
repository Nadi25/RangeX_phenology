
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




















