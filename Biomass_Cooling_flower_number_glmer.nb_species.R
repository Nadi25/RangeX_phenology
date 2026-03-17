
# BIOMASS 7 ---------------------------------------------------------------

# Effect of cooling on flower number glmer.nb adjusted for biomass -------------------------------------
# negative binomial model is better for our data

#  final figure of effect of cooling on flower number: p_bio_nb_species ---------------------

# the best model option for effect of cooling: m_flower_number_cool_bio_nb_species


# uses individual species models for biomass prediction 

# Adding predicted biomass to phenology and fit species models with pred biomass ------------

## Data used: 
## Date:      03.03.26
## Author:    Nadine Arzt
## Purpose:   Does biomass mitigate the effect of cooling on no of flowers?

# load library ------------------------------------------------------------
library(conflicted)
conflict_prefer_all("dplyr", quiet = TRUE)
library(tidyverse)

library(lme4)
library(ggeffects)
library(emmeans)
library(performance)

# script with species models to predict biomass
source("Biomass_traits_correlation_per_species.R")

# script that combines pred biomass with phenology
source("Biomass_phenology_combine_species_models.R")

# add predicted species biomass to pheno data -------------------------------------
bio_flower_species_unique <- bio_flower_species |>
  group_by(unique_plant_ID) |>
  summarise(pred_log_biomass_species = first(pred_log_biomass_species),
            .groups = "drop")


phenology_bio_species <- phenology |>
  left_join(bio_flower_species_unique,
            by = "unique_plant_ID")


# change reference to be low site ----------------------------------------------
# by factor
phenology_bio_species <- phenology_bio_species |>
  mutate(
    site = factor(site),
    site = relevel(site, ref = "low")
  )
# now low site is the reference and will come first in the plot


# center the biomass ------------------------------------------------------
# there are two clusters in biomass
# this centering uses whether a plant is above or below the species average biomass
# this would help when larger individuals within a species produce more flowers
# with this the errorbars in the final plot are shorter
# BUT is that needed?
phenology_bio_species <- phenology_bio_species |>
  group_by(species) |>
  mutate(pred_bio_centered =
           pred_log_biomass_species - mean(pred_log_biomass_species, na.rm = TRUE)) |>
  ungroup()


# fit the glmer.nb with species biomass ------------------------------------
# m_flower_number_cool_bio_nb_species <- glmer.nb(
#   value ~ site * treat_competition + pred_log_biomass_species +
#     (1|species) + (1|block_ID) , 
#   data = phenology_bio_species,
#   subset = treat_warming == "ambient",
#   control = glmerControl(optimizer = "bobyqa")
# )
# summary(m_flower_number_cool_bio_nb_species)

# tried with (1|unique_plant_ID) but doesnt improve

m_flower_number_cool_bio_nb_species <- glmer.nb(
  value ~ site * treat_competition + pred_bio_centered +
    (1|species) + (1|block_ID),
  data = phenology_bio_species,
  subset = treat_warming == "ambient",
  control = glmerControl(optimizer = "bobyqa")
)
summary(m_flower_number_cool_bio_nb_species)

# predict number of flowers ----------------------------------------------------------------------
emm_fl_num_cool_bio_nb_species <- emmeans(
  m_flower_number_cool_bio_nb_species,
  ~ site | treat_competition,
  cov.reduce = mean, #
  type = "response"
)


emm_fl_num_df_cool_bio_nb_species <- as.data.frame(
  summary(emm_fl_num_cool_bio_nb_species, infer = TRUE)
) |>
  rename(
    emmean = response,
    lower.CL = asymp.LCL,
    upper.CL = asymp.UCL
  )


# significances -----------------------------------------------------------
contr_site_bio_nb_species <- contrast(emm_fl_num_cool_bio_nb_species, method = "pairwise") |>
  as.data.frame() |>
  mutate(type = "site",
         sig = case_when(
           p.value < 0.001 ~ "***",
           p.value < 0.01  ~ "**",
           p.value < 0.05  ~ "*",
           TRUE ~ "ns"
         ))

contr_comp_bio_nb_species <- contrast(
  emmeans(m_flower_number_cool_bio_nb_species, ~ treat_competition | site),
  method = "pairwise"
) |>
  as.data.frame() |>
  mutate(type = "competition",
         sig = case_when(
           p.value < 0.001 ~ "***",
           p.value < 0.01  ~ "**",
           p.value < 0.05  ~ "*",
           TRUE ~ "ns"
         ))



# prepare plot, sign, raw data... -----------------------------------------
dodge  <- position_dodge(width = 0.3)
offset <- 0.3 / 2

y_max <- max(emm_fl_num_df_cool_bio_nb_species$upper.CL, na.rm = TRUE)
spacing <- y_max * 0.27


site_levels <- levels(emm_fl_num_df_cool_bio_nb_species$site)

ann_site_bio_nb_species <- contr_site_bio_nb_species |> 
  mutate(
    site1 = site_levels[1],
    site2 = site_levels[2],
    x_shift = ifelse(treat_competition == "with", -offset, +offset),
    xmin = as.numeric(factor(site1, levels = site_levels)) + x_shift,
    xmax = as.numeric(factor(site2, levels = site_levels)) + x_shift,
    y = y_max + row_number() * spacing,
    label = sig
  )

ann_comp_bio_nb_species <- contr_comp_bio_nb_species |> 
  mutate(
    site_x = as.numeric(factor(site, levels = site_levels)),
    xmin = site_x - offset,
    xmax = site_x + offset,
    y = y_max + (nrow(ann_site_bio_nb_species) + row_number()) * spacing,
    label = sig
  )

# raw means for jittered points --------------------------------------
raw_means <- phenology |>
  filter(treat_warming == "ambient") |> 
  group_by(site, treat_competition, species, block_ID) |> 
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")

# helper to draw brackets -------------------------------------------------
bracket_geoms <- function(df) {
  list(
    geom_segment(data = df, aes(x = xmin, xend = xmax, y = y, yend = y)),
    geom_segment(data = df, aes(x = xmin, xend = xmin, y = y - 0.05, yend = y)),
    geom_segment(data = df, aes(x = xmax, xend = xmax, y = y - 0.05, yend = y)),
    geom_text(data = df, aes(x = (xmin + xmax)/2, y = y + 0.05, label = label),
              vjust = 0, size = 5)
  )
}

# final plot ---------------------------------------------------------
p_bio_nb_species <- ggplot() +
  
  # raw data (optional)
  geom_jitter(
    data = raw_means,
    aes(x = site, y = mean_value, color = treat_competition),
    width = 0.1, alpha = 0.3, size = 1.5
  ) +
  
  # EMM points
  geom_point(
    data = emm_fl_num_df_cool_bio_nb_species,
    aes(x = site, y = emmean, color = treat_competition),
    size = 5,
    position = dodge
  ) +
  
  # EMM CI
  geom_errorbar(
    data = emm_fl_num_df_cool_bio_nb_species,
    aes(x = site, ymin = lower.CL, ymax = upper.CL, color = treat_competition),
    width = 0.2,
    linewidth = 0.8,
    position = dodge
  ) +
  
  # EMM lines
  geom_line(
    data = emm_fl_num_df_cool_bio_nb_species,
    aes(x = site, y = emmean, color = treat_competition,
        group = treat_competition),
    position = dodge
  ) +
  
  # competition brackets (within site)
  bracket_geoms(ann_comp_bio_nb_species) +
  
  # site brackets (within competition)
  bracket_geoms(ann_site_bio_nb_species) +
  
  labs(
    x = "Site",
    y = "Predicted mean number of flowers\n(adjusted for biomass)",
    color = "Biotic interactions",
    title = "Effect of site and competition on flower number\nincluding biomass (glmer.nb)\nspecies models"
  ) +
  
  scale_color_manual(values = c("#528B8B", "#CD950C"))
p_bio_nb_species

# ggsave(filename = "Output/Biomass/Cooling_flower_number_NOR_adjusted_biomass(centered)_species_models_glmer.nb.png", 
#        plot = p_bio_nb_species, width = 10, height = 8, units = "in")




