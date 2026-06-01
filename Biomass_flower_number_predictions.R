


# Flower number predictions with biomass -----------------------------------------------

# Effect of transplantation and warming on flower number NOR glmer.nb ---------------------------------

## Data used: 
## Date:      11.03.2026
## Author:    Nadine Arzt
## Purpose:   Effect of transplantation and warming on flower number NOR with glmer.nb


# comments ----------------------------------------------------------------

# NOR.hi.ambi.vege.wf.09.13.1 has NA in value but why

# CHE: assuming that NA means plant not found
# because there are many 0


# load library ------------------------------------------------------------
library(lme4)
library(ggeffects)
library(broom.mixed)
library(emmeans)
library(lubridate)
library(glmmTMB)

# set theme for plots ------------------------------------
theme_set(theme_bw())


# source the phenology data -----------------------------------------------
source("Data_preparation_phenology_NOR_CHE_combined.R")

# use phenology2 which has combined site and temperature treatment

# compare low site ambient with high site ambient = site effect



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


# add site_warming treatment ----------------------------------------------
phenology_bio_species$treatment_site_temp <- paste(phenology_bio_species$site, phenology_bio_species$treat_warming, sep = "_")


phenology_bio_species <- phenology_bio_species |>
  mutate(treatment_site_temp= factor(treatment_site_temp,
                                     levels = c("low_ambient",
                                                "high_ambient",
                                                "high_warmed")))

# exclude cennig and sildio -----------------------------------------------
# cennig is a calculation
# sildio is not correct counts because we also counted stems
# 
phenology_bio_species2 <- phenology_bio_species |> 
  filter(!species %in% c("cennig", "sildio"))



# filter flowers ----------------------------------------------------------
phenology_bio_species2 <- phenology_bio_species2 |>
  filter(phenology_stage == "No_FloOpen")


# get max number of flowers -----------------------------------------------
max_flower_per_plant_bio <- phenology_bio_species2 |>
  group_by(site, treatment_site_temp, species, treat_competition, block_ID, 
           unique_plant_ID, pred_log_biomass_species) |>
  summarise(
    max_flower_number = max(value, na.rm = TRUE),
    .groups = "drop"
  )


# fit the model negative binomial  ------------------------------------
# m_flower_number_bio <- glmer.nb(max_flower_number ~ treatment_site_temp * treat_competition + 
#                                   pred_log_biomass_species + (1|species) + (1|block_ID),
#                             data = max_flower_per_plant_bio,
#                             control = glmerControl(optimizer = "bobyqa"))
# 
# summary(m_flower_number_bio)




m_flower_number_bio2 <- glmer.nb(max_flower_number ~ treatment_site_temp * treat_competition * 
                                  pred_log_biomass_species + (1|species) + (1|block_ID),
                                data = max_flower_per_plant_bio,
                                control = glmerControl(optimizer = "bobyqa"))

summary(m_flower_number_bio2)

#AIC(m_flower_number_bio, m_flower_number_bio2)







performance::model_performance(m_flower_number_bio2)





make_predictions_bio <- function(model, data) {
  
  biomass_mean <- mean(data$pred_log_biomass_species, na.rm = TRUE)
  
  newdat <- expand.grid(
    treatment_site_temp = c("low_ambient", "high_ambient", "high_warmed"),
    treat_competition = c("with", "without"),
    pred_log_biomass_species = biomass_mean
  )|>
    as_tibble()
  
  pred <- predict(
    model,
    newdata = newdat,
    re.form = NA,
    type = "link", # log scale
    se.fit = TRUE
  )
  
  newdat$fit <- pred$fit
  newdat$se.fit <- pred$se.fit
  
  newdat |>
    mutate(
      upper = exp(fit + 1.96 * se.fit),
      lower = exp(fit - 1.96 * se.fit),
      fit = exp(fit)
    )
}

pred_bio <- make_predictions_bio(m_flower_number_bio2, max_flower_per_plant_bio)
pred_bio









# rename treat info -------------------------------------------------------
pred_bio$treatment_site_temp <- 
  recode(pred_bio$treatment_site_temp,
         "low_ambient" = "low ambient",
         "high_ambient" = "high ambient",
         "high_warmed" = "high warm")

max_flower_per_plant_bio$treatment_site_temp <- 
  recode(max_flower_per_plant_bio$treatment_site_temp,
         "low_ambient" = "low ambient",
         "high_ambient" = "high ambient",
         "high_warmed" = "high warm")



# this should be new script
# plot --------------------------------------------------
pd <- position_dodge(width = 0.9) 

flow_no_bio <- ggplot(pred_bio, aes(
  x = treatment_site_temp,
  y = fit,
  color = treat_competition,              
  shape = treatment_site_temp  
)) +
  
  geom_point(position = pd, size = 4, stroke = 1.2) +
  
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = .2,
                position = pd) +
  
  geom_jitter(
    data = max_flower_per_plant_bio,
    aes(
      x = treatment_site_temp,
      y = max_flower_number,
      color = treat_competition
    ),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.9
    ),
    alpha = 0.05,
    size = 2,
    inherit.aes = FALSE
  )+
  
  scale_color_manual(values = c("#528B8B", "#CD950C")) +   
  
  scale_shape_manual(
    values = c(
      "low ambient" = 16,   # circle
      "high ambient" = 17,    # triangle
      "high warm" = 2    # square 
    )
  ) +
  
  labs(
    x = "Site temperature treatment",
    y = "(Predicted) max flower number",
    title = "Effect of transplantation and warming on flower number",
    shape = "Treatment site × warming",
    color = "Biotic interactions"
  )+
  guides(shape = "none")
flow_no_bio

# ggsave(filename = "Output/Biomass/Transplantation_warming_flower_number_predictions_NOR_glmer.nb.png", 
#       plot = flow_no, width = 12, height = 8, units = "in")



# Plot with raw as violin ------------------------------------------------
theme_set(theme_bw())

flow_no_bio2 <- ggplot(pred_bio, aes(
  x = treatment_site_temp,
  y = fit,
  color = treat_competition,
  shape = treatment_site_temp
)) +
  
  # distribution of raw data
  geom_violin(
    data = max_flower_per_plant_bio,
    aes(
      x = treatment_site_temp,
      y = max_flower_number,
      fill = treat_competition
    ),
    position = position_dodge(width = 0.9),
    alpha = 0.25,
    color = NA,
    trim = TRUE # cuts of at 0 because we dont have negative flower no
  ) +
  
  # model estimates
  geom_point(position = pd, size = 4, stroke = 1.2) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = .2,
    position = pd
  ) +
  
  scale_color_manual(values = c("#528B8B", "#CD950C")) +
  scale_fill_manual(values = c("#528B8B", "#CD950C")) +
  
  scale_shape_manual(
    values = c(
      "low ambient" = 16,
      "high ambient" = 17,
      "high warm" = 2
    )
  ) +
  
  labs(
    x = "Site temperature treatment",
    y = "Log (max flower number)",
    title = "Effect of transplantation and warming on flower number adjusted for biomass",
    shape = "Treatment site × warming",
    color = "Biotic interactions",
    fill = "Biotic interactions"
  ) +
  
  guides(shape = "none")+
  coord_transform(y = "log1p")

flow_no_bio2

# ggsave(filename = "Output/Biomass/Transplantation_warming_flower_number_predictions_Biomass_NOR_glmer.nb_violin.png", 
#       plot = flow_no_bio2, width = 12, height = 8, units = "in")


















