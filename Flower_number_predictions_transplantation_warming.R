


# Flower number predictions -----------------------------------------------

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



# factor treatment --------------------------------------------------------
phenology2$treat_competition <- factor(phenology2$treat_competition)


# filter only Norway ------------------------------------------------------
phenology_flower_count_nor <- phenology2 |>
  filter(region == "Norway")


# exclude cennig and sildio -----------------------------------------------
# cennig is a calculation
# sildio is not correct counts because we also counted stems
# 
phenology_flower_count_nor <- phenology_flower_count_nor |> 
  filter(!species %in% c("cennig", "sildio"))



# filter flowers ----------------------------------------------------------
phenology_flower_count_nor <- phenology_flower_count_nor |>
  filter(phenology_stage == "No_FloOpen")


# get max number of flowers -----------------------------------------------
max_flower_per_plant <- phenology_flower_count_nor |>
  group_by(site, treatment_site_temp, species, treat_competition, block_ID, unique_plant_ID) |>
  summarise(
    max_flower_number = max(value, na.rm = TRUE),
    .groups = "drop"
  )


# fit the model negative binomial  ------------------------------------
m_flower_number <- glmer.nb(max_flower_number ~ treatment_site_temp * treat_competition + 
                              (1|species) + (1|block_ID),
                            data = max_flower_per_plant,
                            control = glmerControl(optimizer = "bobyqa"))

summary(m_flower_number)




performance::model_performance(m_flower_number)

pred <- ggpredict(m_flower_number, c("treatment_site_temp", "treat_competition"))
plot(pred)


make_predictions <- function(model) {
  
  newdat <- expand.grid(
    treatment_site_temp = c("lo_ambi", "hi_ambi", "hi_warm"),
    treat_competition = c("with", "without")) |> 
    as_tibble()
  
  pred <- predict(
    model,
    newdat = newdat,
    re.form = NA,
    type = "response", # not log scale
    se.fit = TRUE) |> 
    
    as_tibble() |>
    mutate(upper = fit + 1.96 * se.fit,
           lower = fit - 1.96 * se.fit) |>
    bind_cols(newdat)
  
  return(pred)
}



flower_number_pred <- make_predictions(m_flower_number)
flower_number_pred


# rename treat info -------------------------------------------------------
flower_number_pred$treatment_site_temp <- 
  recode(flower_number_pred$treatment_site_temp,
  "lo_ambi" = "low ambient",
  "hi_ambi" = "high ambient",
  "hi_warm" = "high warm")

max_flower_per_plant$treatment_site_temp <- 
  recode(max_flower_per_plant$treatment_site_temp,
  "lo_ambi" = "low ambient",
  "hi_ambi" = "high ambient",
  "hi_warm" = "high warm")



# this should be new script
# plot --------------------------------------------------
pd <- position_dodge(width = 0.9) 

flow_no <- ggplot(flower_number_pred, aes(
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
    data = max_flower_per_plant,
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
flow_no

# ggsave(filename = "Output/Biomass/Transplantation_warming_flower_number_predictions_NOR_glmer.nb.png", 
#       plot = flow_no, width = 12, height = 8, units = "in")


flow_no +
  geom_line(
    data = flower_number_pred,
    aes(x = treatment_site_temp, y = fit, color = treat_competition,
        group = treat_competition),
        position = pd) +
  geom_violin(trim = FALSE)







# Plot with raw as violin ------------------------------------------------

flow_no2 <- ggplot(flower_number_pred, aes(
  x = treatment_site_temp,
  y = fit,
  color = treat_competition,
  shape = treatment_site_temp
)) +
  
  # model estimates
  geom_point(position = pd, size = 4, stroke = 1.2) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = .2,
    position = pd
  ) +
  
  # distribution of raw data
  geom_violin(
    data = max_flower_per_plant,
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
  
  geom_line(
    data = flower_number_pred,
    aes(x = treatment_site_temp, y = fit, color = treat_competition,
        group = treat_competition),
    position = pd)+
  
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
    y = "(Predicted) max flower number",
    title = "Effect of transplantation and warming on flower number",
    shape = "Treatment site × warming",
    color = "Biotic interactions",
    fill = "Biotic interactions"
  ) +
  
  guides(shape = "none")

flow_no2

# ggsave(filename = "Output/Biomass/Transplantation_warming_flower_number_predictions_NOR_glmer.nb_violin.png", 
#       plot = flow_no2, width = 12, height = 8, units = "in")


















