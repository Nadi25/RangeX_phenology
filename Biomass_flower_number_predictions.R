

# BIOMASS 7 ---------------------------------------------------------------


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
library(multcomp)
library(multcompView)
library(gt)



# source the phenology data -----------------------------------------------
source("Data_preparation_phenology_NOR_CHE_combined.R")

# use phenology2 which has combined site and temperature treatment

# compare low site ambient with high site ambient = site effect



# script with species models to predict biomass
source("Biomass_traits_correlation_per_species.R")

# script that combines pred biomass with phenology
source("Biomass_phenology_combine_species_models.R")

# set theme for plots ------------------------------------
theme_set(theme_bw())

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

phenology_bio_species$treatment_site_temp <- 
  recode(phenology_bio_species$treatment_site_temp,
         "low_ambient" = "lo_ambi",
         "high_ambient" = "hi_ambi",
         "high_warmed" = "hi_warm")

phenology_bio_species <- phenology_bio_species |>
  mutate(treatment_site_temp= factor(treatment_site_temp,
                                     levels = c("lo_ambi",
                                                "hi_warm",
                                                "hi_ambi")))

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
# interactive
m_flower_number_bio <- glmer.nb(max_flower_number ~ treatment_site_temp * treat_competition *
                                  pred_log_biomass_species + (1|species) + (1|block_ID),
                            data = max_flower_per_plant_bio,
                            control = glmerControl(optimizer = "bobyqa"))

summary(m_flower_number_bio)


car::Anova(m_flower_number_bio)

emmeans(m_flower_number_bio,
        pairwise ~ treatment_site_temp * treat_competition)

biomass_mean <- mean(
  max_flower_per_plant_bio$pred_log_biomass_species,
  na.rm = TRUE)
biomass_mean


emm_bio <- emmeans(
  m_flower_number_bio,
  pairwise ~ treatment_site_temp * treat_competition,
  at = list(
    pred_log_biomass_species = biomass_mean),
  type = "response")
emm_bio


small_bio <- quantile(
  max_flower_per_plant_bio$pred_log_biomass_species,
  0.25, na.rm = TRUE)

large_bio <- quantile(
  max_flower_per_plant_bio$pred_log_biomass_species,
  0.75, na.rm = TRUE)


emm_bio_small <- emmeans(
  m_flower_number_bio,
  pairwise ~ treatment_site_temp * treat_competition,
  at = list(pred_log_biomass_species = small_bio),
  type = "response")
emm_bio_small

emm_bio_large <- emmeans(
  m_flower_number_bio,
  pairwise ~ treatment_site_temp * treat_competition,
  at = list(pred_log_biomass_species = large_bio),
  type = "response")
emm_bio_large


























#######################################################

# additive
m_flower_number_bio2 <- glmer.nb(max_flower_number ~ treatment_site_temp * treat_competition + 
                                  pred_log_biomass_species + (1|species) + (1|block_ID),
                                data = max_flower_per_plant_bio,
                                control = glmerControl(optimizer = "bobyqa"))

summary(m_flower_number_bio2)

car::Anova(m_flower_number_bio2)

AIC(m_flower_number_bio, m_flower_number_bio2)


summary(glht(m_flower_number_bio2,
             linfct = mcp(treatment_site_temp = "Tukey")))


emmeans(m_flower_number_bio2,
        pairwise ~ treatment_site_temp * treat_competition)


performance::model_performance(m_flower_number_bio2)





make_predictions_bio <- function(model, data) {
  
  biomass_mean <- mean(data$pred_log_biomass_species, na.rm = TRUE)
  
  newdat <- expand.grid(
    treatment_site_temp = c("lo_ambi", "hi_warm", "hi_ambi"),
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

pred_bio <- make_predictions_bio(m_flower_number_bio, max_flower_per_plant_bio)
pred_bio



pred_bio2 <- make_predictions_bio(m_flower_number_bio2, max_flower_per_plant_bio)
pred_bio2






# rename treat info -------------------------------------------------------
# pred_bio$treatment_site_temp <- 
#   recode(pred_bio$treatment_site_temp,
#          "low_ambient" = "lo_ambi",
#          "high_ambient" = "hi_ambi",
#          "high_warmed" = "hi_warm")
# 
# pred_bio2$treatment_site_temp <- 
#   recode(pred_bio$treatment_site_temp,
#          "low_ambient" = "lo_ambi",
#          "high_ambient" = "hi_ambi",
#          "high_warmed" = "hi_warm")
# 
# max_flower_per_plant_bio$treatment_site_temp <- 
#   recode(max_flower_per_plant_bio$treatment_site_temp,
#          "low_ambient" = "lo_ambi",
#          "high_ambient" = "hi_ambi",
#          "high_warmed" = "hi_warm")
# 


# this should be new script
# plot --------------------------------------------------
pd <- position_dodge(width = 0.9) 

flow_no_bio2 <- ggplot(pred_bio2, aes(
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
      "lo_ambi" = 16,   # circle
      "hi_ambi" = 17,    # triangle
      "hi_warm" = 2    # square 
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
flow_no_bio2

# ggsave(filename = "Output/Biomass/Transplantation_warming_flower_number_predictions_NOR_glmer.nb.png", 
#       plot = flow_no, width = 12, height = 8, units = "in")



# Plot with raw as violin ------------------------------------------------
# additive

flow_no_bio3 <- ggplot(pred_bio2, aes(
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
      "lo_ambi" = 16,
      "hi_ambi" = 17,
      "hi_warm" = 2
    )
  ) +
  
  labs(
    x = "Site temperature treatment",
    y = "Log (max flower number)",
    title = "Effect of transplantation and warming on flower number adjusted for biomass additive",
    shape = "Treatment site × warming",
    color = "Biotic interactions",
    fill = "Biotic interactions"
  ) +
  
  guides(shape = "none")+
  coord_transform(y = "log1p")

flow_no_bio3

# ggsave(filename = "Output/Biomass/Transplantation_warming_flower_number_predictions_Biomass_NOR_glmer.nb_violin_additive.png", 
#       plot = flow_no_bio3, width = 12, height = 8, units = "in")





# plot interactive

flow_no_bio4 <- ggplot(pred_bio, aes(
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
      "lo_ambi" = 16,
      "hi_ambi" = 17,
      "hi_warm" = 2
    )
  ) +
  
  labs(
    x = "Site temperature treatment",
    y = "Log (max flower number)",
    title = "Effect of transplantation and warming on flower number adjusted for biomass interactive",
    shape = "Treatment site × warming",
    color = "Biotic interactions",
    fill = "Biotic interactions") +
  
  guides(shape = "none")+
  coord_transform(y = "log1p")+
  theme(legend.position = "bottom")

flow_no_bio4

# ggsave(filename = "Output/Biomass/Transplantation_warming_flower_number_predictions_Biomass_NOR_glmer.nb_violin_interactive.png", 
#       plot = flow_no_bio4, width = 12, height = 8, units = "in")





# plot with stars as significances ----------------------------------------

emm_bio_contr <- emm_bio$contrasts |> 
  as.data.frame()

emm_bio_contr <- emm_bio_contr |> mutate(stars = case_when(p.value < 0.001 ~ "***",
                                                   p.value < 0.01 ~ "**",
                                                   p.value < 0.05 ~ "*",
                                                   TRUE ~ ""))
emm_bio_contr

emm_bio_contr2 <- emm_bio_contr |> 
  filter(contrast %in% c("lo_ambi with / lo_ambi without",
                         "hi_warm with / hi_warm without",
                         "hi_ambi with / hi_ambi without",
                         "lo_ambi with / hi_ambi with",
                         "lo_ambi without / hi_ambi without",
                         "hi_warm with / hi_ambi with",
                         "hi_warm without / hi_ambi without"))
emm_bio_contr2

pos_df <- tibble(
  group = c(
    "lo_ambi with",
    "lo_ambi without",
    "hi_warm with",
    "hi_warm without",
    "hi_ambi with",
    "hi_ambi without"
  ),
  x = c(
    0.775, 1.225,
    1.775, 2.225,
    2.775, 3.225
  )
)



brackets <- emm_bio_contr2 |>
  separate(
    contrast,
    into = c("group1", "group2"),
    sep = " / "
  )
brackets <- brackets |>
  left_join(pos_df, by = c("group1" = "group")) |>
  rename(xmin = x) |>
  left_join(pos_df, by = c("group2" = "group")) |>
  rename(xmax = x)
brackets


brackets |>
  select(group1, group2, xmin, xmax, stars)

top_y <- max(flower_number_pred$upper)

brackets <- brackets |>
  mutate(
    y = top_y + c(18, 6, 27, 6, 6, 13, 22)
  )
brackets

theme_set(theme_bw(base_size = 20))

flow_no_bio4_sig <- flow_no_bio4 +
  geom_segment(
    data = brackets,
    aes(x = xmin, xend = xmax,
        y = y, yend = y),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = brackets,
    aes(x = xmin, xend = xmin,
        y = y, yend = y - 0.5),
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = brackets,
    aes(x = xmax, xend = xmax,
        y = y, yend = y - 0.5),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = brackets,
    aes(
      x = (xmin + xmax)/2,
      y = y + 0.5,
      label = stars
    ),
    inherit.aes = FALSE)+
  theme(plot.title = element_text(size = 15))
flow_no_bio4_sig

# ggsave(filename = "Output/Biomass/Transplantation_warming_flower_number_predictions_Biomass_NOR_glmer.nb_violin_interactive_sig.png", 
#       plot = flow_no_bio4_sig, width = 12, height = 9, units = "in")


tab <- emm_bio_contr2 |>
  gt() |>
  fmt_number(
    columns = c(ratio, SE),
    decimals = 2
  ) |>
  fmt_number(
    columns = p.value,
    decimals = 3
  ) |>
  cols_label(
    contrast = "Contrast",
    ratio = "Ratio",
    SE = "SE",
    p.value = "p",
    stars = ""
  ) |>
  cols_hide(c(df, null, z.ratio)) |>
  tab_header(
    title = md("**Number of flowers adjusted for biomass comparisons**")
  )
tab

#gtsave(tab, "Output/Biomass/Number_flowers_Biomass_signif_NOR.docx")



t2 <- tbl_regression(
  m_flower_number_bio2,
  exponentiate = TRUE) |>
  as_gt() |>
  tab_header(
    title = md("**Number of flowers adj biomass NOR glm.nb**"))
t2

#gtsave(t2, "Output/Biomass/Number_flowers_Biomass_model_summary_NOR.docx")




